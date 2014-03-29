#include <ctime>
#include <algorithm>
#include <limits>
#include <fstream>
#include <sstream>
#include "fluid.h"
#include "glpointcloud.h"
#include "gltext.h"
#include "settings.h"

// Fluid stuff

template<typename REAL, typename SIZE>
FluidRS<REAL,SIZE>::FluidRS(const PointCloudRS<REAL,SIZE> *pc, FluidParamsPtr params)
  : PointCloudRS<REAL,SIZE>(*pc)
  , m_params(params)
{ }

template<typename REAL, typename SIZE>
FluidRS<REAL,SIZE>::FluidRS(const aiMesh *pc, FluidParamsPtr params)
  : PointCloudRS<REAL,SIZE>(pc)
  , m_params(params)
{ }

template<typename REAL, typename SIZE>
FluidRS<REAL,SIZE>::~FluidRS()
{ }


// the fluid must be initialized before being simulated
template<typename REAL, typename SIZE>
void FluidRS<REAL,SIZE>::init(GLPointCloudRS<REAL, SIZE> *glpc)
{
  m_bmin = m_bbox.corner(Eigen::AlignedBox3f::BottomLeftFloor);
  m_bmax = m_bbox.corner(Eigen::AlignedBox3f::TopRightCeil);

  REAL r = get_radius();
  m_kernel_radius = m_params->kernel_inflation * r;
  m_rest_density = m_params->density;
  m_viscosity = m_params->viscosity;
  m_st = m_params->surface_tension;
  m_mass = m_params->density*8*r*r*r;
  m_recoil_velocity_damping = m_params->recoil_velocity_damping;

  if ( m_params->fluid_type == MCG03 )
    m_c2 = m_params->sound_speed * m_params->sound_speed;
  else if ( m_params->fluid_type == BT07 )
  {
    // a heuristic to determine the speed of sound based on the maximum
    // possible velocity of a particle in the fluid
    float max_velocity2 = 2.0f*M_G*(m_bmax[1] - m_bmin[1]) + m_params->velocity.squaredNorm();
    m_c2 = max_velocity2 / m_params->compressibility;
  }

  m_glpc = glpc;

  m_accel.resizeLike(m_pos);
  m_extern_accel.resizeLike(m_pos);
  m_vel.resizeLike(m_pos);
  m_vel.setZero();
  reset_accel();

  //const clock_t begin_time = clock();
  //float d1 = PointCloudRS<REAL,SIZE>::compute_mindist();
  //float t1 = clock();
  //float d2 = PointCloudRS<REAL,SIZE>::compute_mindist_brute();
  //float t2 = clock();
  //qDebug() << float( t1 - begin_time ) / CLOCKS_PER_SEC << " to get " << d1;
  //qDebug() << float( t2 - t1 ) / CLOCKS_PER_SEC << " to get " << d2;

  // construct output filename (this can be done whenever)
  m_cachefmt = global::dynset.cachedir + "/" + m_params->cacheprefix + "%03d.obj";
}

// clamp value d to min and max boundaries + epsilon,
template<typename REAL, typename SIZE>
inline bool FluidRS<REAL,SIZE>::clamp(REAL &d, REAL min, REAL max)
{
  if ( d < min )
  {
    d = min + 0.005;
    return true;
  }
  else if (d > max)
  {
    d = max - 0.005;
    return true;
  }
  return false;
}

template<typename REAL, typename SIZE>
void FluidRS<REAL,SIZE>::resolve_collisions()
{
  if (this->m_params->fluid_type != MCG03)
    return;

  for (SIZE i = 0; i < this->get_num_vertices(); ++i)
  {
    for (unsigned char j = 0; j < 3; ++j)
    {
      if (clamp(pos_at(i)[j], m_bmin[j], m_bmax[j]))
      {
        vel_at(i)[j] *= -m_recoil_velocity_damping;
      }
    }
  }
}

template<typename REAL, typename SIZE>
inline void FluidRS<REAL,SIZE>::update_data()
{
  if (m_glpc)
    m_glpc->update_data();
}

template<typename REAL, typename SIZE>
inline void FluidRS<REAL,SIZE>::write_to_file(unsigned int frame) 
{
  if (global::dynset.cachedir.empty())
    return;

  // open output file
  FILE * outfile = NULL;
  char buf[64];
  sprintf(buf, m_cachefmt.c_str(), frame);
  outfile = fopen(buf,"w");

  if (!outfile)
  {
    qWarning() << "Output file not opened!";
    return;
  }

  for (SIZE i = 0; i < this->get_num_vertices(); ++i)
    fprintf(outfile, "v %f %f %f\n", pos_at(i)[0], pos_at(i)[1], pos_at(i)[2]);

  for (SIZE i = 1; i <= this->get_num_vertices(); ++i)
    fprintf(outfile, "f %d\n", i);

  fclose(outfile);
}

// return if we successfully read from cached file and have a valid state
template<typename REAL, typename SIZE>
inline bool FluidRS<REAL,SIZE>::read_from_file(unsigned int frame) 
{
  if (global::dynset.cachedir.empty())
    return false;

  // open input file
  std::ifstream infile;
  char buf[64];
  sprintf(buf, m_cachefmt.c_str(), frame);
  infile.open(buf);

  if (!infile.is_open())
  {
    qWarning() << "Input file not opened!";
    return false;
  }

  std::string line;
  SIZE i = 0;
  while ( getline(infile, line) )
  {
    std::istringstream iss(line);
    std::string identifier;
    iss >> identifier;
    if (boost::iequals(identifier, "f"))
      break; // reached faces, done

    if (!boost::iequals(identifier, "v"))
      continue; // skip garbage

    iss >> pos_at(i)[0];
    iss >> pos_at(i)[1];
    iss >> pos_at(i)[2];
    i += 1;
  }

  infile.close();

  if (frame >= global::dynset.frames)
    return true; // last frame, so we're done

  // now we check if there exists the frames cached file, if not we must load
  // a previous cached file and infer the current velocities
  // this allows us to have partially cached data.

  // open input file
  sprintf(buf, m_cachefmt.c_str(), frame+1);
  infile.open(buf);

  if (infile.is_open())
  { // we're good for next frame, just return
    infile.close();
    return true;
  }

  // Oh no! there is no cached file for next frame. Quickly! Figure out what
  // the velocities should be before next simulation step
  // NOTE: this is only an approximation since the actual simulation timestep
  // is often smaller than the time between frames. Alternatively, we may
  // cache velocities separately, but thats expensive
 
  sprintf(buf, m_cachefmt.c_str(), frame-1);
  infile.open(buf);

  if (!infile.is_open())
  { // no previous frame cached, oh well, just leave velocities as is
    // alternatively we may zero out all velocities, not sure what's better
    return true;
  }

  i = 0;
  while ( getline(infile, line) )
  {
    std::istringstream iss(line);
    std::string identifier;
    iss >> identifier;
    if (boost::iequals(identifier, "f"))
      break; // reached faces, done

    if (!boost::iequals(identifier, "v"))
      continue; // skip garbage

    float pos;
    for (unsigned char j = 0; j < 3; ++j)
    {
      iss >> pos;
      vel_at(i)[j] = (pos_at(i)[j] - pos)*global::dynset.fps;
    }
    i += 1;
  }

  infile.close();
  return true;
}


// Typed Fluid Stuff

template<typename REAL, typename SIZE, int FT>
FluidRST<REAL,SIZE,FT>::FluidRST(const PointCloudRS<REAL,SIZE> *pc, FluidParamsPtr params)
  : FluidRS<REAL,SIZE>(pc, params)
{ }

template<typename REAL, typename SIZE, int FT>
FluidRST<REAL,SIZE,FT>::FluidRST(const aiMesh *pc, FluidParamsPtr params)
  : FluidRS<REAL,SIZE>(pc, params)
{ }

template<typename REAL, typename SIZE, int FT>
FluidRST<REAL,SIZE,FT>::~FluidRST()
{ }

template<typename REAL, typename SIZE, int FT>
void FluidRST<REAL,SIZE,FT>::init_processors()
{
  m_fluid_density_proc.init_kernel(m_kernel_radius);
  m_fluid_density_update_proc.init_kernel(m_kernel_radius);
  m_fluid_pressure_proc.init_kernel(m_kernel_radius);
  m_fluid_viscosity_accel_proc.init_kernel(m_kernel_radius);
  m_fluid_pressure_accel_proc.init_kernel(m_kernel_radius);
  m_fluid_surface_tension_accel_proc.init_kernel(m_kernel_radius);

  copy_properties_to_proc(m_fluid_density_proc);
  copy_properties_to_proc(m_fluid_density_update_proc);
  copy_properties_to_proc(m_fluid_pressure_proc);
  copy_properties_to_proc(m_fluid_viscosity_accel_proc);
  copy_properties_to_proc(m_fluid_pressure_accel_proc);
  copy_properties_to_proc(m_fluid_surface_tension_accel_proc);
}

template<typename REAL, typename SIZE, int FT>
template<class OutputType, class KernelType, class ComputeType>
inline void FluidRST<REAL,SIZE,FT>::copy_properties_to_proc(
    CFQ<REAL,SIZE,OutputType,KernelType,ComputeType> &cfq_proc)
{
  cfq_proc.m_mass = m_mass;
  cfq_proc.m_radius = this->get_radius();
  cfq_proc.m_rest_density = m_rest_density;
  cfq_proc.m_viscosity = m_viscosity;
  cfq_proc.m_st = m_st;
  cfq_proc.m_cs2 = m_c2;
  cfq_proc.m_cs = std::sqrt(cfq_proc.m_cs2);
}

template class FluidRS<double, unsigned int>;
template class FluidRST<double, unsigned int, 0>; //MCG03
template class FluidRST<double, unsigned int, 1>; //BT07
template class FluidRST<double, unsigned int, 2>; //AIAST12
template class FluidRST<double, unsigned int, 3>; //DEFAULT
