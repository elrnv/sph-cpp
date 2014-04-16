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
  , m_avg_density(0.0f)
  , m_avg_pressure(0.0f)
  , m_cache(global::dynset.frames+1)
{ }

template<typename REAL, typename SIZE>
FluidRS<REAL,SIZE>::FluidRS(const aiMesh *pc, FluidParamsPtr params)
  : PointCloudRS<REAL,SIZE>(pc)
  , m_params(params)
  , m_avg_density(0.0f)
  , m_avg_pressure(0.0f)
  , m_cache(global::dynset.frames+1)
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

  m_kernel_radius = get_kernel_radius();
  m_rest_density = m_params->density;
  m_viscosity = m_params->viscosity;
  m_st = m_params->surface_tension;
  REAL r = get_radius();
  m_mass = m_params->density*8*r*r*r;
  m_recoil_velocity_damping = m_params->recoil_velocity_damping;

  if ( m_params->fluid_type == MCG03 )
    m_c2 = m_params->sound_speed * m_params->sound_speed;
  else //if ( m_params->fluid_type == BT07 )
  {
    // a heuristic to determine the speed of sound based on the maximum
    // possible velocity of a particle in the fluid
    float max_velocity2 = 2.0f*M_G*(m_bmax[1] - m_bmin[1]) + m_params->velocity.squaredNorm();
    m_c2 = max_velocity2 / m_params->compressibility;
  }

  m_glpc = glpc;

  m_accel.resizeLike(m_pos);
  m_dinv.resize(m_pos.cols());
  m_extern_accel.resizeLike(m_pos);
  m_vel.resizeLike(m_pos);
  for (SIZE i = 0; i < this->get_num_vertices(); ++i)
    m_vel.col(i) = m_params->velocity.template cast<REAL>();
  m_dinv.setConstant(m_rest_density);
  reset_accel();

  //const clock_t begin_time = clock();
  //float d1 = PointCloudRS<REAL,SIZE>::compute_mindist();
  //float t1 = clock();
  //float d2 = PointCloudRS<REAL,SIZE>::compute_mindist_brute();
  //float t2 = clock();
  //qDebug() << float( t1 - begin_time ) / CLOCKS_PER_SEC << " to get " << d1;
  //qDebug() << float( t2 - t1 ) / CLOCKS_PER_SEC << " to get " << d2;

  // construct output filename (this can be done whenever)
  m_savefmt = global::dynset.savedir + "/" + m_params->saveprefix + "%03d.sim";
}

// clamp value d to min and max boundaries + epsilon,
template<typename REAL, typename SIZE>
inline bool FluidRS<REAL,SIZE>::clamp(REAL &d, REAL min, REAL max, REAL tol)
{
  if ( d < min )
  {
    d = min + tol;
    return true;
  }
  else if (d > max)
  {
    d = max - tol;
    return true;
  }
  return false;
}

template<typename REAL, typename SIZE>
inline void FluidRS<REAL,SIZE>::resolve_collisions()
{
  for (SIZE i = 0; i < this->get_num_vertices(); ++i)
  {
    for (unsigned char j = 0; j < 3; ++j)
    {
      if (clamp(pos_at(i)[j], m_bmin[j], m_bmax[j], 0.001))
      {
        vel_at(i)[j] *= -m_recoil_velocity_damping;
      }
    }
  }
}

template<typename REAL, typename SIZE>
inline void FluidRS<REAL,SIZE>::clamp(float tol)
{
  for (SIZE i = 0; i < this->get_num_vertices(); ++i)
    for (unsigned char j = 0; j < 3; ++j)
      clamp(pos_at(i)[j], m_bmin[j], m_bmax[j], tol);
}

template<typename REAL, typename SIZE>
inline void FluidRS<REAL,SIZE>::update_data()
{
  if (m_glpc)
    m_glpc->update_data();
}


// File stuff

template<typename REAL, typename SIZE>
inline void FluidRS<REAL,SIZE>::save(unsigned int frame) 
{
  if (global::dynset.savedir.empty())
    return;

  // open output file
  std::ofstream outfile;
  char buf[128];
  sprintf(buf, m_savefmt.c_str(), frame);
  outfile.open(buf, std::ios::out | std::ios::binary);

  if (!outfile.is_open())
  {
    qWarning() << "Output file" << buf << "not opened!";
    return;
  }

  SIZE num = this->get_num_vertices();
  outfile.write(reinterpret_cast<const char *>(&num), sizeof num);
  float radius = this->get_radius();
  float lev = 0.0f;

  for (SIZE i = 0; i < num; ++i)
  {
    REAL * posptr = m_cache[frame].pos.data() + i*3;
    REAL * velptr = m_cache[frame].vel.data() + i*3;
    float x = posptr[0]; float y = posptr[1]; float z = posptr[2];
    float u = velptr[0]; float v = velptr[1]; float w = velptr[2];
    outfile.write(reinterpret_cast<const char *>(&x), sizeof(float));
    outfile.write(reinterpret_cast<const char *>(&y), sizeof(float));
    outfile.write(reinterpret_cast<const char *>(&z), sizeof(float));
    outfile.write(reinterpret_cast<const char *>(&radius), sizeof(float));
    outfile.write(reinterpret_cast<const char *>(&lev), sizeof(float));
    outfile.write(reinterpret_cast<const char *>(&u), sizeof(float));
    outfile.write(reinterpret_cast<const char *>(&v), sizeof(float));
    outfile.write(reinterpret_cast<const char *>(&w), sizeof(float));
  }

  outfile.close();
}

template<typename REAL, typename SIZE>
inline bool FluidRS<REAL,SIZE>::is_saved(unsigned int frame) 
{
  if (global::dynset.savedir.empty())
    return false;

  // open save file
  FILE * infile = NULL;
  char buf[128];
  sprintf(buf, m_savefmt.c_str(), frame);
  infile = fopen(buf, "r");
  if (!infile)
    return false;

  fclose(infile);
  return true;
}


// return true if loaded, false otherwise
template<typename REAL, typename SIZE>
inline bool FluidRS<REAL,SIZE>::load_saved_cache()
{
  // check if first frame was saved, if so load everything, otherwise ignore
  if (is_saved(0))
  {
    bool read_successful = true;
    for (unsigned int fr = 0; fr <= global::dynset.frames; ++fr)
    {
      read_successful &= read_saved(fr);
      glprintf_trcv(get_color(), "\rLoading cache %d%%", 100*fr/global::dynset.frames);
    }
    glprintf_tr("\r   ");

    if (!read_successful)
      qWarning() << "Some cached files could not be found.";
    return true;
  } 
  return false;
}

// return if we successfully read from saved file and have a valid state
template<typename REAL, typename SIZE>
inline bool FluidRS<REAL,SIZE>::read_saved(unsigned int frame) 
{
  if (global::dynset.savedir.empty())
    return false;

  // open input file
  std::ifstream infile;
  char buf[128];
  sprintf(buf, m_savefmt.c_str(), frame);
  infile.open(buf, std::ios::in | std::ios::binary);

  if (!infile.is_open())
    return false;

  SIZE num;
  infile.read(reinterpret_cast<char *>(&num), sizeof num);

  m_cache[frame].pos.resize(NoChange,num);
  m_cache[frame].vel.resize(NoChange,num);

  for (SIZE i = 0; i < num; ++i)
  {
    float x,y,z,r,l,u,v,w;
    infile.read(reinterpret_cast<char *>(&x), sizeof(float));
    infile.read(reinterpret_cast<char *>(&y), sizeof(float));
    infile.read(reinterpret_cast<char *>(&z), sizeof(float));
    infile.read(reinterpret_cast<char *>(&r), sizeof(float)); // ignore
    infile.read(reinterpret_cast<char *>(&l), sizeof(float)); // ignore
    infile.read(reinterpret_cast<char *>(&u), sizeof(float));
    infile.read(reinterpret_cast<char *>(&v), sizeof(float));
    infile.read(reinterpret_cast<char *>(&w), sizeof(float));
    
    REAL * posptr = m_cache[frame].pos.data() + i*3;
    REAL * velptr = m_cache[frame].vel.data() + i*3;
    posptr[0] = REAL(x); posptr[1] = REAL(y); posptr[2] = REAL(z);
    velptr[0] = REAL(u); velptr[1] = REAL(v); velptr[2] = REAL(w);

    m_cache[frame].valid = true;
  }

  infile.close();

  return true; // last frame, so we're done
}

template<typename REAL, typename SIZE>
inline void FluidRS<REAL,SIZE>::clear_saved() 
{
  if (global::dynset.savedir.empty())
    return;

  for (unsigned int fr = 0; fr <= global::dynset.frames; ++fr)
  {
    char buf[128];
    sprintf(buf, m_savefmt.c_str(), fr);
    remove(buf);
  }
}


// Cache stuff

template<typename REAL, typename SIZE>
inline void FluidRS<REAL,SIZE>::cache(unsigned int frame) 
{
  m_cache[frame].pos = m_pos;
  m_cache[frame].vel = m_vel;
  m_cache[frame].valid = true;

  save(frame);
}

template<typename REAL, typename SIZE>
inline bool FluidRS<REAL,SIZE>::is_cached(unsigned int frame) 
{
  return m_cache[frame].valid;
}

template<typename REAL, typename SIZE>
inline bool FluidRS<REAL,SIZE>::load_cached(unsigned int frame) 
{
  // dont load if not cached or is first frame and next is not cached
  if ( !is_cached(frame) ) 
    return false;

  m_pos = m_cache[frame].pos;
  m_vel = m_cache[frame].vel;

  return true;
}

template<typename REAL, typename SIZE>
inline void FluidRS<REAL,SIZE>::clear_cache() 
{
  bool already_clear = true;
  for (unsigned int fr = 0; fr <= global::dynset.frames; ++fr)
  {
    already_clear &= !m_cache[fr].valid;
    m_cache[fr].valid = false;
  }

  if (already_clear)
    clear_saved();
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
inline void FluidRST<REAL,SIZE,FT>::init_processors()
{
  m_fluid_density_proc.init_kernel(m_kernel_radius);
  m_fluid_density_update_proc.init_kernel(m_kernel_radius);
  m_fluid_accel_proc.init_kernel(m_kernel_radius);
  m_fluid_surface_normal_proc.init_kernel(m_kernel_radius);
  m_fluid_surface_tension_proc.init_kernel(m_kernel_radius);
  m_fluid_prepare_jacobi_proc.init_kernel(m_kernel_radius);
  m_fluid_jacobi_solve1_proc.init_kernel(m_kernel_radius);
  m_fluid_jacobi_solve2_proc.init_kernel(m_kernel_radius);
  m_fluid_pressure_accel_proc.init_kernel(m_kernel_radius);

  m_fluid_density_proc.copy_fluid_params(*this);
  m_fluid_density_update_proc.copy_fluid_params(*this);
  m_fluid_accel_proc.copy_fluid_params(*this);
  m_fluid_surface_normal_proc.copy_fluid_params(*this);
  m_fluid_surface_tension_proc.copy_fluid_params(*this);
  m_fluid_prepare_jacobi_proc.copy_fluid_params(*this);
  m_fluid_jacobi_solve1_proc.copy_fluid_params(*this);
  m_fluid_jacobi_solve2_proc.copy_fluid_params(*this);
  m_fluid_pressure_accel_proc.copy_fluid_params(*this);
}

template class FluidRS<double, unsigned int>;
template class FluidRST<double, unsigned int, 0>; //MCG03
template class FluidRST<double, unsigned int, 1>; //BT07
template class FluidRST<double, unsigned int, 2>; //ICS13
template class FluidRST<double, unsigned int, 3>; //DEFAULT
