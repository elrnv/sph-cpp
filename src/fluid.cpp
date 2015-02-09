#include <ctime>
#include <algorithm>
#include <limits>
#include <fstream>
#include <sstream>
#include "fluid.h"
#include "pointcloud.h"
#include "glpointcloud.h"
#include "gltext.h"
#include "settings.h"
#include "sphgrid.h"

#define M_G 9.81f

// Fluid stuff

Fluid::Fluid(const aiMesh *pc, Index matidx, 
             MaterialManager &matman, FluidParamsPtr params)
  : m_pc(pc, matidx)
  , m_params(params)
  , m_avg_density(0.0f)
  , m_avg_pressure(0.0f)
  , m_color(matman[matidx].kd())
  , m_cache(global::dynset.frames+1)
{ }

Fluid::~Fluid()
{ }

// the fluid must be initialized before being simulated

void 
Fluid::init(const AlignedBox3f &box)
{
  m_bmin = box.corner(AlignedBox3f::BottomLeftFloor);
  m_bmax = box.corner(AlignedBox3f::TopRightCeil);

  m_kernel_radius = get_kernel_radius();
  m_rest_density = m_params->density;
  m_viscosity = m_params->viscosity;
  m_st = m_params->surface_tension;
  Real r = m_pc.get_radius();
  m_mass = m_params->density*8*r*r*r;
  m_recoil_velocity_damping = m_params->recoil_velocity_damping;

  if ( m_params->fluid_type == MCG03 )
    m_c2 = m_params->sound_speed * m_params->sound_speed;
  else //if ( m_params->fluid_type == BT07 )
  {
    // a heuristic to determine the speed of sound based on the maximum
    // possible velocity of a particle in the fluid
    float max_velocity2 =
      2.0f*M_G*(m_bmax[1] - m_bmin[1]) + m_params->velocity.squaredNorm();
    m_c2 = max_velocity2 / m_params->compressibility;
  }

  Matrix3XT<Real> &pos = m_pc.get_pos();
  m_accel.resizeLike(pos);
  m_dinv.resize(pos.cols());
  m_extern_accel.resizeLike(pos);
  m_vel.resizeLike(pos);
  for (Size i = 0; i < get_num_vertices(); ++i)
    m_vel.col(i) = m_params->velocity.template cast<Real>();
  m_dinv.setConstant(m_rest_density);
  reset_accel();

  //const clock_t begin_time = clock();
  //float d1 = PointCloud::compute_mindist();
  //float t1 = clock();
  //float d2 = PointCloud::compute_mindist_brute();
  //float t2 = clock();
  //qDebug() << float( t1 - begin_time ) / CLOCKS_PER_SEC << " to get " << d1;
  //qDebug() << float( t2 - t1 ) / CLOCKS_PER_SEC << " to get " << d2;

  // construct output filename (this can be done whenever)
  m_savefmt = global::dynset.savedir + "/" + m_params->saveprefix + "%03d.sim";
}

// clamp value d to min and max boundaries + epsilon,

bool 
Fluid::clamp(Real &d, Real min, Real max, Real tol)
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


void 
Fluid::resolve_collisions()
{
  for (Size i = 0; i < get_num_vertices(); ++i)
  {
    for (unsigned char j = 0; j < 3; ++j)
    {
      if (clamp(m_pc.pos_at(i)[j], m_bmin[j], m_bmax[j], 0.001))
      {
        vel_at(i)[j] *= -m_recoil_velocity_damping;
      }
    }
  }
}


void
Fluid::clamp(float adjust, float push)
{
  for (Size i = 0; i < get_num_vertices(); ++i)
  {
    for (unsigned char j = 0; j < 3; ++j)
    {
      Real &d = m_pc.pos_at(i)[j];
      Real &a = accel_at(i)[j];
      Real &v = vel_at(i)[j];
      Real min = m_bmin[j]-adjust;
      Real max = m_bmax[j]+adjust;
      if ( d < min )
      {
        d = min + push;
        if (a < 0)
          a = 0;
        if (v < 0)
          v = 0;
      }
      else if (d > max)
      {
        d = max - push;
        if (a > 0)
          a = 0;
        if (v > 0)
          v = 0;
      }
    }
  }
}

// File stuff

void
Fluid::save(unsigned int frame) 
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

  Size num = get_num_vertices();
  outfile.write(reinterpret_cast<const char *>(&num), sizeof num);
  float radius = m_pc.get_radius();
  float lev = 0.0f;

  for (Size i = 0; i < num; ++i)
  {
    Real * posptr = m_cache[frame].pos.data() + i*3;
    Real * velptr = m_cache[frame].vel.data() + i*3;
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


bool
Fluid::is_saved(unsigned int frame) 
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

bool
Fluid::load_saved_cache()
{
  // check if first frame was saved, if so load everything, otherwise ignore
  if (is_saved(0))
  {
    bool read_successful = true;
    for (unsigned int fr = 0; fr <= global::dynset.frames; ++fr)
    {
      read_successful &= read_saved(fr);
      glprintf_trcv(m_color, "\rLoading cache %d%%", 100*fr/global::dynset.frames);
    }
    glprintf_tr("\r   ");

    if (!read_successful)
      qWarning() << "Some cached files could not be found.";
    return true;
  } 
  return false;
}

// return if we successfully read from saved file and have a valid state

bool
Fluid::read_saved(unsigned int frame) 
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

  Size num;
  infile.read(reinterpret_cast<char *>(&num), sizeof num);

  m_cache[frame].pos.resize(NoChange,num);
  m_cache[frame].vel.resize(NoChange,num);

  for (Size i = 0; i < num; ++i)
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
    
    Real * posptr = m_cache[frame].pos.data() + i*3;
    Real * velptr = m_cache[frame].vel.data() + i*3;
    posptr[0] = Real(x); posptr[1] = Real(y); posptr[2] = Real(z);
    velptr[0] = Real(u); velptr[1] = Real(v); velptr[2] = Real(w);

    m_cache[frame].valid = true;
  }

  infile.close();

  return true; // last frame, so we're done
}


void
Fluid::clear_saved() 
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


void
Fluid::cache(unsigned int frame) 
{
  m_cache[frame].pos = m_pc.get_pos();
  m_cache[frame].vel = m_vel;
  m_cache[frame].valid = true;

  save(frame);
}


bool
Fluid::is_cached(unsigned int frame) 
{
  return m_cache[frame].valid;
}


bool
Fluid::load_cached(unsigned int frame) 
{
  // dont load if not cached or is first frame and next is not cached
  if ( !is_cached(frame) ) 
    return false;

  m_pc.set_pos(m_cache[frame].pos);
  m_vel = m_cache[frame].vel;

  return true;
}


void 
Fluid::clear_cache() 
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
