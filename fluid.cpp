#include <ctime>
#include <algorithm>
#include <limits>
#include "fluid.h"
#include "glpointcloud.h"
#include "gltext.h"

// Fluid stuff

// inflating the size of the grid
#define INFLATE 3.0 

// Constructor invoked when point cloud is instantiated as fluid initially
template<typename REAL, typename SIZE>
FluidRS<REAL,SIZE>::FluidRS(
    GLPointCloudRS<REAL, SIZE> *glpc,
    REAL density, REAL viscosity, REAL st)
  : PointCloudRS<REAL,SIZE>(*glpc->get_pointcloud())
  , m_bmin(this->m_bbox.corner(Eigen::AlignedBox3f::BottomLeftFloor))
  , m_bmax(this->m_bbox.corner(Eigen::AlignedBox3f::TopRightCeil))
  , m_kernel_radius(INFLATE*this->get_radius())
  , m_rest_density(density)
  , m_viscosity(4*viscosity)
  , m_st(st) // surface tension
  , m_mass(density*8*this->get_radius()*this->get_radius()*this->get_radius())
  // heuristic to compute speed of sound when objects have 0 initial velocity
  // at the start of simulation.
  // (WARNING: if this assumption changes, the heuristic becomes invalid!)
#ifdef MCG03
  , m_c2(40.0f)
#else
  , m_c2(2.0f*M_G*(m_bmax[1] - m_bmin[1]) / /*% density variation allowed*/ 0.01f)
#endif
  , m_glpc(glpc)
  , m_fluid_density_proc(m_kernel_radius)
  , m_fluid_density_update_proc(m_kernel_radius)
  , m_fluid_pressure_proc(m_kernel_radius)
  , m_fluid_viscosity_accel_proc(m_kernel_radius)
  , m_fluid_pressure_accel_proc(m_kernel_radius)
{
  m_accel.resizeLike(this->m_pos);
  m_extern_accel.resizeLike(this->m_pos);
  m_vel.resizeLike(this->m_pos);
  m_vel.setZero();
  reset_accel();

  copy_properties_to_proc(m_fluid_density_proc);
  copy_properties_to_proc(m_fluid_density_update_proc);
  copy_properties_to_proc(m_fluid_pressure_proc);
  copy_properties_to_proc(m_fluid_viscosity_accel_proc);
  copy_properties_to_proc(m_fluid_pressure_accel_proc);

  //const clock_t begin_time = clock();
  //float d1 = PointCloudRS<REAL,SIZE>::compute_mindist();
  //float t1 = clock();
  //float d2 = PointCloudRS<REAL,SIZE>::compute_mindist_brute();
  //float t2 = clock();
  //qDebug() << float( t1 - begin_time ) / CLOCKS_PER_SEC << " to get " << d1;
  //qDebug() << float( t2 - t1 ) / CLOCKS_PER_SEC << " to get " << d2;
}

// Constructor invoked when point cloud is converted to fluid
template<typename REAL, typename SIZE>
FluidRS<REAL,SIZE>::FluidRS(
    GLPointCloudRS<REAL, SIZE> *glpc,
    REAL density, REAL viscosity, REAL st)
  : PointCloudRS<REAL,SIZE>(*glpc->get_pointcloud())
  , m_bmin(this->m_bbox.corner(Eigen::AlignedBox3f::BottomLeftFloor))
  , m_bmax(this->m_bbox.corner(Eigen::AlignedBox3f::TopRightCeil))
  , m_kernel_radius(INFLATE*this->get_radius())
  , m_rest_density(density)
  , m_viscosity(4*viscosity)
  , m_st(st) // surface tension
  , m_mass(density*8*this->get_radius()*this->get_radius()*this->get_radius())
  // heuristic to compute speed of sound when objects have 0 initial velocity
  // at the start of simulation.
  // (WARNING: if this assumption changes, the heuristic becomes invalid!)
#ifdef MCG03
  , m_c2(40.0f)
#else
  , m_c2(2.0f*M_G*(m_bmax[1] - m_bmin[1]) / /*% density variation allowed*/ 0.01f)
#endif
  , m_glpc(glpc)
  , m_fluid_density_proc(m_kernel_radius)
  , m_fluid_density_update_proc(m_kernel_radius)
  , m_fluid_pressure_proc(m_kernel_radius)
  , m_fluid_viscosity_accel_proc(m_kernel_radius)
  , m_fluid_pressure_accel_proc(m_kernel_radius)
{
  m_accel.resizeLike(this->m_pos);
  m_extern_accel.resizeLike(this->m_pos);
  m_vel.resizeLike(this->m_pos);
  m_vel.setZero();
  reset_accel();

  copy_properties_to_proc(m_fluid_density_proc);
  copy_properties_to_proc(m_fluid_density_update_proc);
  copy_properties_to_proc(m_fluid_pressure_proc);
  copy_properties_to_proc(m_fluid_viscosity_accel_proc);
  copy_properties_to_proc(m_fluid_pressure_accel_proc);

  //const clock_t begin_time = clock();
  //float d1 = PointCloudRS<REAL,SIZE>::compute_mindist();
  //float t1 = clock();
  //float d2 = PointCloudRS<REAL,SIZE>::compute_mindist_brute();
  //float t2 = clock();
  //qDebug() << float( t1 - begin_time ) / CLOCKS_PER_SEC << " to get " << d1;
  //qDebug() << float( t2 - t1 ) / CLOCKS_PER_SEC << " to get " << d2;
}

template<typename REAL, typename SIZE>
FluidRS<REAL,SIZE>::~FluidRS()
{ }

template<typename REAL, typename SIZE>
template<class OutputType, class KernelType, class ComputeType>
inline void FluidRS<REAL,SIZE>::copy_properties_to_proc(
    CFQ<REAL,SIZE,OutputType,KernelType,ComputeType> &cfq_proc)
{
  cfq_proc.m_mass = get_mass();
  cfq_proc.m_radius = this->get_radius();
  cfq_proc.m_rest_density = get_rest_density();
  cfq_proc.m_viscosity = get_viscosity();
  cfq_proc.m_st = get_surface_tension();
  cfq_proc.m_cs2 = get_sound_speed2();
  cfq_proc.m_cs = std::sqrt(cfq_proc.m_cs2);
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
  for (SIZE i = 0; i < this->get_num_vertices(); ++i) // TODO: vectorize this
  {
    for (unsigned char j = 0; j < 3; ++j)
    {
      if (clamp(pos_at(i)[j], m_bmin[j], m_bmax[j]))
      {
        vel_at(i)[j] *= -VELOCITY_DAMPING;
      }
    }
  }
}


template<typename REAL, typename SIZE>
void FluidRS<REAL,SIZE>::update_data()
{
  if (m_glpc)
    m_glpc->update_data();
}

template class FluidRS<double, unsigned int>;
