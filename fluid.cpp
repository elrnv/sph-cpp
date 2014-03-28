#include <ctime>
#include <algorithm>
#include <limits>
#include "fluid.h"
#include "glpointcloud.h"
#include "gltext.h"

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
  qDebug() << "radius is" << r;
  m_kernel_radius = m_params->kernel_inflation * r;
  qDebug() << "kernel radius is" << m_kernel_radius;
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
  if (m_params->fluid_type != MCG03)
    return;

  for (SIZE i = 0; i < this->get_num_vertices(); ++i) // TODO: vectorize this
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
void FluidRS<REAL,SIZE>::update_data()
{
  if (m_glpc)
    m_glpc->update_data();
}


// Typed Fluid Stuff

template<typename REAL, typename SIZE, FluidType FT>
FluidRST<REAL,SIZE,FT>::FluidRST(const PointCloudRS<REAL,SIZE> *pc, FluidParamsPtr params)
  : FluidRS<REAL,SIZE>(*pc, params)
{ }

template<typename REAL, typename SIZE, FluidType FT>
FluidRST<REAL,SIZE,FT>::FluidRST(const aiMesh *pc, FluidParamsPtr params)
  : FluidRS<REAL,SIZE>(pc, params)
{ }

template<typename REAL, typename SIZE, FluidType FT>
FluidRST<REAL,SIZE,FT>::~FluidRST()
{ }

template<typename REAL, typename SIZE, FluidType FT>
void init_processors()
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

template<typename REAL, typename SIZE, FluidType FT>
template<class OutputType, class KernelType, class ComputeType>
inline void FluidRST<REAL,SIZE,FT>::copy_properties_to_proc(
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

template class FluidRST<double, unsigned int, MCG03>;
template class FluidRST<double, unsigned int, BT07>;
template class FluidRST<double, unsigned int, AIAST12>;
template class FluidRST<double, unsigned int, DEFAULT>;
