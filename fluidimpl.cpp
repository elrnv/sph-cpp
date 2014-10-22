#include "fluidimpl.h"

// Typed Fluid Stuff

template<int FT>
FluidT<FT>::FluidT(PointCloud &pc, FluidParamsPtr params)
  : Fluid(pc, params)
{ }

template<int FT>
FluidT<FT>::~FluidT()
{ }

  /*
template<int FT> template<class ComputeType>
inline void 
FluidT<FT>::load_proc_params(CFQ<ComputeType> &cfq)
{
  cfq.m_mass = get_mass();
  cfq.m_radius = get_radius();
  cfq.m_rest_density = get_rest_density();
  cfq.m_viscosity = get_viscosity();
  cfq.m_st = get_surface_tension();
  cfq.m_cs2 = get_sound_speed2();
  cfq.m_cs = std::sqrt(cfq.m_cs2);
  cfq.m_compressibility = get_compressibility();
  cfq.m_friction = get_friction();
  cfq.m_bmin = get_bmin();
  cfq.m_bmax = get_bmax();
}

template<int FT>
inline void 
FluidT<FT>::init_processors()
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

  load_proc_params(m_fluid_density_proc);
  load_proc_params(m_fluid_density_update_proc);
  load_proc_params(m_fluid_accel_proc);
  load_proc_params(m_fluid_surface_normal_proc);
  load_proc_params(m_fluid_surface_tension_proc);
  load_proc_params(m_fluid_prepare_jacobi_proc);
  load_proc_params(m_fluid_jacobi_solve1_proc);
  load_proc_params(m_fluid_jacobi_solve2_proc);
  load_proc_params(m_fluid_pressure_accel_proc);
}
  */

template class FluidT<0>; //MCG03
template class FluidT<1>; //BT07
template class FluidT<2>; //ICS13
template class FluidT<3>; //DEFAULT
