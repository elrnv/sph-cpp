#ifndef FLUIDDATA_H
#define FLUIDDATA_H

#include "fluid.h"
#include "kernel.h"

// These are intermediate structs used to improve cache locality when doing sph
// computations. The values in these structs are used directly by the Particles
struct FluidData
{
  explicit FluidData(Size flidx, FluidVec &fluids)
    : flidx(flidx)
  { 
    Fluid &fl = fluids[flidx];
    m_mass = fl.get_mass();
    m_rest_density = fl.get_rest_density();
    m_viscosity = fl.get_viscosity();
    m_st = fl.get_surface_tension();
    m_cs2 = fl.get_sound_speed2();

    m_friction = fl.get_friction();
    m_cs = std::sqrt(m_cs2);
    m_compressibility = fl.get_compressibility();
    m_bmin = fl.get_bmin();
    m_bmax = fl.get_bmax();
  }

  Size flidx;

  Real m_mass;
  //Real m_radius;
  Real m_rest_density;
  Real m_viscosity;
  Real m_st;
  Real m_cs2; // TODO: test if this is actually needed or just use m_cs

  Real m_friction;
  Real m_cs; // square root of m_cs2
  Real m_compressibility;

  Vector3f m_bmin;
  Vector3f m_bmax;

#ifdef REPORT_DENSITY_VARIATION
  Real *m_max_var;
  Real *m_avg_var;
#endif
}; // FluidT


template<int PT>
struct FluidDataT : public FluidData
{
  explicit FluidDataT(Size flidx, FluidVec &fluids) 
    : FluidData(flidx, fluids) { }

  void init_kernel(float h)
  {
    m_kern.init(h);
  }

  Poly6Kernel m_kern;
}; // FluidDataT

template<>
struct FluidDataT<MCG03> : public FluidData
{
  explicit FluidDataT(Size flidx, FluidVec &fluids) 
    : FluidData(flidx, fluids) { }

  void init_kernel(float h)
  {
    m_kern.init(h);
    m_spikygrad_kern.init(h);
    m_visclap_kern.init(h);
    m_color_kern.init(h);
  }

  Poly6Kernel       m_kern; // for density

  SpikyGradKernel   m_spikygrad_kern; // for pressure
  ViscLapKernel     m_visclap_kern;   // for viscosity
  CubicSplineKernel m_color_kern;     // for surface tension
}; // FluidDataT<MCG03>

template<>
struct FluidDataT<BT07> : public FluidData
{
  explicit FluidDataT(Size flidx, FluidVec &fluids) 
    : FluidData(flidx, fluids) { }

  void init_kernel(float h)
  {
    m_kern.init(h);
    m_grad_kern.init(h);
    m_bound_kern.init(h);
  }

  CubicSplineKernel m_kern;

  CubicSplineGradKernel m_grad_kern;
  MKI04Kernel           m_bound_kern;
}; // FluidDataT<BT07>

// templated "typedef"
template<int PT>
using FluidDataVecT = std::vector< FluidDataT<PT> >;

#endif // FLUIDDATA_H
