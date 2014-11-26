#ifndef FLUIDDATA_H
#define FLUIDDATA_H

#include "fluid.h"

// These are intermediate structs used to improve cache locality when doing sph
// computations. The values in these structs are used directly by the Particles
struct FluidData
{
  explicit FluidData(Fluid &fl)
    : fl(fl)
  { 
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

  Fluid &fl;

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
  explicit FluidDataT(Fluid &fl)
    : fl(fl)
  { }

}; // FluidDataT

#endif // FLUIDDATA_H
