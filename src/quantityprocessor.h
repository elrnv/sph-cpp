#ifndef QUANTITYPROCESSOR_H
#define QUANTITYPROCESSOR_H

#include "types.h"
#include "kernel.h"

class Fluid;
struct Particle;
struct FluidParticle;

// Compute SPH quantity definitions
template<typename ParticleType, class ComputeType>
class ComputeQuantity
{
  public:
    // Compute Interface
    inline void init_particle (ParticleType &p)
    {
      static_cast<ComputeType*>(this)->init_particle(p);
    }
    inline void fluid (ParticleType &p, FluidParticle &near_p)
    {
      static_cast<ComputeType*>(this)->fluid(p, near_p);
    }
    inline void bound (ParticleType &p, Particle &near_p)
    {
      static_cast<ComputeType*>(this)->bound(p, near_p);
    }
    inline void finish_particle (ParticleType &p)
    {
      static_cast<ComputeType*>(this)->finish_particle(p);
    }
};

template<class ComputeType>
class CFQ : 
  public ComputeQuantity<FluidParticle, ComputeType>
{
  // global quantities acquired from the current observed object
public:
  Real m_mass;
  Real m_radius;
  Real m_rest_density;
  Real m_viscosity;
  Real m_st;
  Real m_cs2;
  Real m_cs;
  Real m_compressibility;
  Real m_friction;
  Vector3f m_bmin;
  Vector3f m_bmax;
};

template<class ComputeType>
class CBQ : public ComputeQuantity<Particle, ComputeType>
{ };


// Specific Compute SPH Quantity Processors

template<int FT>
class CFDensityT : public CFQ<CFDensityT<FT> >
{

public:
  void init_kernel(float h);
  void init(Real &mv, Real &av);
  void init_particle(FluidParticle &p);
  void fluid(FluidParticle &p, FluidParticle &near_p);
  void bound(FluidParticle &p, Particle &near_p);
  void finish_particle(FluidParticle &p);

private:
  Real *max_var; // max variation
  Real *avg_var; // average variation
  Poly6Kernel m_kern;
};


template<int FT>
class CFDensityUpdateT : public CFQ<CFDensityUpdateT<FT> >
{
public:
  void init_kernel(float h);
  void init(float ts);
  void init_particle(FluidParticle &p);
  void fluid(FluidParticle &p, FluidParticle &near_p);
  void bound(FluidParticle &p, Particle &near_p);
  void finish_particle(FluidParticle &p);
private:
  float m_timestep;
  CubicSplineGradKernel m_kern;
};

template<int FT>
class CFAccelT : public CFQ<CFAccelT<FT> >
{
public:
  void init_kernel(float h);
  void init_particle(FluidParticle &p);
  void fluid(FluidParticle &p, FluidParticle &near_p);
  void bound(FluidParticle &p, Particle &near_p);
  void finish_particle(FluidParticle &p);

private:
  SpikyGradKernel m_spikygrad_kern; // for pressure
  ViscLapKernel   m_visclap_kern;   // for viscosity
  Poly6GradKernel m_colorgrad_kern;     // surface tension
  Poly6LapKernel  m_colorlap_kern;     // surface tension
};

template<int FT>
class CFSurfaceNormalT
  : public CFQ<CFSurfaceNormalT<FT> >
{
public:
  void init_kernel(float h);
  void init_particle(FluidParticle &p);
  void fluid(FluidParticle &p, FluidParticle &near_p);
  void bound(FluidParticle &p, Particle &near_p);
  void finish_particle(FluidParticle &p);

private:
  CubicSplineGradKernel m_grad_kern;
};

template<int FT>
class CFSurfaceTensionT 
  : public CFQ<CFSurfaceTensionT<FT> >
{
public:
  void init_kernel(float h);
  void init_particle(FluidParticle &p);
  void fluid(FluidParticle &p, FluidParticle &near_p);
  void bound(FluidParticle &p, Particle &near_p);
  void finish_particle(FluidParticle &p);

private:
  CubicSplineGradKernel m_grad_kern;
  CubicSplineKernel  m_kern;
};



// Only used for IISPH
template<int FT>
class CFPressureAccelT : public CFQ<CFPressureAccelT<FT> >
{
public:
  void init_kernel(float h);
  void init_particle(FluidParticle &p);
  void fluid(FluidParticle &p, FluidParticle &near_p);
  void bound(FluidParticle &p, Particle &near_p);
  void finish_particle(FluidParticle &p);

private:
  CubicSplineGradKernel m_grad_kern;
  CubicSplineKernel     m_st_kern;
};

template<int FT>
class CFPrepareJacobiT 
  : public CFQ<CFPrepareJacobiT<FT> >
{
public:
  void init_kernel(float h);
  void init(float ts);
  void init_particle(FluidParticle &p);
  void fluid(FluidParticle &p, FluidParticle &near_p);
  void bound(FluidParticle &p, Particle &near_p);
  void finish_particle(FluidParticle &p);

private:
  float dt;
  CubicSplineGradKernel m_grad_kern;
  CubicSplineKernel  m_kern;
};

template<int FT>
class CFJacobiSolveFirstT 
  : public CFQ<CFJacobiSolveFirstT<FT> >
{
public:
  void init_kernel(float h);
  void init(float ts);
  void init_particle(FluidParticle &p);
  void fluid(FluidParticle &p, FluidParticle &near_p);
  void bound(FluidParticle &p, Particle &near_p);
  void finish_particle(FluidParticle &p);

private:
  float dt;
  CubicSplineGradKernel m_grad_kern;
  CubicSplineKernel  m_kern;
};

template<int FT>
class CFJacobiSolveSecondT
  : public CFQ<CFJacobiSolveSecondT<FT> >
{
public:
  void init_kernel(float h);
  void init(float ts, Real &avg_d, Real &avg_p);
  void init_particle(FluidParticle &p);
  void fluid(FluidParticle &p, FluidParticle &near_p);
  void bound(FluidParticle &p, Particle &near_p);
  void finish_particle(FluidParticle &p);

private:
  float dt;
  Real *avg_density;
  Real *avg_pressure;
  CubicSplineGradKernel m_grad_kern;
  CubicSplineKernel  m_kern;
};

#endif // QUANTITYPROCESSOR_H
