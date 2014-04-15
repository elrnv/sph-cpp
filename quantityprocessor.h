#ifndef QUANTITYPROCESSOR_H
#define QUANTITYPROCESSOR_H

#include "kernel.h"

template<typename REAL>
struct ParticleR;
template<typename REAL>
struct FluidParticleR;

template <typename REAL, typename SIZE>
class FluidRS;

// Compute SPH quantity definitions
template< typename REAL, typename SIZE,
          typename ParticleType,
          class ComputeType>
class ComputeQuantityRS
{
  public:
    // Compute Interface
    inline void init_particle (ParticleType &p)
    {
      static_cast<ComputeType*>(this)->init_particle(p);
    }
    inline void fluid (ParticleType &p, FluidParticleR<REAL> &near_p)
    {
      static_cast<ComputeType*>(this)->fluid(p, near_p);
    }
    inline void bound (ParticleType &p, ParticleR<REAL> &near_p)
    {
      static_cast<ComputeType*>(this)->bound(p, near_p);
    }
    inline void finish_particle (ParticleType &p)
    {
      static_cast<ComputeType*>(this)->finish_particle(p);
    }
};

template<typename REAL, typename SIZE, class ComputeType>
class CFQ : 
  public ComputeQuantityRS<REAL, SIZE, FluidParticleR<REAL>, ComputeType>
{
  // global quantities acquired from the current observed object
public:
  void copy_fluid_params(FluidRS<REAL,SIZE> &fl)
  {
    m_mass = fl.get_mass();
    m_radius = fl.get_radius();
    m_rest_density = fl.get_rest_density();
    m_viscosity = fl.get_viscosity();
    m_st = fl.get_surface_tension();
    m_cs2 = fl.get_sound_speed2();
    m_cs = std::sqrt(m_cs2);
    m_compressibility = fl.get_compressibility();
    m_friction = fl.get_friction();
    m_bmin = fl.get_bmin();
    m_bmax = fl.get_bmax();
  }

protected:
  REAL m_mass;
  REAL m_radius;
  REAL m_rest_density;
  REAL m_viscosity;
  REAL m_st;
  REAL m_cs2;
  REAL m_cs;
  REAL m_compressibility;
  REAL m_friction;
  Vector3f m_bmin;
  Vector3f m_bmax;
};

template<typename REAL, typename SIZE, class ComputeType>
class CBQ : public ComputeQuantityRS<REAL, SIZE, ParticleR<REAL>, ComputeType>
{ };


// Specific Compute SPH Quantity Processors

template<typename REAL, typename SIZE, int FT>
class CFDensityRST : public CFQ< REAL, SIZE, CFDensityRST<REAL,SIZE,FT> >
{

public:
  void init_kernel(float h);
  void init(REAL &mv, REAL &av);
  void init_particle(FluidParticleR<REAL> &p);
  void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p);
  void finish_particle(FluidParticleR<REAL> &p);

private:
  REAL *max_var; // max variation
  REAL *avg_var; // average variation
  Poly6Kernel m_kern;
};


template<typename REAL, typename SIZE, int FT>
class CFDensityUpdateRST : public CFQ<REAL,SIZE,CFDensityUpdateRST<REAL,SIZE,FT> >
{
public:
  void init_kernel(float h);
  void init(float ts);
  void init_particle(FluidParticleR<REAL> &p);
  void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p);
  void finish_particle(FluidParticleR<REAL> &p);
private:
  float m_timestep;
  CubicSplineGradKernel m_kern;
};

template<typename REAL, typename SIZE, int FT>
class CFAccelRST : public CFQ<REAL,SIZE,CFAccelRST<REAL,SIZE,FT> >
{
public:
  void init_kernel(float h);
  void init_particle(FluidParticleR<REAL> &p);
  void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p);
  void finish_particle(FluidParticleR<REAL> &p);

private:
  SpikyGradKernel m_spikygrad_kern; // for pressure
  ViscLapKernel   m_visclap_kern;   // for viscosity
  Poly6GradKernel m_colorgrad_kern;     // surface tension
  Poly6LapKernel  m_colorlap_kern;     // surface tension
};

template<typename REAL, typename SIZE, int FT>
class CFSurfaceNormalRST 
  : public CFQ<REAL,SIZE,CFSurfaceNormalRST<REAL,SIZE,FT> >
{
public:
  void init_kernel(float h);
  void init_particle(FluidParticleR<REAL> &p);
  void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p);
  void finish_particle(FluidParticleR<REAL> &p);

private:
  CubicSplineGradKernel m_grad_kern;
};

template<typename REAL, typename SIZE, int FT>
class CFSurfaceTensionRST 
  : public CFQ<REAL,SIZE,CFSurfaceTensionRST<REAL,SIZE,FT> >
{
public:
  void init_kernel(float h);
  void init_particle(FluidParticleR<REAL> &p);
  void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p);
  void finish_particle(FluidParticleR<REAL> &p);

private:
  CubicSplineGradKernel m_grad_kern;
  CubicSplineKernel  m_kern;
};



// Only used for IISPH
template<typename REAL, typename SIZE, int FT>
class CFPressureAccelRST : public CFQ<REAL,SIZE,CFPressureAccelRST<REAL,SIZE,FT> >
{
public:
  void init_kernel(float h);
  void init_particle(FluidParticleR<REAL> &p);
  void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p);
  void finish_particle(FluidParticleR<REAL> &p);

private:
  CubicSplineGradKernel m_grad_kern;
  CubicSplineKernel     m_st_kern;
};

template<typename REAL, typename SIZE, int FT>
class CFPrepareJacobiRST 
  : public CFQ<REAL,SIZE,CFPrepareJacobiRST<REAL,SIZE,FT> >
{
public:
  void init_kernel(float h);
  void init(float ts);
  void init_particle(FluidParticleR<REAL> &p);
  void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p);
  void finish_particle(FluidParticleR<REAL> &p);

private:
  float dt;
  CubicSplineGradKernel m_grad_kern;
  CubicSplineKernel  m_kern;
};

template<typename REAL, typename SIZE, int FT>
class CFJacobiSolveFirstRST 
  : public CFQ<REAL,SIZE,CFJacobiSolveFirstRST<REAL,SIZE,FT> >
{
public:
  void init_kernel(float h);
  void init(float ts);
  void init_particle(FluidParticleR<REAL> &p);
  void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p);
  void finish_particle(FluidParticleR<REAL> &p);

private:
  float dt;
  CubicSplineGradKernel m_grad_kern;
  CubicSplineKernel  m_kern;
};

template<typename REAL, typename SIZE, int FT>
class CFJacobiSolveSecondRST 
  : public CFQ<REAL,SIZE,CFJacobiSolveSecondRST<REAL,SIZE,FT> >
{
public:
  void init_kernel(float h);
  void init(float ts);
  void init_particle(FluidParticleR<REAL> &p);
  void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p);
  void finish_particle(FluidParticleR<REAL> &p);

private:
  float dt;
  CubicSplineGradKernel m_grad_kern;
  CubicSplineKernel  m_kern;
};

#endif // QUANTITYPROCESSOR_H
