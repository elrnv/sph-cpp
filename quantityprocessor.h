#ifndef QUANTITYPROCESSOR_H
#define QUANTITYPROCESSOR_H

#include "particle.h"
#include "kernel.h"
#include "dynparams.h"

#define M_G 9.81f

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
    inline void fluid (ParticleType &p, FluidParticleRS<REAL,SIZE> &near_p)
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
  public ComputeQuantityRS<REAL, SIZE, FluidParticleRS<REAL,SIZE>, ComputeType>
{
  // global quantities acquired from the current observed object
public:
  REAL m_mass;
  REAL m_radius;
  REAL m_rest_density;
  REAL m_viscosity;
  REAL m_st;
  REAL m_cs2;
  REAL m_cs;
};

template<typename REAL, typename SIZE, class ComputeType>
class CBQ : public ComputeQuantityRS<REAL, SIZE, ParticleR<REAL>, ComputeType>
{ };


// Specific Compute SPH Quantity Processors

template<typename REAL, typename SIZE, int FT>
class CFDensityRST : public CFQ< REAL, SIZE, CFDensityRST<REAL,SIZE,FT> >
{
public:
  inline void init_kernel(float h) { m_kern.init(h); }
  inline void init(REAL &mv, REAL &av);
  inline void init_particle(FluidParticleRS<REAL,SIZE> &p);
  inline void fluid(FluidParticleRS<REAL,SIZE> &p, FluidParticleRS<REAL,SIZE> &near_p);
  inline void bound(FluidParticleRS<REAL,SIZE> &p, ParticleR<REAL> &near_p);
  inline void finish_particle(FluidParticleRS<REAL,SIZE> &p);

private:
  REAL *max_var; // max variation
  REAL *avg_var; // average variation
  Poly6Kernel m_kern;
};


template<typename REAL, typename SIZE, int FT>
class CFDensityUpdateRST : public CFQ<REAL,SIZE,CFDensityUpdateRST<REAL,SIZE,FT> >
{
public:
  inline void init_kernel(float h) { m_kern.init(h); }
  inline void init(float ts);
  inline void init_particle(FluidParticleRS<REAL,SIZE> &p);
  inline void fluid(FluidParticleRS<REAL,SIZE> &p, FluidParticleRS<REAL,SIZE> &near_p);
  inline void bound(FluidParticleRS<REAL,SIZE> &p, ParticleR<REAL> &near_p);
  inline void finish_particle(FluidParticleRS<REAL,SIZE> &p);
private:
  float m_timestep;
  Poly6GradKernel m_kern;
};

template<typename REAL, typename SIZE, int FT>
class CFAccelRST : public CFQ<REAL,SIZE,CFAccelRST<REAL,SIZE,FT> >
{
public:
  inline void init_kernel(float h);
  inline void init();
  inline void init_particle(FluidParticleRS<REAL,SIZE> &p);
  inline void fluid(FluidParticleRS<REAL,SIZE> &p, FluidParticleRS<REAL,SIZE> &near_p);
  inline void bound(FluidParticleRS<REAL,SIZE> &p, ParticleR<REAL> &near_p);
  inline void finish_particle(FluidParticleRS<REAL,SIZE> &p);
};

#endif // QUANTITYPROCESSOR_H
