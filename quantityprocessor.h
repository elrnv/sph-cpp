#ifndef QUANTITYPROCESSOR_H
#define QUANTITYPROCESSOR_H

#include "particle.h"
#include "kernel.h"

template <typename REAL, typename SIZE>
class FluidRS;

// Compute SPH quantity definitions
template< typename REAL,
          typename ParticleType,
          typename OutputType, 
          class KernelType,
          class ComputeType>
class ComputeQuantityR
{
  public:
    ComputeQuantityR(float h) : m_kern(h) { }

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

  public: //TODO: make protected
    // Smoothing kernel used to interpolate data
    Kernel<OutputType, KernelType> m_kern;
};

template<typename REAL, typename SIZE,
         typename OutputType,
         class KernelType,
         class ComputeType>
class CFQ : 
  public ComputeQuantityR<REAL, FluidParticleR<REAL>, OutputType, KernelType, ComputeType>
{
public:
  CFQ(float h) 
    : ComputeQuantityR<REAL, FluidParticleR<REAL>, OutputType, KernelType, ComputeType>(h)
  { }

  friend FluidRS<REAL, SIZE>;
  // global quantities acquired from the current observed object
// TODO: protect the following
  REAL m_mass;
  REAL m_radius;
  REAL m_rest_density;
  REAL m_viscosity;
  REAL m_st;
  REAL m_cs2;
  REAL m_cs;
protected:
};

template<typename REAL, typename OutputType, class KernelType, class ComputeType>
class CBQ : public ComputeQuantityR<REAL, ParticleR<REAL>, OutputType, KernelType, ComputeType>
{
public:
  CBQ(float h) : ComputeQuantityR<REAL, ParticleR<REAL>, OutputType, KernelType, ComputeType>(h) { }
};

// convenience defines used in computing SPH values
#define CFQ_TYPEDEF template<typename REAL, typename SIZE, class ComputeType> using
CFQ_TYPEDEF CFQPoly6RS     = CFQ<REAL, SIZE, double, Poly6Kernel, ComputeType>;
CFQ_TYPEDEF CFQPoly6GradRS = CFQ<REAL, SIZE, Vector3d, Poly6GradKernel, ComputeType>;
CFQ_TYPEDEF CFQSpikyGradRS = CFQ<REAL, SIZE, Vector3d, SpikyGradKernel, ComputeType>;
CFQ_TYPEDEF CFQViscLapRS   = CFQ<REAL, SIZE, double, ViscLapKernel, ComputeType>;

#define CBQ_TYPEDEF template<typename REAL, class ComputeType> using
CBQ_TYPEDEF CBQPoly6R = CBQ<REAL, double, Poly6Kernel, ComputeType>;

#endif // QUANTITYPROCESSOR_H
