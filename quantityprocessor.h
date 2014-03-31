#ifndef QUANTITYPROCESSOR_H
#define QUANTITYPROCESSOR_H

#include "particle.h"
#include "kernel.h"
#include "dynparams.h"

#define M_G 9.81f

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
    inline void init_kernel(float _h) { m_kern.init(_h); }

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

template<typename REAL, typename OutputType, class KernelType, class ComputeType>
class CBQ : public ComputeQuantityR<REAL, ParticleR<REAL>, OutputType, KernelType, ComputeType>
{
};

// convenience defines used in computing SPH values
#define CFQ_TYPEDEF template<typename REAL, typename SIZE, class ComputeType> using
CFQ_TYPEDEF CFQPoly6RS     = CFQ<REAL, SIZE, double, Poly6Kernel, ComputeType>;
CFQ_TYPEDEF CFQPoly6GradRS = CFQ<REAL, SIZE, Vector3d, Poly6GradKernel, ComputeType>;
CFQ_TYPEDEF CFQSpikyGradRS = CFQ<REAL, SIZE, Vector3d, SpikyGradKernel, ComputeType>;
CFQ_TYPEDEF CFQViscLapRS   = CFQ<REAL, SIZE, double, ViscLapKernel, ComputeType>;

#define CBQ_TYPEDEF template<typename REAL, class ComputeType> using
CBQ_TYPEDEF CBQPoly6R = CBQ<REAL, double, Poly6Kernel, ComputeType>;


// Specific Compute SPH Quantity Processors

template<typename REAL, typename SIZE, int FT>
class CFDensityRST :
  public CFQPoly6RS< REAL, SIZE, CFDensityRST<REAL,SIZE,FT> >
{
public:
  inline void init(REAL &mv, REAL &av);
  inline void init_particle(FluidParticleR<REAL> &p);
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p);
  inline void finish_particle(FluidParticleR<REAL> &p);

private:
  REAL *max_var; // max variation
  REAL *avg_var; // average variation
};


template<typename REAL, typename SIZE, int FT>
class CFDensityUpdateRST : 
  public CFQPoly6GradRS<REAL,SIZE,CFDensityUpdateRST<REAL,SIZE,FT> >
{
public:
  inline void init(float ts);
  inline void init_particle(FluidParticleR<REAL> &p);
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p);
  inline void finish_particle(FluidParticleR<REAL> &p);
private:
  float m_timestep;
};

template<typename REAL, typename SIZE, int FT>
class CFPressureRST :
  public CFQPoly6RS<REAL,SIZE,CFPressureRST<REAL,SIZE,FT> >
{
public:
  inline REAL pow7(REAL x);
  inline void init_particle(ParticleR<REAL> &p);
  inline void fluid(ParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  inline void bound(ParticleR<REAL> &p, ParticleR<REAL> &near_p);
  inline void finish_particle(ParticleR<REAL> &p);
};

template<typename REAL, typename SIZE, int FT>
class CFPressureAccelRST :
  public CFQSpikyGradRS<REAL,SIZE,CFPressureAccelRST<REAL,SIZE,FT> >
{
public:
  inline void init();
  inline void init_particle(FluidParticleR<REAL> &p);
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p);
  inline void finish_particle(FluidParticleR<REAL> &p);
private:
  MKI04Kernel m_bound_kern;
};


template<typename REAL, typename SIZE, int FT>
class CFViscosityAccelRST : 
  public CFQSpikyGradRS<REAL,SIZE,CFViscosityAccelRST<REAL,SIZE,FT> >
{
public:
  inline void init_particle(FluidParticleR<REAL> &p);
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p);
  inline void finish_particle(FluidParticleR<REAL> &p);
};

template<typename REAL, typename SIZE, int FT>
class CFSurfaceTensionAccelRST : 
  public CFQPoly6GradRS<REAL,SIZE,CFSurfaceTensionAccelRST<REAL,SIZE,FT> >
{
public:
  inline void init();
  inline void init_particle(FluidParticleR<REAL> &p);
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p);
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p);
  inline void finish_particle(FluidParticleR<REAL> &p);
private:
  Poly6LapKernel m_lap_kern;
};


#endif // QUANTITYPROCESSOR_H
