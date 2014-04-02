#ifndef QUANTITYPROCESSOR_H
#define QUANTITYPROCESSOR_H

#include "particle.h"
#include "kernel.h"
#include "dynparams.h"

#define M_G 9.81f

template<typename REAL>
inline REAL pow7(REAL x) { return (x*x)*(x*x)*(x*x)*x; }

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
  }

protected:
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
  CFDensityRST(float h) { m_kern.init(h); }
  inline void init_kernel(float h) { m_kern.init(h); }
  inline void init(REAL &mv, REAL &av) 
  {
    max_var = &mv; avg_var = &av;
  }
  inline void init_particle(FluidParticleR<REAL> &p)

{ 
  //fprintf(stderr, "density init_particle\n");
  p.dinv = 0.0f; p.vol = 0.0f;
}
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)

{
  p.dinv += this->m_kern[ p.pos - near_p.pos ];
  p.vol += this->m_kern[ p.pos - near_p.pos ];
 // fprintf(stderr, "p.dinv += %.2e\n", p.dinv);
}
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)

{
  if (FT == MCG03)
    return;

  p.dinv += this->m_rest_density * near_p.dinv * this->m_kern[ p.pos - near_p.pos ];
  p.vol += this->m_kern[ p.pos - near_p.pos ];
}
  inline void finish_particle(FluidParticleR<REAL> &p)
{
  if (FT == MCG03)
  {
    p.dinv = 1.0f/(this->m_mass * p.dinv * this->m_kern.coef);
  }
  else
  {
    // By default, use kernel correction
    p.dinv =
      (8*this->m_radius*this->m_radius*this->m_radius*p.vol)/(this->m_mass * p.dinv);
#if 0
  qDebug() << (this->m_mass * p.dinv) << "/"
    << (8*this->m_radius*this->m_radius*this->m_radius*p.vol) << " = " <<
    (this->m_mass * p.dinv)/(8*this->m_radius*this->m_radius*this->m_radius*p.vol);
#endif
  }

  REAL var = std::abs(1.0f/p.dinv - this->m_rest_density);
  if ( var > *max_var )
    *max_var = var;
  *avg_var += var;

  if (FT == MCG03)
  {
    p.pressure = this->m_cs2 * (1.0f/p.dinv - this->m_rest_density);
  }
  else // Enable tait pressure equation by default
  {
    p.pressure =
      this->m_rest_density * this->m_cs2 * 0.14285714285714 * 
      (pow7(1.0f / (p.dinv * this->m_rest_density)) - 1);
  }
//  fprintf(stderr, "p.pressure = %.2e\n", p.pressure);
}

private:
  REAL *max_var; // max variation
  REAL *avg_var; // average variation
  Poly6Kernel m_kern;
};


template<typename REAL, typename SIZE, int FT>
class CFDensityUpdateRST : public CFQ<REAL,SIZE,CFDensityUpdateRST<REAL,SIZE,FT> >
{
public:
  inline void init_kernel(float h) { }
  inline void init(float ts) { }
  inline void init_particle(FluidParticleR<REAL> &p) { }
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p) { }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p) { }
  inline void finish_particle(FluidParticleR<REAL> &p) { }
private:
  float m_timestep;
  Poly6GradKernel m_kern;
};

template<typename REAL, typename SIZE, int FT>
class CFAccelRST : public CFQ<REAL,SIZE,CFAccelRST<REAL,SIZE,FT> >
{
public:
  CFAccelRST(float h)
  {
    m_spikygrad_kern.init(h);
    m_visclap_kern.init(h);
    m_colorgrad_kern.init(h);
  }
  inline void init_kernel(float h)
  {
    m_spikygrad_kern.init(h);
    m_visclap_kern.init(h);
    m_colorgrad_kern.init(h);
  }
  inline void init_particle(FluidParticleR<REAL> &p)
  { Q_UNUSED(p); }

  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;
    Vector3R<REAL> res(
    // pressure
         -p.dinv*this->m_mass*(p.pressure +
          near_p.pressure)*0.5*near_p.dinv*m_spikygrad_kern(p.pos -
            near_p.pos)
         +
    // viscosity
          this->m_mass * p.dinv * this->m_viscosity *
        near_p.dinv * (near_p.vel - p.vel) * this->m_visclap_kern(p.pos - near_p.pos));

    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p) { }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    p.extern_accel[1] -= M_G;
    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += p.extern_accel[i];
    //qDebug() << "p.accel = " << p.accel[0] << p.accel[1] << p.accel[2];
  }

private:
  SpikyGradKernel m_spikygrad_kern; // for pressure
  ViscLapKernel   m_visclap_kern;   // for viscosity
  Poly6GradKernel m_colorgrad_kern;     // surface tension
};

#endif // QUANTITYPROCESSOR_H
