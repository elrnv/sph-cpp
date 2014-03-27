#ifndef QUANTITYPROCESSOR_H
#define QUANTITYPROCESSOR_H

#include "particle.h"
#include "kernel.h"

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
public:
  friend FluidRS<REAL, SIZE>;

  // global quantities acquired from the current observed object
protected:
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

// General Density processor
template<typename REAL, typename SIZE, typename SpecificType>
class CFDensityRS :
  public CFQPoly6RS< REAL, SIZE, CFDensityRS<REAL,SIZE> >
{
public:
  void init(REAL &mv, REAL &av) { max_var = &mv; avg_var = &av; }

  inline void init_particle(FluidParticleR<REAL> &p)
  { 
    p.dinv = 0.0f; p.vol = 0.0f;
  }
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    p.dinv += this->m_kern[ p.pos - near_p.pos ];
    p.vol += this->m_kern[ p.pos - near_p.pos ];
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    static_cast<SpecificType*>(this)->bound(p, near_p);
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    static_cast<SpecificType*>(this)->finish_particle(p);

    REAL var = std::abs(1.0f/p.dinv - this->m_rest_density);
    if ( var > *max_var )
      *max_var = var;
    *avg_var += var;
  }

private:
  REAL *max_var; // max variation
  REAL *avg_var; // average variation
};

template<typename REAL, typename SIZE>
class CFDensityMCG03RS :
  public CFDensityRS<REAL, SIZE, CFDensityMCG03RS<REAL, SIZE> >
{
public:
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p) { }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    p.dinv = 1.0f/(this->m_mass * p.dinv * this->m_kern.coef);
    // qDebug() << (this->m_mass * p.dinv * this->m_kern.coef);
  }
};

template<typename REAL, typename SIZE>
class CFDensityMCG03RS :
  public CFDensityRS<REAL, SIZE, CFDensityMCG03RS<REAL, SIZE> >
{
public:
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
#ifdef BOUNDARY_PARTICLE
    p.dinv += this->m_rest_density * near_p.dinv * this->m_kern[ p.pos - near_p.pos ];
    p.vol += this->m_kern[ p.pos - near_p.pos ];
#endif
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
#ifdef DENSITY_KERNEL_CORRECTION
    p.dinv =
      (8*this->m_radius*this->m_radius*this->m_radius*p.vol)/(this->m_mass * p.dinv);
#if 0
    qDebug() << (this->m_mass * p.dinv) << "/"
    << (8*this->m_radius*this->m_radius*this->m_radius*p.vol) << " = " <<
      (this->m_mass * p.dinv)/(8*this->m_radius*this->m_radius*this->m_radius*p.vol);
#endif
#else
    p.dinv = 1.0f/(this->m_mass * p.dinv * this->m_kern.coef);
#if 0
    qDebug() << (this->m_mass * p.dinv * this->m_kern.coef);
#endif
#endif
  }
};







template<typename REAL, typename SIZE>
class CFDensityUpdateRS : 
  public CFQPoly6GradRS<REAL,SIZE,CFDensityUpdateRS<REAL,SIZE> >
{
public:
  void init(float ts) { m_timestep = ts; }

  inline void init_particle(FluidParticleR<REAL> &p)
  { p.temp = 0.0f; }
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    p.temp += (near_p.vel - p.vel).dot(this->m_kern[ p.pos - near_p.pos ]);
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    p.temp += (-p.vel).dot(this->m_kern[ p.pos - near_p.pos ]);
    //p.dinv += rest_density * near_p.dinv * kern[ p.pos - near_p.pos ];
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    p.dinv += m_timestep * p.dinv * p.dinv * p.temp * this->m_kern.coef;
    qDebug() << "p.dinv: " << p.dinv;
  }

private:
  float m_timestep;
};

template<typename REAL, typename SIZE>
class CFPressureRS :
  public CFQPoly6RS<REAL,SIZE,CFPressureRS<REAL,SIZE> >
{
public:
  inline REAL pow7(REAL x) { return (x*x)*(x*x)*(x*x)*x; }

  inline void init_particle(ParticleR<REAL> &p)
  { 
//    p.dinv = 0.0f;
  }
  inline void fluid(ParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
 //   p.dinv += this->m_kern[ p.pos - near_p.pos ];
  }
  inline void bound(ParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
  //  p.dinv += this->m_rest_density * near_p.dinv * this->m_kern[ p.pos - near_p.pos ];
  }
  inline void finish_particle(ParticleR<REAL> &p)
  {
    //REAL density = this->m_mass * p.dinv * this->m_kern.coef;
    //p.dinv = 1.0f/density;
#if defined TAIT_PRESSURE
    p.pressure =
      this->m_rest_density * this->m_cs2 * 0.14285714285714 * 
      (pow7(1.0f / (p.dinv * this->m_rest_density)) - 1);
#elif defined IDEAL_GAS_PRESSURE
    p.pressure =
      this->m_cs2 *
      (1.0f / p.dinv  -  this->m_rest_density);
#endif
  //  qDebug() << "density:" << 1.0f/p.dinv << " rest_density:" <<
   //   this->m_rest_density << "  pressure:" << p.pressure;
  }
};

#ifdef MCG03_VISCOSITY_FORCE
  CFQ_TYPEDEF CFQViscKernRS = CFQViscLapRS<REAL,SIZE,ComputeType>; // used to compute viscosity force
#else
  CFQ_TYPEDEF CFQViscKernRS = CFQSpikyGradRS<REAL,SIZE,ComputeType>; // used to compute viscosity force
#endif

template<typename REAL, typename SIZE>
class CFViscosityAccelRS : 
  public CFQViscKernRS<REAL,SIZE,CFViscosityAccelRS<REAL,SIZE> >
{
public:
  inline void init_particle(FluidParticleR<REAL> &p)
  { Q_UNUSED(p); }
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;
#ifdef MCG03_VISCOSITY_FORCE
    Vector3R<REAL> res(
          this->m_mass*p.dinv*this->m_viscosity*near_p.dinv*(near_p.vel - p.vel)*this->m_kern(p.pos - near_p.pos) // viscosity contribution
        );
#else
    Vector3R<REAL> x_ab = p.pos - near_p.pos;
    REAL vx = x_ab.dot(p.vel - near_p.vel);
    if (vx >= 0)
      return;

    REAL nu = 2*this->m_viscosity*this->m_radius;
    //REAL pab = -
    Vector3R<REAL> res(
          this->m_mass*this->m_kern(x_ab) // viscosity contribution
        );
#endif

    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    //Vector3R<REAL> res(
    //    + viscosity*near_p.dinv*(near_p.vel - p.vel)*kern(p.pos - near_p.pos) // viscosity contribution
    //    );

    //for (unsigned char i = 0; i < 3; ++i)
    //  p.extern_accel[i] += res[i]; // copy intermediate result
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    //for (unsigned char i = 0; i < 3; ++i)
    //  p.accel[i] = p.accel[i];

//    qDebug() << "p.dinv:" << p.dinv;
    //qDebug(" ac: % 10.5e % 10.5e % 10.5e", p.accel[0], p.accel[1], p.accel[2]);

  }
};

template<typename REAL, typename SIZE>
class CFPressureAccelRS :
  public CFQSpikyGradRS<REAL,SIZE,CFPressureAccelRS<REAL,SIZE> >
{
public:
  inline void init() { m_bound_kern.init(this->h); }

  inline void init_particle(FluidParticleR<REAL> &p)
  { Q_UNUSED(p); }
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;
#ifdef MCG03_PRESSURE_FORCE
    Vector3R<REAL> res(
        -p.dinv*this->m_mass*(p.pressure +
          near_p.pressure)*0.5*near_p.dinv*this->m_kern(p.pos - near_p.pos) // pressure contribution
        );
#else
    Vector3R<REAL> res(
        -this->m_mass*(p.pressure*p.dinv*p.dinv +
          near_p.pressure*near_p.dinv*near_p.dinv)*this->m_kern(p.pos - near_p.pos) // pressure contribution
        );
#endif

    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
   // Vector3R<REAL> res(
   //     -rest_density*near_p.dinv*p.pressure*p.dinv*p.dinv*p_kern[p.pos - near_p.pos] // pressure contribution
   //     );
#ifdef BOUNDARY_PENALTY_FORCE
    float massb = this->m_rest_density * near_p.dinv;
    Vector3R<REAL> res(
        this->m_cs2*(massb/(this->m_mass*(this->m_mass + massb))) * (p.pos - near_p.pos) * m_bound_kern(p.pos - near_p.pos)  // pressure contribution
        );

    for (unsigned char i = 0; i < 3; ++i)
      p.extern_accel[i] += res[i]; // copy intermediate result
#endif
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    p.extern_accel[1] -= M_G;
    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += p.extern_accel[i];

//    qDebug() << "p.dinv:" << p.dinv;
    //qDebug(" ac: % 10.5e % 10.5e % 10.5e", p.accel[0], p.accel[1], p.accel[2]);
  }
private:
  MKI04Kernel m_bound_kern;
};

#endif // QUANTITYPROCESSOR_H
