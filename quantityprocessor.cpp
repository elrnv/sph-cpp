#include "quantityprocessor.h"

// utility
template<typename REAL>
inline REAL pow7(REAL x) { return (x*x)*(x*x)*(x*x)*x; }

// Density

template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::init(REAL &mv, REAL &av)
{
  max_var = &mv; avg_var = &av;
}
template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::init_particle(FluidParticleRS<REAL,SIZE> &p)
{ 
  //fprintf(stderr, "density init_particle\n");
  p.dinv = 0.0f; p.vol = 0.0f;
}
template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::fluid(FluidParticleRS<REAL,SIZE> &p, FluidParticleRS<REAL,SIZE> &near_p)
{
  p.dinv += this->m_kern[ p.pos - near_p.pos ];
  p.vol += this->m_kern[ p.pos - near_p.pos ];
}
template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::bound(FluidParticleRS<REAL,SIZE> &p, ParticleR<REAL> &near_p)
{
  if (FT == MCG03)
    return;

  p.dinv += this->m_rest_density * near_p.dinv * this->m_kern[ p.pos - near_p.pos ];
  p.vol += this->m_kern[ p.pos - near_p.pos ];
}
template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::finish_particle(FluidParticleRS<REAL,SIZE> &p)
{
  //fprintf(stderr, "density finish_particle\n");
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

}


// Density Update

template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::init(float ts) 
{ 
  m_timestep = ts;
}
template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::init_particle(FluidParticleRS<REAL,SIZE> &p)
{
  p.temp = 0.0f;
}
template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::fluid(FluidParticleRS<REAL,SIZE> &p, FluidParticleRS<REAL,SIZE> &near_p)
{
  p.temp += (near_p.vel - p.vel).dot(this->m_kern[ p.pos - near_p.pos ]);
}
template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::bound(FluidParticleRS<REAL,SIZE> &p, ParticleR<REAL> &near_p)
{
  p.temp += (-p.vel).dot(this->m_kern[ p.pos - near_p.pos ]);
  //p.dinv += rest_density * near_p.dinv * kern[ p.pos - near_p.pos ];
}
template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::finish_particle(FluidParticleRS<REAL,SIZE> &p)
{
  p.dinv += m_timestep * p.dinv * p.dinv * p.temp * this->m_kern.coef;

  if (FT == MCG03)
    p.pressure = this->m_cs2 * (1.0f/p.dinv - this->m_rest_density);
  else // Enable tait pressure equation by default
    p.pressure =
      this->m_rest_density * this->m_cs2 * 0.14285714285714 * 
      (pow7(1.0f / (p.dinv * this->m_rest_density)) - 1);
}

// Acceleration

template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::init_kernel(float h)
{ }
template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::init()
{ 
}

template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::init_particle(FluidParticleRS<REAL,SIZE> &p)
{ 
  Q_UNUSED(p);
}

template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::fluid(FluidParticleRS<REAL,SIZE> &p, FluidParticleRS<REAL,SIZE> &near_p)
{
  Q_UNUSED(p);
  Q_UNUSED(near_p);
}

template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::bound(FluidParticleRS<REAL,SIZE> &p, ParticleR<REAL> &near_p)
{ 
  Q_UNUSED(p);
  Q_UNUSED(near_p);
}

template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::finish_particle(FluidParticleRS<REAL,SIZE> &p)
{ }


// Acceleration specialization 
template<typename REAL, typename SIZE>
class CFAccelRST<REAL,SIZE,MCG03> : public CFQ<REAL,SIZE,CFAccelRST<REAL,SIZE,MCG03> >
{
public:
  inline void init_kernel(float h)
  {
    m_spikygrad_kern.init(h);
    m_visclap_kern.init(h);
    m_colorgrad_kern.init(h);
  }
  inline void init_particle(FluidParticleRS<REAL,SIZE> &p)
  { Q_UNUSED(p); }

  inline void fluid(FluidParticleRS<REAL,SIZE> &p, FluidParticleRS<REAL,SIZE> &near_p)
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
  inline void bound(FluidParticleRS<REAL,SIZE> &p, ParticleR<REAL> &near_p) { }
  inline void finish_particle(FluidParticleRS<REAL,SIZE> &p)
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


template<typename REAL, typename SIZE>
class CFAccelRST<REAL,SIZE,BT07> : public CFQ<REAL,SIZE,CFAccelRST<REAL,SIZE,BT07> >
{
public:
  inline void init_kernel(float h)
  {
    m_poly6_kern.init(h);
    m_poly6grad_kern.init(h);
    m_spiky_kern.init(h);
    m_bound_kern.init(h);
  }
  inline void init_particle(FluidParticleRS<REAL,SIZE> &p)
  { Q_UNUSED(p); }

  inline void fluid(FluidParticleRS<REAL,SIZE> &p, FluidParticleRS<REAL,SIZE> &near_p)
  {
    if (&p == &near_p)
      return;

    // pressure

    Vector3R<REAL> res(
        -this->m_mass*(p.pressure*p.dinv*p.dinv +
          near_p.pressure*near_p.dinv*near_p.dinv)*this->m_spiky_kern(p.pos - near_p.pos));

    //qDebug(" ac: % 10.5e % 10.5e % 10.5e", res[0], res[1], res[2]);
    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result

    // viscosity

    Vector3R<REAL> x_ab = p.pos - near_p.pos;
    REAL vx = x_ab.dot(p.vel - near_p.vel);
    if (vx >= 0)
      return;

    REAL nu = 2*this->m_viscosity*this->m_radius;
    //REAL pab = -
    res = this->m_mass*this->m_poly6grad_kern(x_ab); // viscosity contribution

    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result

  }
  inline void bound(FluidParticleRS<REAL,SIZE> &p, ParticleR<REAL> &near_p)
  {
    // pressure contribution from boundary particles
    float massb = this->m_rest_density * near_p.dinv;
    Vector3R<REAL> res(
        this->m_cs2*(massb/(this->m_mass*(this->m_mass + massb))) * (p.pos - near_p.pos) * m_bound_kern(p.pos - near_p.pos)
        );

    for (unsigned char i = 0; i < 3; ++i)
      p.extern_accel[i] += res[i]; // copy intermediate result
  }
  inline void finish_particle(FluidParticleRS<REAL,SIZE> &p)
  {
    p.extern_accel[1] -= M_G;
    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += p.extern_accel[i];
  }
private:
  Poly6Kernel     m_poly6_kern;
  Poly6GradKernel m_poly6grad_kern;
  SpikyGradKernel m_spiky_kern;
  MKI04Kernel     m_bound_kern; // for boundary particles
}; // CFAccel


template class CFDensityRST<double, unsigned int, 0>;
template class CFDensityRST<double, unsigned int, 1>;
template class CFDensityRST<double, unsigned int, 2>;
template class CFDensityRST<double, unsigned int, 3>;
template class CFDensityUpdateRST<double, unsigned int, 0>;
template class CFDensityUpdateRST<double, unsigned int, 1>;
template class CFDensityUpdateRST<double, unsigned int, 2>;
template class CFDensityUpdateRST<double, unsigned int, 3>;
template class CFAccelRST<double, unsigned int, 0>;
template class CFAccelRST<double, unsigned int, 1>;
template class CFAccelRST<double, unsigned int, 2>;
template class CFAccelRST<double, unsigned int, 3>;
