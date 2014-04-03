#include "quantityprocessor.h"
#include "particle.h"
#include "settings.h"

// utility

template<typename REAL>
inline REAL pow7(REAL x) { return (x*x)*(x*x)*(x*x)*x; }


// Density and Pressure

template<typename REAL, typename SIZE>
class CFDensityRST<REAL,SIZE,MCG03> : public CFQ< REAL, SIZE, CFDensityRST<REAL,SIZE,MCG03> >
{
public:
  inline void init_kernel(float h)
  { 
    m_kern.init(h);
  }
  inline void init(REAL &mv, REAL &av)
  {
    max_var = &mv; avg_var = &av;
  }
  inline void init_particle(FluidParticleR<REAL> &p)
  { 
    p.dinv = 0.0f;
  }
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    p.dinv += this->m_kern[ p.pos - near_p.pos ];
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p) { }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    *p._dinv = p.dinv = 1.0f/(this->m_mass * p.dinv * this->m_kern.coef);

    REAL var = std::abs(1.0f/p.dinv - this->m_rest_density);
    if ( var > *max_var )
      *max_var = var;
    *avg_var += var;

    p.pressure = this->m_cs2 * (1.0f/p.dinv - this->m_rest_density);
  }

private:
  REAL *max_var; // max variation
  REAL *avg_var; // average variation
  Poly6Kernel m_kern;
};

template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::init_kernel(float h)
{ 
  m_kern.init(h);
}

template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::init(REAL &mv, REAL &av) 
{
  max_var = &mv; avg_var = &av;
}
template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{ 
  //fprintf(stderr, "density init_particle\n");
  p.dinv = 0.0f; p.vol = 0.0f;
}
template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{
  p.dinv += this->m_kern[ p.pos - near_p.pos ];
  p.vol += this->m_kern[ p.pos - near_p.pos ];
 // fprintf(stderr, "p.dinv += %.2e\n", p.dinv);
}
template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
{
  if (FT == BT07)
    return;

  p.dinv += this->m_rest_density * near_p.dinv * this->m_kern[ p.pos - near_p.pos ];
  p.vol += this->m_kern[ p.pos - near_p.pos ];
}
template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p)
{
#if 1
    *p._dinv = p.dinv = 1.0f/(this->m_mass * p.dinv * this->m_kern.coef);
#else
#if 0
  qDebug() << (this->m_mass * p.dinv) << "/"
    << (8*this->m_radius*this->m_radius*this->m_radius*p.vol) << " = " <<
    (this->m_mass * p.dinv)/(8*this->m_radius*this->m_radius*this->m_radius*p.vol);
#endif
    // kernel correction
    *p._dinv = p.dinv =
      (8*this->m_radius*this->m_radius*this->m_radius*p.vol)/(this->m_mass * p.dinv);
#endif

  REAL var = std::abs(1.0f/p.dinv - this->m_rest_density);
  if ( var > *max_var )
    *max_var = var;
  *avg_var += var;

#if 1
    p.pressure = this->m_cs2 * (1.0f/p.dinv - this->m_rest_density);
#else
  p.pressure =
    this->m_rest_density * this->m_cs2 * 0.14285714285714 * 
    (pow7(1.0f / (p.dinv * this->m_rest_density)) - 1);
#endif
}



// Density Update

template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::init_kernel(float h)
{ 
  m_kern.init(h);
}

template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::init(float ts) 
{
  m_timestep = ts;
}
template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{ 
  p.c = 0.0f; // temp value
  //qDebug() << "init dinv" << *p._dinv;
}
template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{
  if (&p == &near_p)
    return;
  p.c += (p.vel - near_p.vel).dot(this->m_kern[p.pos - near_p.pos]);
}
template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
{
  p.c += (p.vel).dot(this->m_kern[p.pos - near_p.pos]);
}
template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p)
{
  //qDebug() << "kern: " << this->m_kern.coef;
  //qDebug() << "time: " << m_timestep;
  //qDebug() << "p.c: " << p.c;
  //qDebug() << "d increment: " << m_timestep * p.c * this->m_mass * this->m_kern.coef;
  //qDebug() << "dinv inc: " << m_timestep * (*p._dinv) * (*p._dinv) * p.c * this->m_mass * this->m_kern.coef;

  qDebug() << (this->m_mass * p.dinv) << "/"
    << (8*this->m_radius*this->m_radius*this->m_radius*p.vol) << " = " <<
    (this->m_mass * p.dinv)/(8*this->m_radius*this->m_radius*this->m_radius*p.vol);
  //qDebug() << "new dinv with ddinv: " << (*p._dinv) - m_timestep * (*p._dinv)
    //* (*p._dinv) * p.c * this->m_mass * this->m_kern.coef;
    *p._dinv = p.dinv = 1.0f/(1.0f/(*p._dinv) + m_timestep * p.c * this->m_mass * this->m_kern.coef);
    //p.dinv -= m_timestep * p.dinv * p.dinv * p.c * this->m_mass * this->m_kern.coef;

  //qDebug() << "new dinv with d: " << *p._dinv;

#if 1
    p.pressure = this->m_cs2 * (1.0f/p.dinv - this->m_rest_density);
#else
    p.pressure =
      this->m_rest_density * this->m_cs2 * 0.14285714285714 * 
      (pow7(1.0f / (p.dinv * this->m_rest_density)) - 1);
#endif
}




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
    m_colorlap_kern.init(h);
  }
  inline void init_particle(FluidParticleR<REAL> &p)
  { 
    p.c = 0.0f;
    p.n[0] = p.n[1] = p.n[2] = 0.0f;
  }

  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;
    
    Vector3R<REAL> res(
    // pressure
         -0.5 * near_p.dinv * (p.pressure + near_p.pressure) 
              * m_spikygrad_kern(p.pos - near_p.pos)
           +
    // viscosity
          this->m_viscosity * near_p.dinv * (near_p.vel - p.vel) 
              * m_visclap_kern(p.pos - near_p.pos)
    );

    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result

    // surface tension
    p.n = p.n + near_p.dinv * m_colorgrad_kern(p.pos - near_p.pos);
    p.c += near_p.dinv * m_colorlap_kern(p.pos - near_p.pos);
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p) { }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    REAL nnorm2 = p.n.squaredNorm();
    Vector3R<REAL> st(0.0f,0.0f,0.0f);
    if (nnorm2 > 0.01)
    {
      st = -(this->m_st / std::sqrt(nnorm2)) * p.c * p.n;
    }

    p.extern_accel[0] += global::dynset.gravity[0];
    p.extern_accel[1] += global::dynset.gravity[1];
    p.extern_accel[2] += global::dynset.gravity[2];
    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] = (this->m_mass*p.dinv*(p.accel[i] + st[i])) + p.extern_accel[i];

    //qDebug() << "p.accel = " << p.accel[0] << p.accel[1] << p.accel[2];
  }

private:
  SpikyGradKernel m_spikygrad_kern; // for pressure
  Poly6GradKernel m_colorgrad_kern;    // for surface tension
  Poly6LapKernel  m_colorlap_kern;     // for surface tension
  ViscLapKernel   m_visclap_kern;   // for viscosity
};

template<typename REAL, typename SIZE>
class CFAccelRST<REAL,SIZE,BT07> : public CFQ<REAL,SIZE,CFAccelRST<REAL,SIZE,BT07> >
{
public:
  inline void init_kernel(float h)
  {
    m_splinegrad_kern.init(h);
    m_spikygrad_kern.init(h);
    m_bound_kern.init(h);
    m_poly6_kern.init(h);
  }
  inline void init_particle(FluidParticleR<REAL> &p)
  {
    p.n[0] = p.n[1] = p.n[2] = 0.0f;
  }

  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;

    // pressure

    REAL res((p.pressure*p.dinv*p.dinv + near_p.pressure*near_p.dinv*near_p.dinv));

    // viscosity

    Vector3R<REAL> x_ab = p.pos - near_p.pos;
    REAL vx = x_ab.dot(p.vel - near_p.vel);
    if (vx < 0) // compute damping term
    {
      REAL nu = 2*this->m_viscosity*this->m_radius*this->m_cs 
              * p.dinv * near_p.dinv / (p.dinv + near_p.dinv);
      res -= nu*vx / (x_ab.squaredNorm() + 0.01*this->m_radius*this->m_radius);
    }

    p.n = p.n - res * this->m_splinegrad_kern(x_ab);
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    // pressure contribution from boundary particles
    float massb = this->m_rest_density * near_p.dinv;
    Vector3R<REAL> x_ab = p.pos - near_p.pos;
    REAL kern = m_bound_kern((1.0f/this->m_radius)*x_ab);
    Vector3R<REAL> res(
        (this->m_cs2*massb/(this->m_mass + massb)) 
        * x_ab * m_bound_kern(x_ab));

    p.n = p.n + res;
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    p.extern_accel[0] += global::dynset.gravity[0];
    p.extern_accel[1] += global::dynset.gravity[1];
    p.extern_accel[2] += global::dynset.gravity[2];
    for (unsigned char i = 0; i < 3; ++i)
    {
      p.extern_accel[i] += this->m_mass * p.n[i];
      p.accel[i] += p.extern_accel[i];
    }
  }
private:
  SpikyGradKernel m_spikygrad_kern;
  Poly6GradKernel m_splinegrad_kern;
  MKI04Kernel     m_bound_kern; // for boundary particles
  Poly6Kernel     m_poly6_kern; // for boundary particles
}; // CFAccel

template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::init_kernel(float h) { }
template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p) { }
template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p) { }
template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p) { }
template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p) { }

template class CFAccelRST<double, unsigned int, 0>;
template class CFAccelRST<double, unsigned int, 1>;
template class CFAccelRST<double, unsigned int, 2>;
template class CFAccelRST<double, unsigned int, 3>;
template class CFDensityRST<double, unsigned int, 0>;
template class CFDensityRST<double, unsigned int, 1>;
template class CFDensityRST<double, unsigned int, 2>;
template class CFDensityRST<double, unsigned int, 3>;
template class CFDensityUpdateRST<double, unsigned int, 0>;
template class CFDensityUpdateRST<double, unsigned int, 1>;
template class CFDensityUpdateRST<double, unsigned int, 2>;
template class CFDensityUpdateRST<double, unsigned int, 3>;
