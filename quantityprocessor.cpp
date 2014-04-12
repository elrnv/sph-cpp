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
  //if (FT == BT07)
  //  return;

  //p.dinv += this->m_rest_density * near_p.dinv * this->m_kern[ p.pos - near_p.pos ];
  //p.vol += this->m_kern[ p.pos - near_p.pos ];
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

#if 0
    p.pressure = this->m_cs2 * (1.0f/p.dinv - this->m_rest_density);
#else
  p.pressure =
    this->m_rest_density * this->m_cs2 * 0.14285714285714 * 
    (pow7(1.0f / (p.dinv * this->m_rest_density)) - 1.0f);
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
    
    Vector3R<REAL> v_ab(p.vel - near_p.vel);
    Vector3R<REAL> x_ab(p.pos - near_p.pos);

    Vector3R<REAL> res(
    // pressure
         0.5 * near_p.dinv * (p.pressure + near_p.pressure) 
              * m_spikygrad_kern(x_ab)
           +
    // viscosity
//          (v_ab.dot(x_ab) < 0 ? 
          this->m_viscosity * near_p.dinv * v_ab
              * m_visclap_kern(x_ab) 
             // : Vector3R<REAL>(0.0,0.0,0.0))
    );

    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] -= res[i]; // copy intermediate result

    // surface tension
    p.n = p.n + near_p.dinv * m_colorgrad_kern(p.pos - near_p.pos);
    p.c += near_p.dinv * m_colorlap_kern(p.pos - near_p.pos);
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p) { }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    REAL nnorm = p.n.norm();
    Vector3R<REAL> st(0.0f,0.0f,0.0f);
    if (nnorm > 0.05)
      st = -(this->m_st / nnorm) * p.c * p.n;

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
    m_maingrad_kern.init(h);
    m_spikygrad_kern.init(h);
    m_bound_kern.init(h);
    m_st_kern.init(h);
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

    REAL res(p.pressure*p.dinv*p.dinv + near_p.pressure*near_p.dinv*near_p.dinv);

    // viscosity

    Vector3R<REAL> x_ab = p.pos - near_p.pos;
    REAL vx = x_ab.dot(p.vel - near_p.vel);
    if (vx < 0)
    {
      REAL nu = 2*this->m_viscosity*m_maingrad_kern.h*this->m_cs 
              * p.dinv * near_p.dinv / (p.dinv + near_p.dinv);
      res -= nu*vx / (x_ab.squaredNorm() + 0.01*m_maingrad_kern.h2);
    }

    p.n = p.n - res * m_maingrad_kern(x_ab);
    
    // surface tension 
    // (apply surface tension only to particles from the same phase)
    if (p.id == near_p.id)
      p.n = p.n - this->m_st * m_st_kern(x_ab) * x_ab;
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    // BT07 fluids dont interact with static boundary particles, instead they
    // mirror fluid particles on the boundary, and compute repulsion forces
    // this way (in the finish_particle method)
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    // if the particle is close to the boundary and add a repulsive force
    // (assume boundary particle is infinitely more massive than the fluid
    // particle, such that m_b / (m_a + m_b) == 1 where m_a and m_b are the
    // masses of the fluid and boundary particles respectively
    Vector3R<REAL> hvec(m_bound_kern.h,m_bound_kern.h,m_bound_kern.h);
    Vector3R<REAL> res(0.0f, 0.0f, 0.0f);
    Vector3R<REAL> dmin =  // distances to min boundaries
       hvec + p.pos - this->m_bmin.template cast<REAL>();
    Vector3R<REAL> dmax =  // distances to max boundaries
      hvec + this->m_bmax.template cast<REAL>() - p.pos;
    if ( dmin[0] < m_bound_kern.h )
      res[0] += m_bound_kern(Vector3d(dmin[0],0,0));// / dmin[0];
    else if ( dmax[0] < m_bound_kern.h)
      res[0] -= m_bound_kern(Vector3d(dmax[0],0,0));// / dmax[0];

    if ( dmin[1] < m_bound_kern.h )
      res[1] += m_bound_kern(Vector3d(dmin[1],0,0));// / dmin[1];
    else if ( dmax[1] < m_bound_kern.h)
      res[1] -= m_bound_kern(Vector3d(dmax[1],0,0));// / dmax[1];

    if ( dmin[2] < m_bound_kern.h )
      res[2] += m_bound_kern(Vector3d(dmin[2],0,0));// / dmin[2];
    else if ( dmax[2] < m_bound_kern.h)
      res[2] -= m_bound_kern(Vector3d(dmax[2],0,0));// / dmax[2];

    p.n = p.n + this->m_compressibility*this->m_cs2*res;

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
  CubicSplineGradKernel m_maingrad_kern;
  MKI04Kernel     m_bound_kern;
  CubicSplineKernel m_st_kern;
}; // CFAccel


template<typename REAL, typename SIZE>
class CFAccelRST<REAL,SIZE,ICS13> : public CFQ<REAL,SIZE,CFAccelRST<REAL,SIZE,ICS13> >
{
public:
  inline void init_kernel(float h)
  {
    m_maingrad_kern.init(h);
    m_spikygrad_kern.init(h);
    m_bound_kern.init(h);
    m_st_kern.init(h);
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

    REAL res(p.pressure*p.dinv*p.dinv + near_p.pressure*near_p.dinv*near_p.dinv);

    // viscosity

    Vector3R<REAL> x_ab = p.pos - near_p.pos;
    REAL vx = x_ab.dot(p.vel - near_p.vel);
    if (vx < 0)
    {
      REAL nu = 2*this->m_viscosity*m_maingrad_kern.h*this->m_cs 
              * p.dinv * near_p.dinv / (p.dinv + near_p.dinv);
      res -= nu*vx / (x_ab.squaredNorm() + 0.01*m_maingrad_kern.h2);
    }

    p.n = p.n - res * m_maingrad_kern(x_ab);
    
    // surface tension 
    // (apply surface tension only to particles from the same phase)
    if (p.id == near_p.id)
      p.n = p.n - this->m_st * m_st_kern(x_ab) * x_ab;
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    // pressure contribution from boundary particles
    //float massb = this->m_rest_density * near_p.dinv;
    //Vector3R<REAL> x_ab = p.pos - near_p.pos;
    //Vector3d kern = m_bound_kern(x_ab);

    //if (0 && x_ab.norm() < 1.5f*m_bound_kern.h)
    //{
    //  qDebug() << "xab = "  << x_ab[0] << x_ab[1] << x_ab[2];
    //  qDebug() << "xab.norm = "  << x_ab.norm();
    //  qDebug() << "h = "  << m_bound_kern.h;
    //  qDebug() << "massb = "  << massb;
    //  qDebug() << "kern = " << kern[0] << kern[1] << kern[2];
    //}
    //Vector3R<REAL> res(
    //    (this->m_cs2*massb/(this->m_mass + massb)) * kern / x_ab.norm());

    //p.n = p.n + res;
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
  CubicSplineGradKernel m_maingrad_kern;
  MKI04Kernel     m_bound_kern;
  CubicSplineKernel m_st_kern;
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
