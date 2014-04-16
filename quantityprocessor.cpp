#include <cmath>
#include "quantityprocessor.h"
#include "particle.h"
#include "settings.h"

#define BT07_BOUNDARY_PARTICLES

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
    *p._dinv = p.dinv = 1.0f / (this->m_mass * p.dinv * this->m_kern.coef);
    p.vol = p.dinv * this->m_mass;

    p.pressure =
      std::max(REAL(0), this->m_cs2 * (1.0f/p.dinv - this->m_rest_density));

    REAL var = std::abs(1.0f/p.dinv - this->m_rest_density);
    if ( var > *max_var )
      *max_var = var;
    *avg_var += var;
  }

private:
  REAL *max_var; // max variation
  REAL *avg_var; // average variation
  Poly6Kernel m_kern;
};

template<typename REAL, typename SIZE>
class CFDensityRST<REAL,SIZE,BT07> : public CFQ< REAL, SIZE, CFDensityRST<REAL,SIZE,BT07> >
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
    p.vol = 0.0f;
  }
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    p.dinv += this->m_mass * this->m_kern[ p.pos - near_p.pos ];
    p.vol += this->m_kern[ p.pos - near_p.pos ];
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p) 
  {
    // near_p.dinv is inverse number density (aka volume) of the boundary paticle
#ifdef BT07_BOUNDARY_PARTICLES
    p.dinv += this->m_rest_density * near_p.dinv * this->m_kern[ p.pos - near_p.pos ];
    p.vol += this->m_kern[ p.pos - near_p.pos ];
#endif
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    *p._dinv = p.dinv = 1.0f/(p.dinv * this->m_kern.coef);
    p.vol = 1.0f / (p.vol * this->m_kern.coef);

    p.pressure = std::max(REAL(0),
      this->m_rest_density * this->m_cs2 * 0.14285714285714 * 
      (pow7(1.0f / (p.dinv * this->m_rest_density)) - 1.0f));

    REAL var = std::abs(1.0f/p.dinv - this->m_rest_density);
    if ( var > *max_var )
      *max_var = var;
    *avg_var += var;
  }

private:
  REAL *max_var; // max variation
  REAL *avg_var; // average variation
  Poly6Kernel m_kern;
};


template<typename REAL, typename SIZE, int FT>
inline void CFDensityRST<REAL,SIZE,FT>::init_kernel(float h)
{ 
  m_kern.init(h);
}

template<typename REAL, typename SIZE, int FT>
inline void CFDensityRST<REAL,SIZE,FT>::init(REAL &mv, REAL &av) 
{
  max_var = &mv; avg_var = &av;
}
template<typename REAL, typename SIZE, int FT>
inline void CFDensityRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{ 
  p.dinv = 0.0f; p.vol = 0.0f;
}
template<typename REAL, typename SIZE, int FT>
inline void CFDensityRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{
  p.dinv += this->m_mass * this->m_kern[ p.pos - near_p.pos ];
  p.vol += this->m_kern[ p.pos - near_p.pos ];
}
template<typename REAL, typename SIZE, int FT>
inline void CFDensityRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
{
  // near_p.dinv is inverse number density (aka volume) of the boundary paticle
  //p.dinv += this->m_rest_density * near_p.dinv * this->m_kern[ p.pos - near_p.pos ];
  //p.vol += this->m_kern[ p.pos - near_p.pos ];
}
template<typename REAL, typename SIZE, int FT>
inline void CFDensityRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p)
{
  *p._dinv = p.dinv = 1.0f/(p.dinv * this->m_kern.coef);
  p.vol = 1.0f / (p.vol * this->m_kern.coef);

  p.pressure = 0.0f;

  REAL var = std::abs(1.0f/p.dinv - this->m_rest_density);
  if ( var > *max_var )
    *max_var = var;
  *avg_var += var;
}



// Density Update

template<typename REAL, typename SIZE, int FT>
inline void CFDensityUpdateRST<REAL,SIZE,FT>::init_kernel(float h)
{ 
  m_kern.init(h);
}
template<typename REAL, typename SIZE, int FT>
inline void CFDensityUpdateRST<REAL,SIZE,FT>::init(float ts) 
{
  m_timestep = ts;
}
template<typename REAL, typename SIZE, int FT>
inline void CFDensityUpdateRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{ 
  p.c = 0.0f; // temp value
}
template<typename REAL, typename SIZE, int FT>
inline void CFDensityUpdateRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{
  if (&p == &near_p)
    return;
  p.c += (p.vel - near_p.vel).dot(this->m_kern[p.pos - near_p.pos]);
}
template<typename REAL, typename SIZE, int FT>
inline void CFDensityUpdateRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
{
  p.c += (p.vel).dot(this->m_kern[p.pos - near_p.pos]);
}
template<typename REAL, typename SIZE, int FT>
inline void CFDensityUpdateRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p)
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
    m_color_kern.init(h);
    m_colorgrad_kern.init(h);
  }
  inline void init_particle(FluidParticleR<REAL> &p)
  { 
    p.c = p.vol * REAL(p.color) * m_color_kern(Vector3R<REAL>(0.0f,0.0f,0.0f));
    p.n[0] = p.n[1] = p.n[2] = 0.0f;
    p.m[0] = p.vol * m_color_kern(Vector3R<REAL>(0.0f,0.0f,0.0f));
  }

  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;
    
    Vector3R<REAL> v_ab(p.vel - near_p.vel);
    Vector3R<REAL> x_ab(p.pos - near_p.pos);

    Vector3R<REAL> res(
    // pressure
        near_p.vol * 0.5*(p.pressure + near_p.pressure)*m_spikygrad_kern(x_ab)
        +
    // viscosity
        this->m_viscosity * near_p.vol * v_ab * m_visclap_kern(x_ab) 
    );

    for (unsigned char i = 0; i < 3; ++i)
      p._accel[i] -= res[i]; // copy intermediate result

    // surface tension based on color field (intermediate values)
    p.c += near_p.vol * REAL(near_p.color) * m_color_kern(x_ab);
    p.m[0] += near_p.vol * m_color_kern(x_ab); // denom for smoothed color
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p) { }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    //qDebug() << "p.c = " << p.c;
    p.c = p.c / p.m[0];

    //qDebug() << "p.c/normal = " << p.c;

    p._extern_accel[0] += global::dynset.gravity[0];
    p._extern_accel[1] += global::dynset.gravity[1];
    p._extern_accel[2] += global::dynset.gravity[2];
    for (unsigned char i = 0; i < 3; ++i)
      p._accel[i] = ((p.vol/this->m_mass)*p._accel[i]) + p._extern_accel[i];

    //qDebug() << "p._accel = " << p._accel[0] << p._accel[1] << p._accel[2];
  }

private:
  SpikyGradKernel m_spikygrad_kern; // for pressure
  ViscLapKernel   m_visclap_kern;   // for viscosity
  CubicSplineKernel m_color_kern;   // for surface tension
  Poly6GradKernel   m_colorgrad_kern; // for surface tension
};



// Computing the surface normal for color based surface tension
template<typename REAL, typename SIZE, int FT>
inline void CFSurfaceNormalRST<REAL,SIZE,FT>::init_kernel(float h)
{ m_grad_kern.init(h); }
template<typename REAL, typename SIZE, int FT>
inline void CFSurfaceNormalRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{
  p.n[0] = p.n[1] = p.n[2] = 0.0f;
}
template<typename REAL, typename SIZE, int FT>
inline void CFSurfaceNormalRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{ 
  if (&p == &near_p)
    return;
  p.n = p.n + near_p.vol * (near_p.c - p.c) * m_grad_kern(p.pos - near_p.pos);
}
template<typename REAL, typename SIZE, int FT>
inline void CFSurfaceNormalRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
{ /* assume non particle boundaries */ }
template<typename REAL, typename SIZE, int FT>
inline void CFSurfaceNormalRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p)
{  
  // store normal norms
  REAL norm = p.n.norm();
  p.m[1] = norm > (0.01/m_grad_kern.h) ? 1.0f/norm : 0.0f;
}


// surface tension from curvature based on color values
template<typename REAL, typename SIZE, int FT>
inline void CFSurfaceTensionRST<REAL,SIZE,FT>::init_kernel(float h)
{ 
  m_grad_kern.init(h);
  m_kern.init(h);
}
template<typename REAL, typename SIZE, int FT>
inline void CFSurfaceTensionRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{
  p.m[2] = 0.0f;
  p.m[0] = 0.0f; // recompute denom
}
template<typename REAL, typename SIZE, int FT>
inline void CFSurfaceTensionRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{ 
  // compute curvature
  if (p.m[1] && near_p.m[1]) // if normals are in good shape
  {
    Vector3R<REAL> x_ab = p.pos - near_p.pos;
    p.m[2] -= near_p.vol * 
      (near_p.n*near_p.m[1] - p.n*p.m[1]).dot(m_grad_kern(x_ab));
    p.m[0] += near_p.vol * m_kern(x_ab);     // denom for color
  }
}
template<typename REAL, typename SIZE, int FT>
inline void CFSurfaceTensionRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
{ }
template<typename REAL, typename SIZE, int FT>
inline void CFSurfaceTensionRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p)
{ 
  //qDebug() << "pvol = " << p.vol;
  //qDebug() << "st = " << this->m_st;
  //qDebug() << "p.n = " << p.n[0] << p.n[1] << p.n[2];
  //qDebug() << "k = " << p.m[2]/p.m[0];
  if (p.m[0] && p.m[2])
    for (unsigned char i = 0; i < 3; ++i)
      p._accel[i] += p.vol*this->m_st * p.n[i] * (p.m[2]/p.m[0]);
}



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
    p.m[0] = p.m[1] = p.m[2] = 0.0f;
  }

  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;

    // pressure

    REAL res = (p.pressure*p.vol*p.vol + near_p.pressure*near_p.vol*near_p.vol)/this->m_mass;

    // viscosity

    Vector3R<REAL> x_ab = p.pos - near_p.pos;
    REAL vx = x_ab.dot(p.vel - near_p.vel);
    if (vx < 0)
    {
      REAL nu = 2*this->m_viscosity*m_maingrad_kern.h*this->m_cs *
              near_p.vol * p.dinv / (p.dinv + near_p.dinv);
      res -= nu*vx / (x_ab.squaredNorm() + 0.01*m_maingrad_kern.h2);
    }

    p.n = p.n - res * m_maingrad_kern(x_ab);
    
    // surface tension 
    // (apply surface tension only to particles from the same phase)
    if (p.color == near_p.color)
      p.n = p.n - (near_p.vol/(near_p.dinv*this->m_mass))*this->m_st * m_st_kern(x_ab) * x_ab;
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    // BT07 fluids dont interact with static boundary particles, instead they
    // mirror fluid particles on the boundary, and compute repulsion forces
    // this way (in the finish_particle method)

#ifdef BT07_BOUNDARY_PARTICLES // the following handles boundaries as in AIAST12
    Vector3R<REAL> x_ab = p.pos - near_p.pos;
    float massb = this->m_rest_density * near_p.dinv;

    // pressure force contribution from boundary particles
    REAL res( p.pressure * p.dinv * p.dinv);

    REAL vx = x_ab.dot(p.vel);
    // viscosity foces contribution from boundary particles
    
    if (vx < 0)
    {
      REAL nu = 0.5*this->m_friction*m_maingrad_kern.h*this->m_cs*p.dinv;
      res -= nu*vx / (x_ab.squaredNorm() + 0.01*m_maingrad_kern.h2);
    }

    p.m = p.m - massb * res * m_maingrad_kern(x_ab);
#endif
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    // if the particle is close to the boundary and add a repulsive force
    // (assume boundary particle is infinitely more massive than the fluid
    // particle, such that m_b / (m_a + m_b) == 1 where m_a and m_b are the
    // masses of the fluid and boundary particles respectively
#ifndef BT07_BOUNDARY_PARTICLES
    Vector3R<REAL> hvec(m_bound_kern.h,m_bound_kern.h,m_bound_kern.h);
    Vector3R<REAL> res(0.0f, 0.0f, 0.0f);
    Vector3R<REAL> dmin =  // distances to min boundaries
      0.5*hvec + p.pos - this->m_bmin.template cast<REAL>();
    Vector3R<REAL> dmax =  // distances to max boundaries
      0.5*hvec + this->m_bmax.template cast<REAL>() - p.pos;
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

    // acceleration due to boundary collisions
    Vector3R<REAL> bdry_accel( this->m_compressibility*this->m_cs2*res );
#endif
    for (unsigned char i = 0; i < 3; ++i)
    {
      p._extern_accel[i] += global::dynset.gravity[i] +
#ifdef BT07_BOUNDARY_PARTICLES
        p.m[i];
#else
        bdry_accel[i];
#endif
      p._accel[i] += p._extern_accel[i] + p.n[i];
    }
  }
private:
  SpikyGradKernel       m_spikygrad_kern;
  CubicSplineGradKernel m_maingrad_kern;
  MKI04Kernel           m_bound_kern;
  CubicSplineKernel     m_st_kern;
}; // CFAccel


template<typename REAL, typename SIZE>
class CFAccelRST<REAL,SIZE,ICS13> : public CFQ<REAL,SIZE,CFAccelRST<REAL,SIZE,ICS13> >
{
public:
  inline void init_kernel(float h)
  {
    m_grad_kern.init(h);
    m_st_kern.init(h);
  }
  inline void init_particle(FluidParticleR<REAL> &p)
  {
    p.n[0] = p.n[1] = p.n[2] = 0.0f;
    p.m[0] = p.m[1] = p.m[2] = 0.0f; // d_ii
  }

  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  { // non-pressure accelerations
    if (&p == &near_p)
      return;

    // viscosity

    REAL res = 0.0f;
    Vector3R<REAL> x_ab = p.pos - near_p.pos;
    REAL vx = x_ab.dot(p.vel - near_p.vel);
    if (vx < 0)
    {
      REAL nu = 2*this->m_viscosity*m_grad_kern.h*this->m_cs 
              * near_p.vol * p.dinv / (p.dinv + near_p.dinv);
      res -= nu*vx / (x_ab.squaredNorm() + 0.01*m_grad_kern.h2);
    }

    //p.n = p.n - res * m_grad_kern(x_ab);
    
    // surface tension 
    // (apply surface tension only to particles from the same phase)
    //if (p.color == near_p.color)
    //  p.n = p.n - (near_p.vol/(near_p.dinv*this->m_mass))*this->m_st * m_st_kern(x_ab) * x_ab;

    p.m = p.m - (p.vol*p.vol/this->m_mass)*m_grad_kern(x_ab); // d_ii
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    //p.m = p.m - (this->m_rest_density*near_p.dinv*p.dinv*p.dinv)*m_grad_kern(p.pos - near_p.pos); // d_ii
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    //p.m = p.m * (1.0f/this->m_mass); // d_ii

    for (unsigned char i = 0; i < 3; ++i)
    {
      p._extern_accel[i] += global::dynset.gravity[i];
      p._accel[i] += p._extern_accel[i] + p.n[i];
    }
  }
private:
  CubicSplineGradKernel m_grad_kern;
  CubicSplineKernel m_st_kern;
  SpikyGradKernel m_spikygrad_kern;
  Poly6Kernel m_kern;
}; // CFAccel

// default implementation
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


template<typename REAL, typename SIZE, int FT>
inline void CFPrepareJacobiRST<REAL,SIZE,FT>::init_kernel(float h)
{ 
  m_grad_kern.init(h);
}
template<typename REAL, typename SIZE, int FT>
inline void CFPrepareJacobiRST<REAL,SIZE,FT>::init(float ts)
{ dt = ts; }
template<typename REAL, typename SIZE, int FT>
inline void CFPrepareJacobiRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{
  p.b = 0.0f; // intermediate value for density guess
  p.c = 0.0f; // a_ii
}
template<typename REAL, typename SIZE, int FT>
inline void CFPrepareJacobiRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{ 
  if (&p == &near_p)
    return;

  Vector3R<REAL> x_ab(p.pos - near_p.pos);
  Vector3R<REAL> v_ab(p.vel - near_p.vel);
  p.b += p.dinv*near_p.vol*v_ab.dot(m_grad_kern(x_ab)); // density guess
  p.c += near_p.vol*((-dt*dt*p.vol*p.vol*near_p.dinv/near_p.vol)*m_grad_kern(-x_ab) 
      - dt*dt*p.m).dot(m_grad_kern(x_ab));
}
template<typename REAL, typename SIZE, int FT>
inline void CFPrepareJacobiRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
{
//  Vector3R<REAL> x_ab(p.pos - near_p.pos);
//  p.b += p.dinv*p.dinv*this->m_rest_density*near_p.dinv*p.vel.dot(m_grad_kern(x_ab)); // density guess
}
template<typename REAL, typename SIZE, int FT>
inline void CFPrepareJacobiRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p)
{  
  p.c *= p.dinv; // a_ii
  p.dinv = p.dinv - dt*p.b; // density guess
  p.pressure = 0.5*p.pressure;     // initial pressure guess
}


template<typename REAL, typename SIZE, int FT>
inline void CFJacobiSolveFirstRST<REAL,SIZE,FT>::init_kernel(float h)
{ m_grad_kern.init(h); }
template<typename REAL, typename SIZE, int FT>
inline void CFJacobiSolveFirstRST<REAL,SIZE,FT>::init(float ts)
{ dt = ts; }
template<typename REAL, typename SIZE, int FT>
inline void CFJacobiSolveFirstRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{
  p.n[0] = p.n[1] = p.n[2] = 0.0f; // sum of d_ij
}
template<typename REAL, typename SIZE, int FT>
inline void CFJacobiSolveFirstRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{ 
  if (&p == &near_p)
    return;
  
  Vector3R<REAL> x_ab(p.pos - near_p.pos);
  p.n = p.n - (near_p.pressure*near_p.vol*near_p.vol)*m_grad_kern(x_ab);
}
template<typename REAL, typename SIZE, int FT>
inline void CFJacobiSolveFirstRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
{ }
template<typename REAL, typename SIZE, int FT>
inline void CFJacobiSolveFirstRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p)
{  
  p.n = p.n * dt*dt/this->m_mass;
}


template<typename REAL, typename SIZE, int FT>
inline void CFJacobiSolveSecondRST<REAL,SIZE,FT>::init_kernel(float h)
{ m_grad_kern.init(h); }
template<typename REAL, typename SIZE, int FT>
inline void CFJacobiSolveSecondRST<REAL,SIZE,FT>::init(float ts, REAL &avg_d, REAL &avg_p)
{ dt = ts; avg_density = &avg_d; avg_pressure = &avg_p;  }
template<typename REAL, typename SIZE, int FT>
inline void CFJacobiSolveSecondRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{
  p.b = 0.0f;
}
template<typename REAL, typename SIZE, int FT>
inline void CFJacobiSolveSecondRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{ 
  if (&p == &near_p)
    return;
  
  Vector3R<REAL> x_ab(p.pos - near_p.pos);
  p.b += near_p.vol*(near_p.n - p.n + dt*dt*near_p.m*near_p.pressure).dot(m_grad_kern(x_ab));
}
template<typename REAL, typename SIZE, int FT>
inline void CFJacobiSolveSecondRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
{
 // Vector3R<REAL> x_ab(p.pos - near_p.pos);
 // p.b -= (p.vol/this->m_mass)*(this->m_rest_density*near_p.vol*p.n).dot(m_grad_kern(x_ab));
}
template<typename REAL, typename SIZE, int FT>
inline void CFJacobiSolveSecondRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p)
{  
  REAL A = p.b*p.vol/this->m_mass;
  *avg_density += 1.0f/(p.dinv + p.pressure*p.c + A);
 // REAL avg = 1.0f/p.dinv + p.pressure*p.c + A;

  if (p.pos[1] < -0.90)
  {
    //qDebug() << "p = " << &p;
    //qDebug() << "A = " << A;
    //qDebug() << "a_ii = " << p.c;
    //qDebug() << "cur pressure = " << p.pressure;
    //qDebug() << "avg_Density = " << avg;
//    qDebug() << "suggestsed pressure = " << (1.0f/this->m_rest_density - p.dinv - A)/p.c;
  }
  *avg_pressure += (1.0f/this->m_rest_density - p.dinv - A)/p.c;
  
  p.pressure = //std::max(REAL(0), // clamp negative pressures
      0.5*(p.pressure + (1.0f/this->m_rest_density - p.dinv - A)/p.c);//);
}

// update pressure accelerations
template<typename REAL, typename SIZE, int FT>
inline void CFPressureAccelRST<REAL,SIZE,FT>::init_kernel(float h)
{ m_grad_kern.init(h); }
template<typename REAL, typename SIZE, int FT>
inline void CFPressureAccelRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{  
  p.n[0] = p.n[1] = p.n[2] = 0.0f;
  //p.m[0] = p.m[1] = p.m[2] = 0.0f;
}
template<typename REAL, typename SIZE, int FT>
inline void CFPressureAccelRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{ 
  if (&p == &near_p)
    return;
  
  Vector3R<REAL> x_ab(p.pos - near_p.pos);
  p.n = p.n - (p.pressure*p.vol*p.vol + near_p.pressure*near_p.vol*near_p.vol) * m_grad_kern(x_ab);
}
template<typename REAL, typename SIZE, int FT>
inline void CFPressureAccelRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
{
  Vector3R<REAL> x_ab = p.pos - near_p.pos;
  float massb = this->m_rest_density * near_p.dinv;

  // pressure force contribution from boundary particles
  REAL res( p.pressure * p.dinv * p.dinv);

  REAL vx = x_ab.dot(p.vel);
  // viscosity foces contribution from boundary particles

  if (vx < 0)
  {
    REAL nu = 0.5*this->m_friction*m_grad_kern.h*this->m_cs*p.dinv;
    res -= nu*vx / (x_ab.squaredNorm() + 0.01*m_grad_kern.h2);
  }

  //p.m = p.m - massb * res * m_grad_kern(x_ab);
}
template<typename REAL, typename SIZE, int FT>
inline void CFPressureAccelRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p)
{  
  for (unsigned char i = 0; i < 3; ++i)
  {
    p._extern_accel[i] += global::dynset.gravity[i];// + p.m[i];
    p._accel[i] += p._extern_accel[i] + p.n[i]/this->m_mass;
  }
}

#define INSTANTIATE_PROC_TYPE( proc, type ) \
  template class proc<double, unsigned int, type>
  
#define INSTANTIATE_PROC( proc ) \
  INSTANTIATE_PROC_TYPE( proc, 0 ); \
  INSTANTIATE_PROC_TYPE( proc, 1 ); \
  INSTANTIATE_PROC_TYPE( proc, 2 ); \
  INSTANTIATE_PROC_TYPE( proc, 3 ); \

INSTANTIATE_PROC(CFAccelRST)
INSTANTIATE_PROC(CFSurfaceNormalRST)
INSTANTIATE_PROC(CFSurfaceTensionRST)

INSTANTIATE_PROC(CFPrepareJacobiRST)
INSTANTIATE_PROC(CFJacobiSolveFirstRST)
INSTANTIATE_PROC(CFJacobiSolveSecondRST)
INSTANTIATE_PROC(CFPressureAccelRST)
INSTANTIATE_PROC(CFDensityRST)
INSTANTIATE_PROC(CFDensityUpdateRST)
