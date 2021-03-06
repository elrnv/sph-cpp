#include "types.h"
#include "fluiddata.h"
#include "settings.h"
#include "particle.h"

// Particle Type specializations implementing the interface given above

// MCG03
// Density
template<> template<> void 
ParticleT<MCG03>::init<Density>()
{
  dinv = 0.0;
}

template<> template<> void 
ParticleT<MCG03>::neigh<Density>(FluidParticle &neigh)
{
  dinv += fl->m_kern[ pos - neigh.pos ];
}
template<> template<> void 
ParticleT<MCG03>::neigh<Density>(Particle &neigh) { (void) neigh; }
template<> template<> void 
ParticleT<MCG03>::finish<Density>()
{
  *_dinv = dinv = 1.0 / (fl->m_mass * dinv * fl->m_kern.coef);
  vol = dinv * fl->m_mass;

  pressure = std::max(Real(0), fl->m_cs2 * (1.0/dinv - fl->m_rest_density));

  // TODO: implement this
  //Real var = std::abs(1.0/dinv - fl->m_rest_density);
  //if ( var > *max_var )
  //  *max_var = var;
  //*avg_var += var;
}

// Accel
template<> template<> void 
ParticleT<MCG03>::init<Accel>()
{ 
  c = vol * Real(color) * fl->m_color_kern(Vector3T<Real>(0.0,0.0,0.0));
  n[0] = n[1] = n[2] = 0.0;
  m[0] = vol * fl->m_color_kern(Vector3T<Real>(0.0,0.0,0.0));
}

template<> template<> void 
ParticleT<MCG03>::neigh<Accel>(FluidParticle &neigh)
{
  if (this == &neigh)
    return;
  
  Vector3T<Real> v_ab(vel - neigh.vel);
  Vector3T<Real> x_ab(pos - neigh.pos);

  Vector3T<Real> res(
  // pressure
      neigh.vol * 0.5*(pressure + neigh.pressure)*fl->m_spikygrad_kern(x_ab)
      +
  // viscosity
      fl->m_viscosity * neigh.vol * v_ab * fl->m_visclap_kern(x_ab) 
  );

  for (unsigned char i = 0; i < 3; ++i)
    _accel[i] -= res[i]; // copy intermediate result

  // surface tension of current phase
  // (apply surface tension only to particles from the same phase)
  if (color == neigh.color)
    n = n - (neigh.vol/(neigh.dinv * fl->m_mass))*fl->m_st * fl->m_color_kern(x_ab) * x_ab;

  // surface tension based on color field (intermediate values)
  c += neigh.vol * Real(neigh.color) * fl->m_color_kern(x_ab);
  m[0] += neigh.vol * fl->m_color_kern(x_ab); // denom for smoothed color
}

template<> template<> void 
ParticleT<MCG03>::neigh<Accel>(Particle &neigh) { (void) neigh; }
template<> template<> void
ParticleT<MCG03>::finish<Accel>()
{
  //qDebug() << "c = " << c;
  //c = c / m[0];

  //qDebug() << "c/normal = " << c;

  _extern_accel[0] += global::dynset.gravity[0];
  _extern_accel[1] += global::dynset.gravity[1];
  _extern_accel[2] += global::dynset.gravity[2];
  for (unsigned char i = 0; i < 3; ++i)
    _accel[i] = ((vol / fl->m_mass) * _accel[i]) + _extern_accel[i];

  //qDebug() << "_accel = " << _accel[0] << _accel[1] << _accel[2];
}


// BT07
#define BT07_BOUNDARY_PARTICLES

// Density
template<> template<> void
ParticleT<BT07>::init<Density>()
{ 
  dinv = 0.0;
  vol = 0.0;
}
template<> template<> void
ParticleT<BT07>::neigh<Density>(FluidParticle &neigh)
{
  dinv += fl->m_mass * fl->m_kern[ pos - neigh.pos ];
  //std::cerr << "dist = " << (pos - neigh.pos).transpose() << std::endl;
  //std::cerr << "kern = " << fl->m_kern[ pos - neigh.pos ] << std::endl;
  vol += fl->m_kern[ pos - neigh.pos ];
}

template<> template<> void
ParticleT<BT07>::neigh<Density>(Particle &neigh) 
{
  // neigh.dinv is inverse number density (aka volume) of the boundary paticle
#ifdef BT07_BOUNDARY_PARTICLES
  dinv += fl->m_rest_density * neigh.dinv * fl->m_kern[ pos - neigh.pos ];
  vol += fl->m_kern[ pos - neigh.pos ];
#endif
}
template<> template<> void
ParticleT<BT07>::finish<Density>()
{
  *_dinv = dinv = 1.0 / (dinv * fl->m_kern.coef);
  vol = 1.0 / (vol * fl->m_kern.coef);
  //std::cerr << "dinv = " << dinv << std::endl;

  pressure = std::max(Real(0),
    fl->m_rest_density * fl->m_cs2 * 0.14285714285714 * 
    (pow7(1.0 / (dinv * fl->m_rest_density)) - 1.0));

#ifdef REPORT_DENSITY_VARIATION
  Real var = std::abs(1.0/dinv - fl->m_rest_density);
  if ( var > *fl->m_max_var )
    *fl->m_max_var = var;
  *fl->m_avg_var += var;
#endif
}

// Accel
template<> template<> void
ParticleT<BT07>::init<Accel>()
{
  n[0] = n[1] = n[2] = 0.0;
  m[0] = m[1] = m[2] = 0.0;
}

template<> template<> void
ParticleT<BT07>::neigh<Accel>(FluidParticle &neigh)
{
  if (this == &neigh)
    return;

  // pressure

  Real res = (pressure*vol*vol + neigh.pressure*neigh.vol*neigh.vol)/fl->m_mass;

  // viscosity

  Vector3T<Real> x_ab = pos - neigh.pos;
  Real vx = x_ab.dot(vel - neigh.vel);
  if (vx < 0)
  {
    Real nu = 2*fl->m_viscosity*fl->m_grad_kern.h*fl->m_cs *
            neigh.vol * dinv / (dinv + neigh.dinv);
    res -= nu*vx / (x_ab.squaredNorm() + 0.01*fl->m_grad_kern.h2);
  }

  n = n - res * fl->m_grad_kern(x_ab);
  
  // surface tension 
  // (apply surface tension only to particles from the same phase)
  if (color == neigh.color)
    n = n - (neigh.vol/(neigh.dinv*fl->m_mass))*fl->m_st * fl->m_kern(x_ab) * x_ab;
}
template<> template<> void
ParticleT<BT07>::neigh<Accel>(Particle &neigh)
{
  // BT07 fluids dont interact with static boundary particles, instead they
  // mirror fluid particles on the boundary, and compute repulsion forces
  // this way (in the finish_particle method)

#ifdef BT07_BOUNDARY_PARTICLES // the following handles boundaries as in AIAST12
  Vector3T<Real> x_ab = pos - neigh.pos;
  float massb = fl->m_rest_density * neigh.dinv;

  // pressure force contribution from boundary particles
  Real res( pressure * dinv * dinv);

  Real vx = x_ab.dot(vel);
  // viscosity foces contribution from boundary particles
  
  if (vx < 0)
  {
    Real nu = 0.5*fl->m_friction*fl->m_grad_kern.h*fl->m_cs*dinv;
    res -= nu*vx / (x_ab.squaredNorm() + 0.01*fl->m_grad_kern.h2);
  }

  m = m - massb * res * fl->m_grad_kern(x_ab);
#endif
}

template<> template<> void
ParticleT<BT07>::finish<Accel>()
{
  // if the particle is close to the boundary and add a repulsive force
  // (assume boundary particle is infinitely more massive than the fluid
  // particle, such that m_b / (m_a + m_b) == 1 where m_a and m_b are the
  // masses of the fluid and boundary particles respectively
#ifndef BT07_BOUNDARY_PARTICLES
  Vector3T<Real> hvec(fl->m_bound_kern.h,fl->m_bound_kern.h,fl->m_bound_kern.h);
  Vector3T<Real> res(0.0, 0.0, 0.0);
  Vector3T<Real> dmin =  // distances to min boundaries
    0.5*hvec + pos - fl->m_bmin.template cast<Real>();
  Vector3T<Real> dmax =  // distances to max boundaries
    0.5*hvec + fl->m_bmax.template cast<Real>() - pos;
  if ( dmin[0] < fl->m_bound_kern.h )
    res[0] += fl->m_bound_kern(Vector3d(dmin[0],0,0));// / dmin[0];
  else if ( dmax[0] < fl->m_bound_kern.h)
    res[0] -= fl->m_bound_kern(Vector3d(dmax[0],0,0));// / dmax[0];

  if ( dmin[1] < fl->m_bound_kern.h )
    res[1] += fl->m_bound_kern(Vector3d(dmin[1],0,0));// / dmin[1];
  else if ( dmax[1] < fl->m_bound_kern.h)
    res[1] -= fl->m_bound_kern(Vector3d(dmax[1],0,0));// / dmax[1];

  if ( dmin[2] < fl->m_bound_kern.h )
    res[2] += fl->m_bound_kern(Vector3d(dmin[2],0,0));// / dmin[2];
  else if ( dmax[2] < fl->m_bound_kern.h)
    res[2] -= fl->m_bound_kern(Vector3d(dmax[2],0,0));// / dmax[2];

  // acceleration due to boundary collisions
  Vector3T<Real> bdry_accel( fl->m_compressibility*fl->m_cs2*res );
#endif
  for (unsigned char i = 0; i < 3; ++i)
  {
    _extern_accel[i] += global::dynset.gravity[i] +
#ifdef BT07_BOUNDARY_PARTICLES
      m[i];
#else
      bdry_accel[i];
#endif
    _accel[i] += _extern_accel[i] + n[i];
  }
}

// Boundary

// From [Solenthaler and Pajarola 2008], an alternative density
// (number_density * mass) is used
template<> void
ParticleT<STATIC>::init<Volume>() 
{ 
  dinv = 0.0; 
}
template<> void 
ParticleT<STATIC>::neigh<Volume>(FluidParticle &neigh) { (void) neigh; }
template<> void
ParticleT<STATIC>::neigh<Volume>(Particle &neigh)
{
  dinv += CubicSplineKernel::compute( pos - neigh.pos, radius );
}
template<> void
ParticleT<STATIC>::finish<Volume>()
{
  dinv = 1.0 / (CubicSplineKernel::compute_coef(radius) * dinv);
}
