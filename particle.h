#ifndef PARTICLE_H
#define PARTICLE_H

#include "types.h"
#include "eigen.h"
#include "dynparams.h"

struct __attribute__ ((__packed__)) Particle
{
  Particle(Vector3R<Real> p) : pos(p){ }
  ~Particle() { }

  Vector3R<Real> pos;
  Real dinv;
  Real d;
  Real vol; // volume estimate used for shepard filter
  Real pressure;
};

struct __attribute__ ((__packed__)) FluidParticle : public Particle
{
  explicit FluidParticle(Vector3R<Real> p, Vector3R<Real> v,
      Real *a, Real *ea, Real *dinv, unsigned int i, Real c)
    : Particle(p)
    , vel(v)
    , _accel(a)
    , _extern_accel(ea)
    , _dinv(dinv)
    , color(c)
    , id(i) { }
  ~FluidParticle() { }

  Vector3R<Real> vel;
  Vector3R<Real> n;    // for various computations
  Vector3R<Real> m;    // for various computations
  Real *_accel;        // 3-array to which we write total acceleration
  Real *_extern_accel; // 3-array to which we write acceleration from external forces
  Real *_dinv; // pointer to fluid entry
  Real c;      // intermediate for various computations
  Real b;      // intermediate for various computations
  Real color;  // used for surface tension computations
  unsigned int id;  // index of the fluid this particle belongs to in the grid
};

template<int FT>
struct FluidParticleT : public FluidParticle
{ 
  explicit FluidParticleT(Vector3R<Real> p, Vector3R<Real> v,
      Real *a, Real *ea, Real *dinv, unsigned int i, Real c)
    : FluidParticle(p, v, a, ea, dinv, i, c) { }
};

// would like to infer the Type of fluid for a particular particle type:
template <typename>
struct extract_fluid_type;

template <int FT>
struct extract_fluid_type< FluidParticleT<FT> >
{
  static const int type = FT;
};



#endif // PARTICLE_H
