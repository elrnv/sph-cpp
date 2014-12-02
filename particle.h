#ifndef PARTICLE_H
#define PARTICLE_H

#include "eigen.h"
#include "math.h"

template<int PT>
struct FluidDataT;

struct __attribute__ ((__packed__)) Particle
{
  Particle(const Vector3T<Real> &p) : pos(p){ }
  ~Particle() { }

  Vector3T<Real> pos;
  Real dinv;
  Real d;
  Real vol; // volume estimate used for shepard filter
  Real pressure;
};

struct __attribute__ ((__packed__)) FluidParticle : public Particle
{
  explicit FluidParticle(
      const Vector3T<Real> &p, const Vector3T<Real> &v,
      Real *a, Real *ea, Real *dinv, Real c)
    : Particle(p)
    , vel(v)
    , _accel(a)
    , _extern_accel(ea)
    , _dinv(dinv)
    , color(c) { }
  ~FluidParticle() { }

  Vector3T<Real> vel;
  Vector3T<Real> n;    // for various computations
  Vector3T<Real> m;    // for various computations
  Real *_accel;        // 3-array to which we write total acceleration
  Real *_extern_accel; // 3-array to which we write acceleration from external forces
  Real *_dinv; // pointer to fluid entry
  Real color;  // used for surface tension computations
  Real c;      // intermediate for various computations
};

struct __attribute__ ((__packed__)) ImplicitFluidParticle : public FluidParticle
{
  explicit ImplicitFluidParticle(
      const Vector3T<Real> &p, const Vector3T<Real> &v,
      Real *a, Real *ea, Real *dinv, Real c)
    : FluidParticle(p,v,a,ea,dinv,c) { }
  ~ImplicitFluidParticle() { }

  Real b;      // intermediate for various computations
};


// This describes a general fluid particle
// an explicit specialization must be given to non-fluid particles
template<int PT>
struct ParticleT : public FluidParticle
{ 
  explicit ParticleT(Vector3T<Real> p, Vector3T<Real> v,
      Real *a, Real *ea, Real *dinv, Real c, FluidDataT<PT> &fldata)
    : FluidParticle(p, v, a, ea, dinv, c), fl(&fldata) { }

  FluidDataT<PT> *fl;

  template<int T>
  void init() {  std::cerr << PT << "  default init " << T << std::endl; }
  template<int T>
  void neigh(Particle &neigh) { (void) neigh; std::cerr << PT << " default neigh p " << T << std::endl;  }
  template<int T>
  void neigh(FluidParticle &neigh) { (void) neigh; std::cerr << PT << " default neigh fp " << T << std::endl; }
  template<int T>
  void finish() { std::cerr << PT << " default finish " << T << std::endl;}
};

// particle vector (templated "typedef")
template<int PT>
using ParticleVecT = std::vector< ParticleT<PT> >;

// boundary particles look different than fluid particles:
template<>
struct ParticleT<STATIC> : public Particle
{ 
  explicit ParticleT(const Vector3T<Real> &p, Real h)
    : Particle(p), radius(h) { }

  Real radius;

  template<int T>
  void init() { }
  template<int T>
  void neigh(Particle &neigh) { (void) neigh; }
  template<int T>
  void neigh(FluidParticle &neigh) { (void) neigh; }
  template<int T>
  void finish() { }
};
#endif // PARTICLE_H
