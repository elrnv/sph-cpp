#ifndef PARTICLE_H
#define PARTICLE_H

#include "eigen.h"

template<typename REAL>
struct ParticleR
{
  ParticleR(Vector3R<REAL> p) : pos(p){ }
  ~ParticleR() { }

  Vector3R<REAL> pos;
  REAL dinv;
  REAL pressure;
};

template<typename REAL>
struct FluidParticleR : public ParticleR<REAL>
{
  FluidParticleR(Vector3R<REAL> p, Vector3R<REAL> v, REAL *a, REAL *ea,
      unsigned short i)
    : ParticleR<REAL>(p), vel(v), accel(a), extern_accel(ea), id(i) { }
  ~FluidParticleR() { }

  Vector3R<REAL> vel;
  REAL *accel;        // 3-array to which we write total acceleration
  REAL *extern_accel; // 3-array to which we write acceleration from external forces
  REAL temp;          // store temporary values at the particle during computation
  REAL vol;           // volume estimate
  unsigned short id;  // index of the fluid this particle belongs to in the grid
};

#endif // PARTICLE_H