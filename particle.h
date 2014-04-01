#ifndef PARTICLE_H
#define PARTICLE_H

#include "eigen.h"
#include "dynparams.h"

template<typename REAL,typename SIZE>
class FluidRS;

template<typename REAL>
struct ParticleR
{
  ParticleR(Vector3R<REAL> p) : pos(p){ }
  ~ParticleR() { }

  Vector3R<REAL> pos;
  REAL dinv;
  REAL pressure;
};

template<typename REAL, typename SIZE>
struct FluidParticleRS : public ParticleR<REAL>
{
  explicit FluidParticleRS(Vector3R<REAL> p, Vector3R<REAL> v, REAL *a, REAL *ea,
      FluidRS<REAL,SIZE> *f)
    : ParticleR<REAL>(p), vel(v), accel(a), extern_accel(ea), fl(f) { }
  ~FluidParticleRS() { }

  Vector3R<REAL> vel;
  REAL *accel;        // 3-array to which we write total acceleration
  REAL *extern_accel; // 3-array to which we write acceleration from external forces
  REAL temp;          // store temporary values at the particle during computation
  REAL vol;           // volume estimate
  FluidRS<REAL,SIZE> *fl; // fluid to which this particle belongs to
};


#endif // PARTICLE_H
