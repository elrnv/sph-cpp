#ifndef PARTICLE_H
#define PARTICLE_H

#include "eigen.h"
#include "dynparams.h"

template<typename REAL>
struct ParticleR
{
  ParticleR(Vector3R<REAL> p) : pos(p){ }
  ~ParticleR() { }

  Vector3R<REAL> pos;
  REAL dinv;
  REAL vol;           // volume estimate
  REAL pressure;
};

template<typename REAL>
struct FluidParticleR : public ParticleR<REAL>
{
  explicit FluidParticleR(Vector3R<REAL> p, Vector3R<REAL> v,
      REAL *a, REAL *ea, REAL *dinv, unsigned short i)
    : ParticleR<REAL>(p), vel(v), accel(a), extern_accel(ea), _dinv(dinv), id(i) { }
  ~FluidParticleR() { }

  Vector3R<REAL> vel;
  REAL *accel;        // 3-array to which we write total acceleration
  REAL *extern_accel; // 3-array to which we write acceleration from external forces
  REAL *_dinv; // pointer to fluid entry
  REAL c;         // for surface tension computations
  Vector3R<REAL> n;         // for surface tension computations
  unsigned short id;  // index of the fluid this particle belongs to in the grid
};

template<typename REAL, int FT>
struct FluidParticleRT : public FluidParticleR<REAL>
{ 
  explicit FluidParticleRT(Vector3R<REAL> p, Vector3R<REAL> v,
      REAL *a, REAL *ea, REAL *dinv, unsigned short i)
    : FluidParticleR<REAL>(p, v, a, ea, dinv, i) { }
};

// would like to infer the Type of fluid for a particular particle type:
template <typename>
struct extract_fluid_type;

template <typename REAL, int FT>
struct extract_fluid_type< FluidParticleRT<REAL, FT> >
{
  static const int type = FT;
};



#endif // PARTICLE_H
