#ifndef TYPES_H
#define TYPES_H

#include "eigen.h"

// This file contains type defines used throughout the program. Most importantly,
// floating point types, integral types and some constants

// One radian measured in degrees (conversion constant)
#define RADIAN 0.017453292519943

typedef double Real;
typedef unsigned int Size;
typedef Size Index;

#define INVALID_INDEX Index(-1)

#define ALL_FLUID_PARTICLE_TYPES DEFAULT, MCG03, BT07
#define NUMFLUIDTYPES 3
#define ALL_PARTICLE_TYPES ALL_FLUID_PARTICLE_TYPES, STATIC
#define NUMTYPES 4

// SPH particle types
enum ParticleType
{
  ALL_PARTICLE_TYPES
};

// SPH compute types
enum ComputeType
{
  Density,
  Accel,
  Volume
};

extern AlignedBox3f UnitBox;

#endif // TYPES_H
