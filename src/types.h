#ifndef TYPES_H
#define TYPES_H

#include <boost/preprocessor/seq/size.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/cat.hpp>
#include <boost/preprocessor.hpp>
#include "eigen.h"

// This file contains type defines used throughout the program. Most importantly,
// floating point types, integral types and some constants

// One radian measured in degrees (conversion constant)
#define RADIAN 0.017453292519943

typedef double Real;
typedef unsigned int Size;
typedef Size Index;
#define INVALID_INDEX Index(-1)

////////////////////////////////////////////////////////////////////////////////
// When adding new SPH particle types, simply add to these two defines.
// The enum and strings will be generated for you.
#define SPH_FLUID_TYPES (MCG03)(BT07)
#define SPH_RIGID_TYPES (STATIC)
////////////////////////////////////////////////////////////////////////////////

#define SPH_TYPES SPH_FLUID_TYPES SPH_RIGID_TYPES
#define SPH_TYPE_ENUM(seq) BOOST_PP_SEQ_ENUM(SPH_TYPES)

////////////////////////////////////////////////////////////////////////////////
// The following structures are automatically generated from the given SPH types
// above for use throughout the code
enum SPHParticleType { BOOST_PP_SEQ_ENUM(SPH_TYPES) };
extern const char * SPHParticleTypeString[];
#define NUMFLUIDSPHTYPES BOOST_PP_SEQ_SIZE(SPH_FLUID_TYPES)
#define NUMSPHTYPES BOOST_PP_SEQ_SIZE(SPH_TYPES)
#define ALL_FLUID_SPH_TYPES BOOST_PP_SEQ_ENUM(SPH_FLUID_TYPES)
#define ALL_SPH_TYPES BOOST_PP_SEQ_ENUM(SPH_TYPES)
////////////////////////////////////////////////////////////////////////////////


// SPH compute types
enum ComputeType
{
  Density,
  Accel,
  Volume
};

// globally defined unit box
extern AlignedBox3f UnitBox;


#endif // TYPES_H
