#ifndef DYNPARAMS_H
#define DYNPARAMS_H

#include <boost/shared_ptr.hpp>
#include <string>
#include "eigen.h"

// General object specific parameters

struct DynParams
{
  std::string cacheprefix;

  enum Type
  {
    NONE,
    FLUID,
    RIGID,
    STATIC
  } type;

  Vector3f velocity;
  Vector3f angular_velocity;

  DynParams() // Default values:
   : type(NONE)
   , velocity(0,0,0)
   , angular_velocity(0,0,0)
  { }

  ~DynParams() { }
};


// Fluid specific parameters

// When adding fluid types, make sure you add them here, in the macros below and:
// fluid.cpp:
//   add a relevant fluid template instantiation
enum __attribute__ ((__packed__)) FluidType
{
  NOTFLUID = -1,
  DEFAULT = 0,
  MCG03 = 1,
  BT07 = 2,
  AIAST12 = 3,
  NUMTYPES = 4
};


struct FluidParams : public DynParams
{
  FluidType fluid_type;

  float density;
  float viscosity;
  float surface_tension;
  float sound_speed;
  float compressibility;
  float kernel_inflation;
  float recoil_velocity_damping;

  FluidParams() // Default values:
   : fluid_type(DEFAULT)
   , density(1000.0f)
   , viscosity(0.5f)
   , surface_tension(0.0728f)
   , sound_speed(8.0f)
   , compressibility(0.01f)
   , kernel_inflation(3.0f)
   , recoil_velocity_damping(0.4f)
  {
    this->type = FLUID;
  }

  ~FluidParams() { }
};

// DynParams should not need to be owned by anybody
typedef boost::shared_ptr<DynParams> DynParamsPtr;
typedef boost::shared_ptr<FluidParams> FluidParamsPtr;

#endif // DYNPARAMS_H
