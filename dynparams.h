#ifndef DYNPARAMS_H
#define DYNPARAMS_H

#include <boost/shared_ptr.hpp>
#include <boost/functional/hash.hpp>
#include <string>
#include "eigen.h"

// General object specific parameters

struct DynParams
{
  std::string saveprefix;

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

  friend std::size_t hash_value( const DynParams &dp )
  {
    std::size_t seed = 0;
    boost::hash_combine(seed, int(dp.type));
    boost::hash_combine(seed, dp.velocity[0]);
    boost::hash_combine(seed, dp.velocity[1]);
    boost::hash_combine(seed, dp.velocity[2]);
    boost::hash_combine(seed, dp.angular_velocity[0]);
    boost::hash_combine(seed, dp.angular_velocity[1]);
    boost::hash_combine(seed, dp.angular_velocity[2]);
    return seed;
  }
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
  ICS13 = 3,
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
  float friction;
  float kernel_inflation;
  float recoil_velocity_damping;

  FluidParams() // Default values:
   : fluid_type(DEFAULT)
   , density(1000.0f)
   , viscosity(0.5f)
   , surface_tension(0.0728f)
   , sound_speed(8.0f)
   , compressibility(0.01f)
   , friction(0.5f)
   , kernel_inflation(3.0f)
   , recoil_velocity_damping(0.4f)
  {
    this->type = FLUID;
  }

  ~FluidParams() { }

  friend std::size_t hash_value( const FluidParams &fp )
  {
    std::size_t seed = 0;
    boost::hash_combine(seed, int(fp.fluid_type));
    boost::hash_combine(seed, fp.density);
    boost::hash_combine(seed, fp.viscosity);
    boost::hash_combine(seed, fp.surface_tension);
    boost::hash_combine(seed, fp.sound_speed);
    boost::hash_combine(seed, fp.compressibility);
    boost::hash_combine(seed, fp.kernel_inflation);
    boost::hash_combine(seed, fp.recoil_velocity_damping);
    boost::hash_combine(seed, hash_value(static_cast<const DynParams &>(fp)));
    return seed;
  }
};

// DynParams should not need to be owned by anybody
typedef boost::shared_ptr<DynParams>   DynParamsPtr;
typedef boost::shared_ptr<FluidParams> FluidParamsPtr;

#endif // DYNPARAMS_H
