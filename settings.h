#ifndef SETTINGS_H
#define SETTINGS_H
#include <string>
#include <boost/functional/hash.hpp>
#include "eigen.h"

// data structure containing a global scene configuration

struct DynSettings
{
  std::string hashfile;
  unsigned int frames;
  unsigned int fps;
  unsigned int substeps;
  unsigned int init_steps;
  Vector3f gravity;
  std::string savedir;

  friend std::size_t hash_value( DynSettings const& set )
  {
    std::size_t seed = 0;
    boost::hash_combine(seed, set.frames);
    boost::hash_combine(seed, set.fps);
    boost::hash_combine(seed, set.substeps);
    boost::hash_combine(seed, set.init_steps);
    boost::hash_combine(seed, set.gravity[0]);
    boost::hash_combine(seed, set.gravity[1]);
    boost::hash_combine(seed, set.gravity[2]);
    boost::hash_combine(seed, set.savedir);
    return seed;
  }
};

struct SceneSettings
{
  Vector2f padx;
  Vector2f pady;
  Vector2f padz;
  float rotx;
  float roty;
  bool normalize;

  friend std::size_t hash_value( SceneSettings const& set )
  {
    std::size_t seed = 0;
    boost::hash_combine(seed, set.padx[0]);
    boost::hash_combine(seed, set.padx[1]);
    boost::hash_combine(seed, set.pady[0]);
    boost::hash_combine(seed, set.pady[1]);
    boost::hash_combine(seed, set.padz[0]);
    boost::hash_combine(seed, set.padz[1]);
    boost::hash_combine(seed, set.rotx);
    boost::hash_combine(seed, set.roty);
    boost::hash_combine(seed, set.normalize);
    return seed;
  }
};

namespace global
{
  extern DynSettings dynset;
  extern SceneSettings sceneset;
};

#endif // SETTINGS_H
