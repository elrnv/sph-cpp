#ifndef SETTINGS_H
#define SETTINGS_H
#include <string>
#include "eigen.h"

// data structure containing a global scene configuration

struct DynSettings
{
  unsigned int frames;
  unsigned int fps;
  unsigned int substeps;
  unsigned int init_steps;
  Vector3f gravity;
  std::string savedir;
};

struct SceneSettings
{
  Vector2f padx;
  Vector2f pady;
  Vector2f padz;
  float rotx;
  float roty;
  bool normalize;
};

namespace global
{
  extern DynSettings dynset;
  extern SceneSettings sceneset;
};

#endif // SETTINGS_H
