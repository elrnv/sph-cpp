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
  std::string cachedir;
};

struct SceneSettings
{
  Vector2f padx;
  Vector2f pady;
  Vector2f padz;
  float rotx;
  float roty;
};

namespace global
{
  extern DynSettings dynset;
  extern SceneSettings sceneset;
};

#endif // SETTINGS_H
