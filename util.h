#ifndef UTIL_H
#define UTIL_H
#include <QDebug>
#include <vector>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/shared_ptr.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <libconfig.h++>
#include "scene.h"
#include "glmesh.h"
#include "glpointcloud.h"
#include "dynparams.h"
#include "settings.h"

namespace Util
{

SceneNode *
processScene( const aiScene *scene, const aiNode *node, DynParamsPtr dyn_params )
{
  SceneNode *scene_node = new SceneNode(node);

  unsigned int num_meshes = node->mNumMeshes;

  // copy the meshes if any
  unsigned int *meshidx = node->mMeshes;
  for (unsigned int i = 0; i < num_meshes; ++i)
  {
    aiMesh *aimesh = scene->mMeshes[meshidx[i]];
    aiMaterial *aimat = scene->mMaterials[aimesh->mMaterialIndex];
    scene_node->add_child(new GeometryNode(aimesh, aimat, dyn_params));
  }

  // Continue for all child nodes
  unsigned int num_children = node->mNumChildren;
  for (unsigned int i = 0; i < num_children; ++i)
    scene_node->add_child(processScene(scene, node->mChildren[i], dyn_params));

  return scene_node;
}

DynParamsPtr
loadDynamics( const std::string &filename )
{
  // parse the dynamics file
  std::ifstream file(filename, std::ifstream::in);

  if (!file.is_open())
  {
    qWarning() << "Could not open dynamics file:" << filename.c_str();
    return DynParamsPtr(NULL);
  }

  DynParamsPtr params_ptr;

  // determine type of dynamics
  std::string line;
  while (getline(file, line))
  {
    std::istringstream iss(line);
    std::string var_name;
    iss >> var_name;

    if (boost::iequals(var_name, std::string("fluid")))
    {
      params_ptr = FluidParamsPtr(new FluidParams());
      continue;
    }

    if (!params_ptr)
      continue;

    if (boost::iequals(var_name, "vel"))
    {
      iss >> params_ptr->velocity[0];
      iss >> params_ptr->velocity[1];
      iss >> params_ptr->velocity[2];
    }
    else if (boost::iequals(var_name, "angvel"))
    {
      iss >> params_ptr->angular_velocity[0];
      iss >> params_ptr->angular_velocity[1];
      iss >> params_ptr->angular_velocity[2];
    }
  }

  if (params_ptr && params_ptr->type == DynParams::FLUID)
  {
    FluidParams &fparams = static_cast<FluidParams&>(*params_ptr);

    file.clear();
    file.seekg(0); // rewind

    while (getline(file, line))
    {
      std::istringstream iss(line);
      std::string var_name;
      iss >> var_name;

      if (params_ptr->type == DynParams::FLUID)
      {
        if (boost::iequals(var_name, "t"))
        {
          std::string fluid_type;
          iss >> fluid_type;
          if (boost::iequals(fluid_type, "mcg03"))
            fparams.fluid_type = MCG03;
          else if (boost::iequals(fluid_type, "bt07"))
            fparams.fluid_type = BT07;
          else if (boost::iequals(fluid_type, "aiast12"))
            fparams.fluid_type = AIAST12;
        }
        else if (boost::iequals(var_name, "d"))
          iss >> fparams.density;
        else if (boost::iequals(var_name, "v"))
          iss >> fparams.viscosity;
        else if (boost::iequals(var_name, "st"))
          iss >> fparams.surface_tension;
        else if (boost::iequals(var_name, "cs"))
          iss >> fparams.sound_speed;
        else if (boost::iequals(var_name, "c"))
          iss >> fparams.compressibility;
        else if (boost::iequals(var_name, "ki"))
          iss >> fparams.kernel_inflation;
        else if (boost::iequals(var_name, "rvd"))
          iss >> fparams.recoil_velocity_damping;
      }
    }
  }

  file.close();
  return params_ptr;
}

SceneNode *
loadObject( const std::string &filename, DynParamsPtr params_ptr )
{
  Assimp::Importer importer;

  importer.SetPropertyInteger(AI_CONFIG_PP_FD_REMOVE, 1);

  // Read the file into an assimp scene data structure
  const aiScene *scene = 
      importer.ReadFile( filename, 0
        | aiProcess_CalcTangentSpace     
        | aiProcess_GenSmoothNormals // or GenNormals
        | aiProcess_JoinIdenticalVertices
        | aiProcess_Triangulate
        | aiProcess_GenUVCoords
        | aiProcess_SortByPType
        | aiProcess_FindDegenerates
          );

  if (!scene)
  {
    qWarning() << importer.GetErrorString();
    return NULL;
  }

  return processScene(scene, scene->mRootNode, params_ptr);
}

SceneNode *
loadScene( const std::string &filename )
{
  SceneNode *root = NULL;
  libconfig::Config cfg;
  try
  {
    cfg.readFile(filename.c_str());

    // Dynamic Settings
    global::dynset.fps = 60;
    global::dynset.frames = 250;
    global::dynset.substeps = 10;
    global::dynset.cachedir = "";
    cfg.lookupValue("dynamics.fps", global::dynset.fps);
    cfg.lookupValue("dynamics.frames", global::dynset.frames);
    cfg.lookupValue("dynamics.substeps", global::dynset.substeps);
    cfg.lookupValue("dynamics.cachedir", global::dynset.cachedir);

    global::sceneset.padx = Vector2f(0,0);
    global::sceneset.pady = Vector2f(0,0);
    global::sceneset.padz = Vector2f(0,0);

    try
    {
      libconfig::Setting &padset = cfg.lookup("scene.padding");
      global::sceneset.padx = Vector2f(padset["x"][0], padset["x"][1]);
      global::sceneset.pady = Vector2f(padset["y"][0], padset["y"][1]);
      global::sceneset.padz = Vector2f(padset["z"][0], padset["z"][1]);
    }
    catch (libconfig::SettingTypeException &te)
    {
      qWarning() << "Setting " << te.getPath() << "has the wrong type!";
    }
    catch (libconfig::SettingNotFoundException &nfe) { }
    
    global::sceneset.rotx = 0.0f;
    global::sceneset.roty = 0.0f;
    cfg.lookupValue("scene.rotation.x", global::sceneset.rotx);
    cfg.lookupValue("scene.rotation.y", global::sceneset.roty);
    
    std::string objdir = cfg.lookup("scene.objdir");
    libconfig::Setting& objset = cfg.lookup("scene.objects");

    int num = objset.getLength();
    if (!num)
      return root;

    // get scene objects
    for (int i = 0; i < num; ++i)
    {
      try
      {
        std::string dynfile = objset[i]["dynfile"];

        DynParamsPtr params_ptr(loadDynamics( objdir + "/" + dynfile ));

        std::string objfile = objset[i]["objfile"];
        
        if (params_ptr)
        {
          size_t start = objfile.find_last_of("/");
          size_t end = objfile.find_last_of(".");
          if( start == std::string::npos )
            start = -1;
          params_ptr->cacheprefix = objfile.substr(start+1, end - start - 1);
        }

        SceneNode *obj = loadObject( objdir + "/" + objfile, params_ptr );
        if(obj)
        {
          if (!root)
            root = new SceneNode("root");

          root->add_child(obj);
        }
      }
      catch (libconfig::SettingTypeException &te)
      {
        qWarning() << "Setting " << te.getPath() << "has the wrong type!";
      }
      catch (libconfig::SettingNotFoundException &nfe)
      {
        qWarning() << "Setting " << nfe.getPath() << "not found!";
      }
    }
  }
  catch (libconfig::ParseException &pe) 
  {
    qWarning() << "Configuration file:" << filename.c_str() << ", could not be parsed:";
    qWarning() << "libconfig:" << pe.getError();
  }
  catch (libconfig::FileIOException &ioe) 
  {
    qWarning() << "Configuration file:" << filename.c_str() << ", could not be opened!";
  }
  return root;
}

// PRE: we must have the current GL context
void loadGLData(
    const SceneNode *node,
    std::vector< GLPrimitive * > &gl_prims,
    UniformBuffer &ubo,
    ShaderManager &shaderman)
{
  if (node->is_geometry())
  {
    const GeometryNode *geonode = static_cast<const GeometryNode*>(node);
    PrimitivePtr prim = geonode->get_primitive();
    if (prim)
    {
      if (prim->is_mesh())
      {
        gl_prims.push_back(new GLMesh(boost::static_pointer_cast<Mesh>(prim),
                                  geonode->get_material(), ubo, shaderman));
      }
      else if (prim->is_pointcloud())
      {
        gl_prims.push_back(new GLPointCloud(boost::static_pointer_cast<PointCloud>(prim),
                                        geonode->get_material(), ubo, shaderman));
      
      }
    }
  }

  for ( const SceneNode *child : node->get_children())
    loadGLData(child, gl_prims, ubo, shaderman);
}

};

#endif // UTIL_H
