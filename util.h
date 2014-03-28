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
#include "scene.h"
#include "glmesh.h"
#include "glpointcloud.h"
#include "dynparams.h"

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
    GeometryNode *geo_node = new GeometryNode(aimesh, aimat, dyn_params);
    scene_node->add_child(geo_node);
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
  std::string line;
  std::ifstream file(filename, std::ifstream::in);

  if (!file.is_open())
  {
    qWarning() << "Could not open dynamics file:" << QString(filename.c_str());
    return DynParamsPtr(NULL);
  }

  DynParamsPtr params_ptr;

  // determine type of dynamics
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

    if (!params_ptr.get())
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

  if (params_ptr.get() && params_ptr->type == DynParams::FLUID)
  {
    FluidParamsPtr fparams_ptr = boost::static_pointer_cast<FluidParams>(params_ptr);

    file.clear();
    file.seekg(0); // rewind

    while (getline(file, line))
    {
      std::istringstream iss(line);
      std::string var_name;
      iss >> var_name;

      if (params_ptr->type == DynParams::FLUID)
      {
        if (boost::iequals(var_name, "t") == 0)
        {
          std::string fluid_type;
          iss >> fluid_type;
          if (boost::iequals(fluid_type, "mcg03"))
            fparams_ptr->fluid_type = MCG03;
          else if (boost::iequals(fluid_type, "bt07"))
            fparams_ptr->fluid_type = BT07;
          else if (boost::iequals(fluid_type, "aiast12"))
            fparams_ptr->fluid_type = AIAST12;
        }
        else if (boost::iequals(var_name, "d"))
          iss >> fparams_ptr->density;
        else if (boost::iequals(var_name, "v"))
          iss >> fparams_ptr->viscosity;
        else if (boost::iequals(var_name, "st"))
          iss >> fparams_ptr->surface_tension;
        else if (boost::iequals(var_name, "cs"))
          iss >> fparams_ptr->sound_speed;
        else if (boost::iequals(var_name, "c"))
          iss >> fparams_ptr->compressibility;
        else if (boost::iequals(var_name, "ki"))
          iss >> fparams_ptr->kernel_inflation;
        else if (boost::iequals(var_name, "rvd"))
          iss >> fparams_ptr->recoil_velocity_damping;
      }
    }
  }

  file.close();
  return params_ptr;
}

SceneNode *
loadScene( const std::string &filename )
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

  // check for dynamics
  std::string line;
  std::ifstream file;
  file.open(filename);

  if (!file.is_open())
  {
    qWarning() << "Could not open input file:" << QString(filename.c_str());
    return NULL;
  }

  DynParamsPtr params_ptr;

  while (getline(file, line))
  {
    if (line.compare(0, 6, "dynlib") == 0)
    {
      // search for dynamics file in the same directory
      std::string dyn_filename = boost::trim_copy(line.substr(7));
      size_t pos = filename.find_last_of("/");
      if (pos == std::string::npos)
        params_ptr = loadDynamics( dyn_filename );
      else
        params_ptr = loadDynamics( filename.substr(0,pos+1) + dyn_filename );
      break;
    }
  }

  file.close();

  return processScene(scene, scene->mRootNode, params_ptr);
}

// PRE: we must have the current GL context
void loadGLData(
    const SceneNode *node,
    std::vector< boost::shared_ptr<GLPrimitive> > &gl_prims,
    UniformBuffer &ubo,
    ShaderManager &shaderman)
{
  if (node->is_geometry())
  {
    const GeometryNode *geonode = static_cast<const GeometryNode*>(node);
    Primitive *prim = geonode->get_primitive();
    if (prim)
    {
      if (prim->is_mesh())
      {
        boost::shared_ptr<GLPrimitive> mesh(
            new GLMesh(static_cast<Mesh*>(prim),
                       static_cast<const Material*>(geonode->get_material()),
                       ubo, shaderman));
      
        gl_prims.push_back(mesh);
      }
      else if (prim->is_pointcloud())
      {
        boost::shared_ptr<GLPrimitive> pc(
            new GLPointCloud(static_cast<PointCloud *>(prim),
                             static_cast<const Material*>(geonode->get_material()),
                             ubo, shaderman));
      
        gl_prims.push_back(pc);
      }
    }
  }

  for ( const SceneNode *child : node->get_children())
    loadGLData(child, gl_prims, ubo, shaderman);
}

};

#endif // UTIL_H
