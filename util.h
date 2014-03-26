#ifndef UTIL_H
#define UTIL_H
#include <QDebug>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <string>
#include <iostream>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include "scene.h"
#include "glmesh.h"
#include "glpointcloud.h"

namespace Util
{

struct DynParams
{
  enum Type
  {
    NONE,
    FLUID,
    RIGID,
    STATIC
  } type;

  float density;
  float viscosity;
  float surface_tension;
  float sound_speed;
  Vector3f velocity;
  Vector3f angular_velocity;

  DynParams() // Default values:
   : type(NONE)
   , density(1000.0f)
   , viscosity(0.5f)
   , surface_tension(0.0728f)
   , sound_speed(8.0f)
   , velocity(0,0,0)
   , angular_velocity(0,0,0)
  { }

  ~DynParams() { }
};

SceneNode *
processScene( const aiScene *scene, const aiNode *node, const DynParams &dyn_params )
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

DynParams
loadDynamics( const std::string &filename )
{
  // parse the dynamics file
  string line;
  ifstream file;
  file.open(filename);
  if (!file.is_open())
  {
    qWarning() << "Could not open dynamics file:" << filename;
    return DynParams();
  }

  DynParams params;

  while (getline(file, line))
  {
    if (boost::trim(line).compare("fluid") == 0)
    {
      params.type = DynParams::FLUID;
      continue;
    }

    istringstream iss(line);
    string var_name;
    iss >> var_name;

    if (var_name.compare("d") == 0)
      iss >> params.density;
    else if (var_name.compare("v") == 0)
      iss >> params.viscosity;
    else if (var_name.compare("s") == 0)
      iss >> params.surface_tension;
    else if (var_name.compare("c") == 0)
      iss >> params.sound_speed;
    else if (var_name.compare("vel") == 0)
    {
      iss >> params.velocity[0];
      iss >> params.velocity[1];
      iss >> params.velocity[2];
    }
    else if (var_name.compare("angvel") == 0)
    {
      iss >> params.angular_velocity[0];
      iss >> params.angular_velocity[1];
      iss >> params.angular_velocity[2];
    }
  }

  file.close();
  return params;
}


SceneNode *
loadScene( const std::string &filename )
{
  Assimp::Importer importer;

  importer.SetPropertyInteger(AI_CONFIG_PP_FD_REMOVE, 1);

  // Read the file into an assimp scene data structure
  const aiScene *scene = 
    importer.ReadFile( "data/" + filename, 0
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
  string line;
  ifstream file;
  file.open(filename);
  if (!file.is_open())
  {
    qWarning() << "Could not open input file:" << filename;
    return NULL;
  }

  DynParams params;

  while (getline(file, line))
  {
    if (line.compare(0, 6, "dynlib") == 0)
    {
      params = loadDynamics( boost::trim(line.substr(7)) );
      break;
    }
  }

  file.close();

  return processScene(scene, scene->mRootNode, params);
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
