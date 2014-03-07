#ifndef UTIL_H
#define UTIL_H
#include <QDebug>
#include <vector>
#include <boost/shared_ptr.hpp>
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

SceneNode *
processScene( const aiScene *scene, const aiNode *node )
{
  SceneNode *scene_node = new SceneNode(node);

  unsigned int num_meshes = node->mNumMeshes;

  // copy the meshes if any
  unsigned int *meshidx = node->mMeshes;
  for (unsigned int i = 0; i < num_meshes; ++i)
  {
    aiMesh *aimesh = scene->mMeshes[meshidx[i]];
    aiMaterial *aimat = scene->mMaterials[aimesh->mMaterialIndex];
    GeometryNode *geo_node = new GeometryNode(aimesh, aimat);
    scene_node->add_child(geo_node);
  }

  // Continue for all child nodes
  unsigned int num_children = node->mNumChildren;
  for (unsigned int i = 0; i < num_children; ++i)
    scene_node->add_child(processScene(scene, node->mChildren[i]));

  return scene_node;
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

  return processScene(scene, scene->mRootNode);
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
