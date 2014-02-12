#ifndef UTIL_H
#define UTIL_H
#include <QDebug>
#include <vector>
#include <string>
#include <iostream>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include "scene.h"
#include "mesh.h"
#include "util.h"

namespace Util
{

SceneNode *
processScene( const aiScene *scene, const aiNode *node )
{
  unsigned int num_meshes = node->mNumMeshes;
  SceneNode *scene_node = new SceneNode(node);

  // copy the meshes if any
  unsigned int *meshidx = node->mMeshes;
  for (unsigned int i = 0; i < num_meshes; ++i)
  {
    GeometryNode *geo_node = new GeometryNode(scene->mMeshes[meshidx[i]]);
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

  // Read the file into an assimp scene data structure
  const aiScene *scene = 
    importer.ReadFile( filename,
        aiProcess_CalcTangentSpace      |
        aiProcess_Triangulate           |
        aiProcess_JoinIdenticalVertices |
        aiProcess_SortByPType);

  if (!scene)
  {
    qWarning() << importer.GetErrorString();
    return NULL;
  }

  if (!scene->HasMeshes())
    return NULL;

  return processScene(scene, scene->mRootNode);
}



// Hidden functions used to extract vertex position and index data from a mesh
// for OpenGL rendering
template <typename FLOAT, typename INDEX>
static void extractMeshData(const Mesh *mesh, FLOAT *&vertices, FLOAT *&normals, INDEX *&indices)
{
  const VertexVec &verts = mesh->get_verts();

  VertexVec::const_iterator v_it = verts.begin();
  for ( ; v_it != verts.end(); ++v_it)
    for (short j = 0; j < 3; ++j)
    {
      *(vertices++) = FLOAT(v_it->pos[j]);
      *(normals++) = FLOAT(v_it->nml[j]);
    }


  const FaceVec &faces = mesh->get_faces();

  FaceVec::const_iterator f_it = faces.begin();
  for ( ; f_it != faces.end(); ++f_it)
    for (short j = 0; j < 3; ++j)
      *(indices++) = INDEX((*f_it)[j]);

//  FLOAT colors[] = {
//    1.0f, 1.0f, 1.0f,
//    1.0f, 1.0f, 0.0f,
//    1.0f, 0.0f, 1.0f,
//    1.0f, 0.0f, 0.0f,
//    0.0f, 1.0f, 1.0f,
//    0.0f, 1.0f, 0.0f,
//    0.0f, 0.0f, 1.0f,
//    0.0f, 0.0f, 0.0f
//  };
}

template <typename FLOAT, typename INDEX>
static void extractSceneData(SceneNode *node, FLOAT *&vertices, FLOAT *normals, INDEX *&indices)
{
  if (node->is_geometry())
  {
    const Primitive *prim = ((GeometryNode *)node)->get_primitive();
    if (prim->is_mesh())
      extractMeshData((Mesh *) prim, vertices, normals, indices);
  }

  for (NodeList::const_iterator it = node->begin(); it != node->end(); ++it)
    extractSceneData(*it, vertices, normals, indices);
}


// Data extraction interface
template <typename FLOAT, typename INDEX>
void extractGLSceneData(SceneNode *node, FLOAT *vertices, FLOAT *normals, INDEX *indices)
{
  extractSceneData<FLOAT, INDEX>(node, vertices, normals, indices);
}

// Count the number of vertices and indices in the scene
// PRE: this function assumes that num_vertices and num_indices are intialized
template <typename SIZE>
void extractGLSceneDataSize(
    SceneNode *node,
    SIZE &num_vertices, 
    SIZE &num_indices)
{
  if (node->is_geometry())
  {
    const Primitive *prim = ((GeometryNode *)node)->get_primitive();
    if (prim->is_mesh())
      num_vertices += ((Mesh *) prim)->get_verts().size() * 3;
      num_indices  += ((Mesh *) prim)->get_faces().size() * 3;
  }

  for (NodeList::const_iterator it = node->begin(); it != node->end(); ++it)
    extractGLSceneDataSize<SIZE>(*it, num_vertices, num_indices);
}

};

#endif // UTIL_H
