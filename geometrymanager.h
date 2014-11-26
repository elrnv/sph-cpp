#ifndef GEOMETRYMANAGER_H
#define GEOMETRYMANAGER_H

#include "pointcloud.h"
#include "mesh.h"
#include "types.h"

// GeometryManager keeps track of geometry not being simulated

class GeometryManager
{
public:
  GeometryManager() { }
  ~GeometryManager() { } 

  /// Manager interface
  void add_geometry_object(const aiMesh *mesh, Index mat_idx)
  {
    if ( mesh->mPrimitiveTypes & aiPrimitiveType_POINT )
    {
      add_pointcloud(mesh, mat_idx);
    }
    else if ( mesh->mPrimitiveTypes & aiPrimitiveType_TRIANGLE )
    {
      add_mesh(mesh, mat_idx);
    }
    // else ignore. TODO: extend this to include other dynamic objects
  }

  void add_pointcloud(const aiMesh *mesh, Index mat_idx)
  {
    m_pointclouds.push_back(PointCloud(mesh, mat_idx));
  }

  void add_mesh(const aiMesh *mesh)
  {
    m_meshes.push_back(Mesh(mesh, mat_idx));
  }
  
  typedef std::vector< PointCloud > PointCloudVec;
  typedef std::vector< Mesh >       MeshVec;

  PointCloudVec &get_pointclouds() { return m_pointclouds; }
  MeshVec       &get_meshes()      { return m_meshes; }

private:
  PointCloudVec m_pointclouds;
  MeshVec       m_meshes;
};

#endif // GEOMETRYMANAGER_H
