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

  // normalization routines used to fit the dynamic objects within a unit cube
  void normalize_models()
  {
    normalize_models(Vector3f(0.0f, 0.0f, 0.0f),Vector3f(0.0f, 0.0f, 0.0f));
  }
  void normalize_models(
      const Vector2f &ext_x,
      const Vector2f &ext_y,
      const Vector2f &ext_z)
  {
    normalize_models(Vector3f(ext_x[0], ext_y[0], ext_z[0]),
                     Vector3f(ext_x[1], ext_y[1], ext_z[1]));
  }

  void normalize_models(const Vector3f &ext_blf, const Vector3f &ext_trc)
  {
    compute_bbox();
    Vector3f blf( m_bbox.corner(Eigen::AlignedBox3f::BottomLeftFloor) );
    Vector3f trc( m_bbox.corner(Eigen::AlignedBox3f::TopRightCeil) );
    m_bbox.extend( blf - ext_blf );
    m_bbox.extend( trc + ext_trc );

    Vector3f sizevec = m_bbox.sizes();
    transform_models(Affine3f::Identity() * Scaling(2.0f/sizevec.maxCoeff()));
    Vector3f box_center = m_bbox.center();
    transform_models(Affine3f::Identity() * Translation3f(-box_center));
  }

  void transform_models(const Affine3f &trans)
  {
    for ( auto &mesh : m_meshes )
      mesh.transform_in_place(trans);
    for ( auto &pc : m_pointclouds )
      pc.transform_in_place(trans);
  }

  AlignedBox3f &compute_bbox()
  {
    m_bbox.setEmpty();
    for ( auto &mesh : m_meshes )
      m_bbox.extend(mesh.compute_bbox());
    for ( auto &pc : m_pointclouds )
      m_bbox.extend(pc.compute_bbox());
    return m_bbox;
  }

  inline void set_bbox(const AlignedBox3f &bbox) { m_bbox = bbox; }

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

  void add_mesh(const aiMesh *mesh, Index mat_idx)
  {
    m_meshes.push_back(Mesh(mesh, mat_idx));
  }

  void clear()
  {
    m_meshes.clear();
    m_pointclouds.clear();
  }

  inline Size get_num_meshes() { return m_meshes.size(); }
  inline Size get_num_pointclouds() { return m_pointclouds.size(); }
  
  PointCloudVec &get_pointclouds() { return m_pointclouds; }
  MeshVec       &get_meshes()      { return m_meshes; }

private:
  PointCloudVec m_pointclouds;
  MeshVec       m_meshes;
  AlignedBox3f  m_bbox; // bounding box of all stored models
};

#endif // GEOMETRYMANAGER_H
