#ifndef GLPOINTCLOUD_H
#define GLPOINTCLOUD_H

#include <iostream>
#include <assimp/mesh.h>
#include "pointcloud.h"
#include "glprimitive.h"
#include "dynparams.h"

class UniformGrid;
class Fluid;

// A point cloud representation for OpenGL applications

class GLPointCloud : public GLPrimitive
{
public:
  explicit GLPointCloud(
      PointCloud &pc,
      bool dynamic,
      MaterialManager &matman,
      UniformBuffer &ubo,
      ShaderManager &shaderman,
      float halo_radius = 0.0f);
  ~GLPointCloud();

  inline bool is_pointcloud() const { return true; }
  inline bool is_dynamic()    const { return m_isdynamic; }
  inline Real get_radius()    const { return m_radius; }
  inline Real get_halo_radius() const { return m_halo_radius; }

  inline Size get_num_indices()  const { return get_num_vertices(); }
  inline Size get_num_vertices() const { return m_pc.get_num_vertices(); }

  //void sort_by_z(Matrix3XT<GLfloat> &vispos);
  void sort_by_z(const Matrix3XT<GLfloat> &pos, Matrix3XT<GLfloat> &outpos);
  void update_glbuf_withsort(const AffineCompact3f &trans,
                             const AffineCompact3f &nmltrans);
  void update_glbuf_nosort();
  void update_shader(ShaderManager::ShaderType t, ShaderManager &sm);

  void clear_cache();

  void print() const { std::cerr << m_pc << std::endl; }

  friend UniformGrid;

private:  // helper functions
  void update_glbuf(const Matrix3XT<GLfloat> &vispos);

private:
  PointCloud &m_pc; // reference to the native point cloud object

  Real m_radius;      // particle radius inherited from point cloud
  Real m_halo_radius; // particle influence radius inherited from point cloud
  bool m_insync;
  bool m_isdynamic;
};

typedef boost::shared_ptr< GLPointCloud > GLPointCloudPtr;

#endif // GLPOINTCLOUD_H
