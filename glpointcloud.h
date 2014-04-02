#ifndef GLPOINTCLOUD_H
#define GLPOINTCLOUD_H

#include <iostream>
#include <assimp/mesh.h>
#include "pointcloud.h"
#include "glprimitive.h"
#include "dynparams.h"
#include "dynamics.h"

template<typename REAL, typename SIZE>
class FluidRS;

// A point cloud representation for OpenGL applications
template<typename REAL, typename SIZE>
class GLPointCloudRS : public GLPrimitiveS<SIZE>
{
public:
  explicit GLPointCloudRS(
      PointCloudPtrRS<REAL, SIZE> pc,
      MaterialConstPtr mat,
      UniformBuffer &ubo,
      ShaderManager &shaderman);
  ~GLPointCloudRS();

  inline bool is_pointcloud() const { return true; }
  inline bool is_dynamic()    const { return m_isdynamic; }
  inline REAL get_radius()    const { return m_radius; }
  inline REAL get_halo_radius() const { return m_halo_radius; }

  inline SIZE get_num_indices()  const { return get_num_vertices(); }
  inline SIZE get_num_vertices() const { return m_pc->get_num_vertices(); }

  Vector3f get_closest_pt() const;
  void sort_by_depth(const AffineCompact3f &mvmtx);
  void update_glbuf();
  void update_shader(ShaderManager::ShaderType type);

  //inline PointCloudRS<REAL,SIZE> *init_dynamics();
  //FluidRS<REAL,SIZE> *make_dynamic(FluidParamsPtr params);

  void print() const { std::cerr << *m_pc << std::endl; }

  void update_data();

  bool is_halos() const { return m_halos; }
  void toggle_halos() { m_halos = !m_halos; }

  friend UniformGridRS<REAL,SIZE>;

protected:
  PointCloudPtrRS<REAL, SIZE> m_pc; // reference to the native point cloud object

  // intermediate buffer between pc and glbuffer
  Matrix3XR<GLfloat> m_vertices;

  REAL m_radius;      // particle radius inherited from point cloud
  REAL m_halo_radius; // particle influence radius inherited from point cloud
  bool m_insync;
  bool m_halos; // show kernel influence visually
  bool m_isdynamic;
};

// defaults
typedef GLPointCloudRS<double, unsigned int> GLPointCloud;

#endif // GLPOINTCLOUD_H
