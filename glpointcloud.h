#ifndef GLPOINTCLOUD_H
#define GLPOINTCLOUD_H

#include <iostream>
#include <assimp/mesh.h>
#include "pointcloud.h"
#include "glprimitive.h"

template<typename REAL, typename SIZE>
class FluidRS;

// A point cloud representation for OpenGL applications
template<typename REAL, typename SIZE>
class GLPointCloudRS : public GLPrimitiveS<SIZE>
{
public:
  explicit GLPointCloudRS(
      PointCloudRS<REAL, SIZE> *pc,
      const Material *mat,
      UniformBuffer &ubo,
      ShaderManager &shaderman);
  ~GLPointCloudRS();

  PointCloudRS<REAL, SIZE> *get_pointcloud() { return m_pc; }

  inline bool is_pointcloud() const { return true; }

  inline SIZE get_num_indices()  const { return get_num_vertices(); }
  inline SIZE get_num_vertices() const { return m_pc->get_num_vertices(); }

  void sort_by_depth(const AffineCompact3f &mvmtx);
  void update_glbuf();
  void update_shader(ShaderManager::ShaderType type);

  FluidRS<REAL,SIZE> *make_dynamic(REAL d, REAL v, REAL st);

  void print() const { std::cerr << *m_pc << std::endl; }

  void update_data();

  bool is_halos() const { return m_halos; }
  void toggle_halos() { m_halos = !m_halos; }

protected:
  PointCloudRS<REAL, SIZE> *m_pc; // reference to the native mesh object

  // intermediate buffer between pc and glbuffer
  Matrix3XR<GLfloat> m_vertices;

  bool m_insync;

  bool m_halos; // show kernel influence visually
};

// defaults
typedef GLPointCloudRS<double, unsigned int> GLPointCloud;

#endif // GLPOINTCLOUD_H
