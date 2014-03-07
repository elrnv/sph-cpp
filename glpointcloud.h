#ifndef POINTCLOUD_H
#define POINTCLOUD_H

#include <iostream>
#include <assimp/mesh.h>
#include "primitive.h"

// Partially resolve matrix template for convenience
template<typename REAL>
using Matrix3XR = Matrix<REAL, 3, Dynamic>;

template<typename REAL, typename SIZE>
class PointCloudRS;

template<typename REAL, typename SIZE>
std::ostream& operator<<(std::ostream& out, const PointCloudRS<REAL,SIZE>& pc);

// A cloud of points
template<typename REAL, typename SIZE>
class PointCloudRS : public Primitive 
{
public:
  explicit PointCloudRS(const aiMesh *pc);
  ~PointCloudRS();

  Matrix3XR<REAL> &get_pos() { return m_pos; }
  inline SIZE get_num_vertices() const { return m_pos.cols(); }
  inline bool is_pointcloud() const { return true; }

  friend std::ostream& operator<< <>(std::ostream& out, const PointCloudRS<REAL,SIZE>& pc);

protected:
  void compute_bbox();

  Matrix3XR<REAL> m_pos;
}; // class PointCloudRS

template<typename REAL, typename SIZE>
class DynamicPointCloudRS;

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

  void update_data();
  void update_glbuf();
  void update_shader(ShaderManager::ShaderType type);

  DynamicPointCloudRS<REAL,SIZE> *make_dynamic(REAL mass);

  void print() const { std::cerr << *m_pc << std::endl; }

protected:
  PointCloudRS<REAL, SIZE> *m_pc; // reference to the native mesh object

  // intermediate buffer between pc and glbuffer
  std::vector<GLfloat> m_vertices;
  bool                 m_insync;
};

// defaults
typedef PointCloudRS<double, unsigned int> PointCloud;
typedef GLPointCloudRS<double, unsigned int> GLPointCloud;

#endif // POINTCLOUD_H
