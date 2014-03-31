#ifndef GLPRIMITIVE_H
#define GLPRIMITIVE_H

#include <QtGui/QOpenGLBuffer>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLShaderProgram>
#include <thread>
#include "uniformbuffer.h"
#include "shadermanager.h"
#include "material.h"
#include "eigen.h"

// A primitive representation for OpenGL applications (abstract)
template<typename SIZE>
class GLPrimitiveS : public QObject
{
public:
  explicit GLPrimitiveS(
      const Material *mat,
      UniformBuffer &ubo,
      ShaderManager &shaderman)
    : m_ubo(ubo)
    , m_vao(this)
    , m_prog(NULL)
    , m_shaderman(shaderman)
    , m_pos(QOpenGLBuffer::VertexBuffer)
    , m_nml(QOpenGLBuffer::VertexBuffer)
    , m_idx(QOpenGLBuffer::IndexBuffer)
    , m_specular_power(mat->get_shininess())
    , m_opacity(mat->get_opacity())
    , m_reflectivity(mat->get_reflectivity())
  { 
    const Vector3f &ka = mat->ka();
    const Vector3f &kd = mat->kd();
    const Vector3f &ks = mat->ks();
    m_ambient_color = QVector3D(ka[0], ka[1], ka[2]);
    m_diffuse_color = QVector3D(kd[0], kd[1], kd[2]);
    m_specular_color = QVector3D(ks[0], ks[1], ks[2]);
  }
  virtual ~GLPrimitiveS() { }

  QOpenGLVertexArrayObject &get_vao()     { return m_vao; }
  QOpenGLShaderProgram     *get_program() { return m_prog; }
  QOpenGLBuffer &get_pos() { return m_pos; }
  QOpenGLBuffer &get_nml() { return m_nml; }
  QOpenGLBuffer &get_idx() { return m_idx; }

  const QVector3D &get_ambient()  const { return m_ambient_color; }
  const QVector3D &get_diffuse()  const { return m_diffuse_color; }
  const QVector3D &get_specular() const { return m_specular_color; }
  float get_specpow()  const { return m_specular_power; }
  float get_opacity()  const { return m_opacity; }
  float get_reflectivity() const { return m_reflectivity; }

  virtual bool is_pointcloud() const { return false; }
  virtual bool is_mesh()       const { return false; }
  virtual bool is_dynamic()    const { return false; }

  virtual inline SIZE get_num_indices()  const = 0;
  virtual inline SIZE get_num_vertices() const = 0;

  // return the coordinates of the point closest to the camera (needed for sorting)
  virtual Vector3f get_closest_pt() const = 0;

  virtual void sort_by_depth(const AffineCompact3f &mvtrans) = 0;
  virtual void update_glbuf() = 0;
  virtual void update_shader(ShaderManager::ShaderType type) = 0;

  virtual void print() const = 0;

  std::mutex &get_lock() { return m_lock; }

protected:
  UniformBuffer &m_ubo;
  QOpenGLVertexArrayObject m_vao; // vertex array object for this mesh
  QOpenGLShaderProgram *m_prog;  // glsl programs used for rendering
  ShaderManager &m_shaderman;

  QOpenGLBuffer m_pos; // position buffer
  QOpenGLBuffer m_nml; // normal buffer
  QOpenGLBuffer m_idx; // index buffer

  QVector3D m_ambient_color;
  QVector3D m_diffuse_color;
  QVector3D m_specular_color;
  float    m_specular_power; // shininess
  float    m_opacity;
  float    m_reflectivity;

  std::mutex m_lock; // lock used for transferring dynamic data to GL
};

typedef GLPrimitiveS<unsigned int> GLPrimitive;

#endif // GLPRIMITIVE_H
