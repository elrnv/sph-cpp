#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include <QtGui/QOpenGLBuffer>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLShaderProgram>
#include "uniformbuffer.h"
#include "shadermanager.h"
#include "material.h"
#include "eigen.h"

class Primitive {
public:
  Primitive()
    : m_bbox(Vector3f(-1.0f,-1.0f,-1.0f), Vector3f(1.0f,1.0f,1.0f))
  { }
  virtual ~Primitive() { }
  virtual inline bool is_mesh()       const { return false; }
  virtual inline bool is_pointcloud() const { return false; }
  virtual inline bool is_dynamic()    const { return false; }
  AlignedBox3f &get_bbox() { return m_bbox; }

protected:
  AlignedBox3f m_bbox;
};

class Sphere : public Primitive 
{
public:
  Sphere()
  { }
};

class Cylinder : public Primitive
{
public:
  Cylinder()
    : m_normal_t(0,1,0)
  { }
  
private:
  // The default cylinder is hierarchical and is positioned at the origin with
  // radius 1 and height 2
  Vector3d m_normal_t; // normal of the top cap (bottom cap is the negative)

};

class Cone : public Primitive
{
public:
  Cone()
    : m_normal_b(0,-1,0) 
  { }

private:
  // The default cone is hierarchical with its base positioned 1 below the origin
  // and the tip distance 1 above the origin, that is pointing up
  Vector3d m_normal_b;
};



// A primitive representation for OpenGL applications (abstract)
template<typename SIZE>
class GLPrimitiveS : public QObject
{
public:
  explicit GLPrimitiveS(
      const PhongMaterial *mat,
      UniformBuffer &ubo,
      ShaderManager &shaderman)
    : m_ubo(ubo)
    , m_vao(this)
    , m_prog(NULL)
    , m_shaderman(shaderman)
    , m_pos(QOpenGLBuffer::VertexBuffer)
    , m_nml(QOpenGLBuffer::VertexBuffer)
    , m_idx(QOpenGLBuffer::IndexBuffer)
    , m_diffuse_color(mat->get_kd())
    , m_specular_color(mat->get_ks())
    , m_specular_power(mat->get_shininess()) { }
  ~GLPrimitiveS() { }

  QOpenGLVertexArrayObject &get_vao()     { return m_vao; }
  QOpenGLShaderProgram     *get_program() { return m_prog; }
  QOpenGLBuffer &get_pos() { return m_pos; }
  QOpenGLBuffer &get_nml() { return m_nml; }
  QOpenGLBuffer &get_idx() { return m_idx; }

  const Vector4f &get_diffuse()  const { return m_diffuse_color; }
  const Vector4f &get_specular() const { return m_specular_color; }
        double    get_specpow()  const { return m_specular_power; }

  virtual bool is_pointcloud() const { return false; }
  virtual bool is_mesh() const { return false; }

  virtual inline SIZE get_num_indices()  const = 0;
  virtual inline SIZE get_num_vertices() const = 0;

  virtual void update_shader(ShaderManager::ShaderType type) = 0;

  virtual void print() const = 0;

protected:
  UniformBuffer &m_ubo;
  QOpenGLVertexArrayObject m_vao; // vertex array object for this mesh
  QOpenGLShaderProgram *m_prog;  // glsl programs used for rendering
  ShaderManager &m_shaderman;

  QOpenGLBuffer m_pos; // position buffer
  QOpenGLBuffer m_nml; // normal buffer
  QOpenGLBuffer m_idx; // index buffer

  Vector4f m_diffuse_color;
  Vector4f m_specular_color;
  float    m_specular_power; // shininess
};

typedef GLPrimitiveS<unsigned int> GLPrimitive;

#endif // PRIMITIVE_H
