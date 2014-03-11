#include <QDebug>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLBuffer>
#include <assimp/scene.h>
#include "glpointcloud.h"
#include "dynamics.h"

// GLPointCloud stuff
template<typename REAL, typename SIZE>
GLPointCloudRS<REAL,SIZE>::GLPointCloudRS(
    PointCloudRS<REAL,SIZE> *pc,
    const Material *mat,
    UniformBuffer &ubo,
    ShaderManager &shaderman)
  : GLPrimitiveS<SIZE>(mat, ubo, shaderman)
  , m_pc(pc)
  , m_insync(true)
{
  // collect vertex attributes
  const Matrix3XR<REAL> &pts = pc->get_pos();
  GLfloat vertices[pts.size()];
  const REAL *pt_data = pts.data();
  for (SIZE i = 0; i < pts.size(); ++i)
    vertices[i] = GLfloat(pt_data[i]);

  this->m_vao.create();
  this->m_vao.bind();

  this->m_pos.create();
  this->m_pos.setUsagePattern( QOpenGLBuffer::StaticDraw );
  this->m_pos.bind();
  this->m_pos.allocate( vertices, sizeof( vertices ) );
  //qDebug() << "#" << vertices[0] << vertices[1] << vertices[2];

  update_shader(ShaderManager::PARTICLE);
}

template<typename REAL, typename SIZE>
GLPointCloudRS<REAL,SIZE>::~GLPointCloudRS()
{
}

template<typename REAL, typename SIZE>
void GLPointCloudRS<REAL,SIZE>::update_data()
{
  std::lock_guard<std::mutex> guard(this->m_lock);
  if (!m_insync)
    return;
  
  const Matrix3XR<REAL> &pts = m_pc->get_pos();
  m_vertices.resize(pts.size());
  const REAL *pt_data = pts.data();
  //qDebug() << "numpts written:" << pts.size();
  for (SIZE i = 0; i < pts.size(); ++i)
    m_vertices[i] = GLfloat(pt_data[i]);

  //qDebug() << "#" << m_vertices[0] << m_vertices[1] << m_vertices[2];

  m_insync = false;
}

template<typename REAL, typename SIZE>
void GLPointCloudRS<REAL,SIZE>::update_glbuf()
{
  std::lock_guard<std::mutex> guard(this->m_lock); // prevent others from reading buffers

  if (m_insync)
    return;

  //qDebug() << "num bytes read:" << sizeof(GLfloat) * m_vertices.size();
  this->m_pos.bind();
  this->m_pos.write( 0, m_vertices.data(), sizeof( GLfloat ) * m_vertices.size() );
  this->m_pos.release();
  //qDebug() << "-" << m_vertices[0] << m_vertices[1] << m_vertices[2];

  m_insync = true;
}


template<typename REAL, typename SIZE>
void GLPointCloudRS<REAL,SIZE>::update_shader(ShaderManager::ShaderType type)
{
  if (this->m_prog)
  {
    this->m_prog->disableAttributeArray( "pos" );
  }

  Q_UNUSED(type);
  //if (type == ShaderManager::PHONG)
  //  this->m_prog = this->m_shaderman.get_phong_shader();
  //else if (type == ShaderManager::PARTICLE)
    this->m_prog = this->m_shaderman.get_particle_shader();
  //else
  //  this->m_prog = this->m_shaderman.get_wireframe_shader();

  this->m_vao.bind();

  this->m_pos.bind();
  this->m_prog->enableAttributeArray( "pos" );
  this->m_prog->setAttributeBuffer( "pos", GL_FLOAT, 0, 3 );

  this->m_ubo.bindToProg(this->m_prog->programId(), "Globals");
}

template<typename REAL, typename SIZE>
DynamicPointCloudRS<REAL,SIZE> *GLPointCloudRS<REAL,SIZE>::make_dynamic(
    REAL density, REAL viscosity, REAL st)
{
  this->m_pos.setUsagePattern( QOpenGLBuffer::StreamDraw );

  DynamicPointCloudRS<REAL,SIZE> *dpc =
    new DynamicPointCloudRS<double, SIZE>(this, density, viscosity, st);

  delete m_pc; // delete the old point cloud
  m_pc = dpc;  // set new dynamic point cloud as a member
  return dpc;  // return new dynamic point cloud
}

template class GLPointCloudRS<double, unsigned int>;
