#include <QDebug>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLBuffer>
#include <assimp/scene.h>
#include <thread>
#include "pointcloud.h"
#include "dynamics.h"

// PointCloud stuff

template<typename REAL, typename SIZE>
PointCloudRS<REAL,SIZE>::PointCloudRS(const aiMesh *mesh)
{
  // First copy all the vertices
  if (!mesh->HasPositions())
    return;

  aiVector3D *verts = mesh->mVertices;
  SIZE num_verts = mesh->mNumVertices;
  m_pos.resize(NoChange, num_verts);

  for (SIZE i = 0; i < num_verts; ++i)
    m_pos.col(i) << verts[i].x, verts[i].y, verts[i].z;

  compute_bbox(); // compute the bounded box if requested
}

template<typename REAL, typename SIZE>
PointCloudRS<REAL,SIZE>::~PointCloudRS()
{
}

template<typename REAL, typename SIZE>
void PointCloudRS<REAL,SIZE>::compute_bbox()
{
  m_bbox.setEmpty();
  SIZE num_verts = m_pos.cols();
  for (SIZE i = 0; i < num_verts; ++i)
    m_bbox.extend(Vector3d(m_pos.col(i)).cast<float>());
}

template<typename REAL, typename SIZE>
std::ostream& operator<<(std::ostream& out, const PointCloudRS<REAL,SIZE>& pc)
{
  SIZE num_verts = pc.get_num_vertices();
  out << "pc({ ";
  for (SIZE i = 0; i < num_verts; ++i)
  {
    if (i > 0) out << ",\n     ";
    Vector3d v = pc.m_pos.col(i);
    out << v[0] << " " << v[1] << " " << v[2];
  }
  out << "}";
  return out;
}



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
  std::cerr << "YAY" << std::endl;
  if (!m_insync)
    return;

  const Matrix3XR<REAL> &pts = m_pc->get_pos();
  m_vertices.resize(pts.size());
  const REAL *pt_data = pts.data();
  for (SIZE i = 0; i < pts.size(); ++i)
    m_vertices[i] = GLfloat(pt_data[i]);

  m_insync = false;
}

template<typename REAL, typename SIZE>
void GLPointCloudRS<REAL,SIZE>::update_glbuf()
{
  std::lock_guard<std::mutex> guard(this->m_lock); // prevent others from reading buffers

  if (m_insync)
    return;

  this->m_pos.bind();
  this->m_pos.write( 0, m_vertices.data(), sizeof( m_vertices ) );
  this->m_pos.release();

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
DynamicPointCloudRS<REAL,SIZE> *GLPointCloudRS<REAL,SIZE>::make_dynamic(REAL mass)
{
  this->m_pos.setUsagePattern( QOpenGLBuffer::StreamDraw );
  DynamicPointCloudRS<REAL,SIZE> *dpc = new DynamicPointCloudRS<double, SIZE>(*m_pc, mass, &GLPointCloudRS<REAL,SIZE>::update_data);
  delete m_pc; // delete the old point cloud
  m_pc = dpc;  // set new dynamic point cloud as a member
  return dpc;  // return new dynamic point cloud
}

template class PointCloudRS<double, unsigned int>;
template class GLPointCloudRS<double, unsigned int>;

template std::ostream& operator<<(std::ostream& out, const PointCloud& pc);
