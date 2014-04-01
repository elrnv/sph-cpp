#include <QDebug>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLBuffer>
#include <assimp/scene.h>
#include <limits>
#include "glpointcloud.h"
#include "dynamics.h"
#include "fluid.h"

// GLPointCloud stuff
template<typename REAL, typename SIZE>
GLPointCloudRS<REAL,SIZE>::GLPointCloudRS(
    PointCloudPtrRS<REAL,SIZE> pc,
    MaterialConstPtr mat,
    UniformBuffer &ubo,
    ShaderManager &shaderman)
  : GLPrimitiveS<SIZE>(mat, ubo, shaderman)
  , m_pc(pc)
  , m_vertices(3, get_num_vertices())
  , m_radius(pc->get_radius())
  , m_halo_radius(pc->get_halo_radius())
  , m_insync(true)
  , m_halos(false)
  , m_isdynamic(pc->is_dynamic())
{
  m_vertices = m_pc->get_pos().template cast<float>(); // copy position data

  this->m_vao.create();
  this->m_vao.bind();

  this->m_pos.create();
  this->m_pos.setUsagePattern( 
      is_dynamic() ? QOpenGLBuffer::StreamDraw : QOpenGLBuffer::StaticDraw );
  this->m_pos.bind();
  this->m_pos.allocate( m_vertices.data(), sizeof( GLfloat ) * m_vertices.size() );

  this->m_vao.release();

  if (is_dynamic())
  {
    FluidRS<REAL,SIZE> &fl = static_cast<FluidRS<REAL,SIZE> &>(*pc);
    fl.init(this);
  }

  update_shader(ShaderManager::PARTICLE);
}

template<typename REAL, typename SIZE>
GLPointCloudRS<REAL,SIZE>::~GLPointCloudRS()
{ }

template<typename REAL, typename SIZE>
void GLPointCloudRS<REAL,SIZE>::update_data()
{
  std::lock_guard<std::mutex> guard(this->m_lock);

  if (!m_insync)
    return;
  
  m_vertices = m_pc->get_pos().template cast<float>();

  m_insync = false;
}

template<typename REAL, typename SIZE>
Vector3f GLPointCloudRS<REAL,SIZE>::get_closest_pt() const
{
  return Vector3f(m_vertices.col(get_num_vertices()-1));
}

template<typename REAL, typename SIZE>
void GLPointCloudRS<REAL,SIZE>::sort_by_depth(const AffineCompact3f &mvtrans)
{
  std::lock_guard<std::mutex> guard(this->m_lock); // prevent others from reading buffers

  // Sort all vertices by the z value
  Matrix3XR<GLfloat> mvpos = mvtrans * m_vertices; // TODO: mem alloc expensive?
  SIZE num_verts = get_num_vertices();
  VectorXT<SIZE> perm_vec(num_verts);
  for (SIZE i = 0; i < num_verts; ++i)
    perm_vec[i] = i;

  SIZE *perm_data = perm_vec.data();

  std::sort(perm_data, perm_data + num_verts,
      [mvpos](SIZE i, SIZE j) { return mvpos.col(i)[2] < mvpos.col(j)[2]; });

  m_vertices = m_vertices * PermutationMatrix<Dynamic, Dynamic, SIZE>(perm_vec);

  m_insync = false;
}

template<typename REAL, typename SIZE>
void GLPointCloudRS<REAL,SIZE>::update_glbuf()
{
  // even though this routine doesn't access m_pc, another thread may be
  // writing to m_vertices so we still need to lock it
  std::lock_guard<std::mutex> guard(this->m_lock);

  if (m_insync)
    return;

  // write to gl buffer object
  this->m_pos.bind();
  this->m_pos.write( 0, m_vertices.data(), sizeof( GLfloat ) * m_vertices.size() );
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

  this->m_vao.release();

  this->m_ubo.bindToProg(this->m_prog->programId(), "Globals");
}

//template<typename REAL, typename SIZE>
//FluidRST<REAL,SIZE,FT> *GLPointCloudRS<REAL,SIZE>::make_dynamic(FluidParamsPtr params)
//{
//  this->m_pos.setUsagePattern( QOpenGLBuffer::StreamDraw );
//
//  FluidRST<REAL,SIZE> *dpc =
//    new FluidRST<REAL, SIZE>(m_pc, params);
//
//  delete m_pc; // delete the old point cloud
//  m_pc = dpc;  // set new dynamic point cloud as a member
//  return dpc;  // return new dynamic point cloud
//}
//
//template<typename REAL, typename SIZE>
//PointCloudRS<REAL,SIZE> *GLPointCloudRS<REAL,SIZE>::init_dynamics()
//{ 
//  if (m_pc->is_dynamic())
//  {
//    static_cast<FluidRST<REAL,SIZE,FT> *>(m_pc)->init(this);
//    return m_pc;
//  }
//
//  return NULL;
//}

template class GLPointCloudRS<double, unsigned int>;
