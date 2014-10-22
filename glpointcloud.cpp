#include <QDebug>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLBuffer>
#include <assimp/scene.h>
#include <limits>
#include "glpointcloud.h"
#include "dynamics.h"
#include "fluid.h"

// GLPointCloud stuff

GLPointCloud::GLPointCloud(
    PointCloudPtr pc,
    bool dynamic,
    MaterialConstPtr mat,
    UniformBuffer &ubo,
    ShaderManager &shaderman)
  : GLPrimitive(mat, ubo, shaderman)
  , m_pc(pc)
  , m_vertices(3, get_num_vertices())
  , m_radius(pc->get_radius())
  , m_halo_radius(pc->get_halo_radius())
  , m_insync(true)
  , m_halos(false)
  , m_isdynamic(dynamic)
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
    Fluid &fl = static_cast<Fluid &>(*pc);
    fl.init(this);
  }

  update_shader(ShaderManager::ADDITIVE_PARTICLE);
}


GLPointCloud::~GLPointCloud()
{ }


void GLPointCloud::update_data()
{
  std::lock_guard<std::mutex> guard(this->m_lock);

  if (!m_insync)
    return;
  
  m_vertices = m_pc->get_pos().template cast<float>();

  m_insync = false;
}


void GLPointCloud::clear_cache()
{
  std::lock_guard<std::mutex> guard(this->m_lock);
  if (is_dynamic())
  {
    Fluid &fl = static_cast<Fluid &>(*m_pc);
    fl.clear_cache();
  }
}


Vector3f GLPointCloud::get_closest_pt() const
{
  return Vector3f(m_vertices.col(get_num_vertices()-1));
}


void GLPointCloud::sort_by_depth(const AffineCompact3f &mvtrans)
{
  std::lock_guard<std::mutex> guard(this->m_lock); // prevent others from reading buffers

  // Sort all vertices by the z value
  clock_t s = clock();
  Matrix3XR<GLfloat> mvpos = mvtrans * m_vertices; // TODO: mem alloc expensive?
  Size num_verts = get_num_vertices();
  VectorXT<Size> perm_vec(num_verts);
  for (Size i = 0; i < num_verts; ++i)
    perm_vec[i] = i;

  Size *perm_data = perm_vec.data();

  std::sort(perm_data, perm_data + num_verts,
      [mvpos](Size i, Size j) { return mvpos.col(i)[2] < mvpos.col(j)[2]; });
  clock_t sort_t = clock();

  m_vertices = m_vertices * PermutationMatrix<Dynamic, Dynamic, Size>(perm_vec);

  clock_t mult_t = clock();
  fprintf(stderr, "\rsort: %05.2e  mult: %05.2e", float(sort_t - s), float(mult_t - sort_t));

  m_insync = false;
}


void GLPointCloud::update_glbuf()
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



void GLPointCloud::update_shader(ShaderManager::ShaderType type)
{
  if (this->m_prog)
  {
    this->m_prog->disableAttributeArray( "pos" );
  }

  if (type == ShaderManager::PARTICLE)
    this->m_prog = this->m_shaderman.get_particle_shader();
  else 
    this->m_prog = this->m_shaderman.get_additive_particle_shader();

  this->m_vao.bind();

  this->m_pos.bind();
  this->m_prog->enableAttributeArray( "pos" );
  this->m_prog->setAttributeBuffer( "pos", GL_FLOAT, 0, 3 );

  this->m_vao.release();

  this->m_ubo.bindToProg(this->m_prog->programId(), "Globals");
}
