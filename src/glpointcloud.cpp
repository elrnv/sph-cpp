#include <QDebug>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLBuffer>
#include <assimp/scene.h>
#include <limits>
#include "materialmanager.h"
#include "glpointcloud.h"
#include "fluid.h"

// GLPointCloud stuff

GLPointCloud::GLPointCloud(
    PointCloud &pc,
    bool dynamic,
    MaterialManager &matman,
    UniformBuffer &ubo,
    ShaderManager &shaderman,
    float halo_radius)
  : GLPrimitive(matman[pc.get_material_idx()], ubo)
  , m_pc(pc)
  , m_radius(pc.get_radius())
  , m_halo_radius(halo_radius)
  , m_isdynamic(dynamic)
{
  this->m_vao.create();
  this->m_vao.bind();

  this->m_pos.create();
  this->m_pos.setUsagePattern( 
      is_dynamic() ? QOpenGLBuffer::StreamDraw : QOpenGLBuffer::StaticDraw );

  const Matrix3XT<GLfloat> &vispos = m_pc.get_vispos();
  this->m_pos.bind();
  this->m_pos.allocate( vispos.data(), sizeof( GLfloat ) * vispos.size() );

  this->m_vao.release();

  update_shader(ShaderManager::ADDITIVE_PARTICLE, shaderman);
}

GLPointCloud::~GLPointCloud()
{ }

typedef PermutationMatrix<Dynamic,Dynamic,Size> SortMatrix;

//inline void GLPointCloud::sort_by_z(Matrix3XT<GLfloat> &pos)
inline void GLPointCloud::sort_by_z(const Matrix3XT<GLfloat> &pos,
                                    Matrix3XT<GLfloat> &outpos)
{
  // Sort all vertices by the z value
  clock_t s = clock();
  Size num_verts = get_num_vertices();
  VectorXT<Size> perm_vec(num_verts);
  for (Size i = 0; i < num_verts; ++i)
    perm_vec[i] = i;

  Size *perm_data = perm_vec.data();

  std::sort(perm_data, perm_data + num_verts,
      [pos](Size i, Size j) { return pos.col(i)[2] < pos.col(j)[2]; });
  clock_t sort_t = clock();

  outpos = outpos * SortMatrix(perm_vec);

  clock_t mult_t = clock();
  fprintf(stderr, "\rsort: %05.2e  mult: %05.2e", float(sort_t - s), float(mult_t - sort_t));
}

// private helper for the function below
inline void 
GLPointCloud::update_glbuf(const Matrix3XT<GLfloat> &vispos)
{
  this->m_pos.bind();
  this->m_pos.write( 0, vispos.data(), sizeof( GLfloat ) * vispos.size() );
  this->m_pos.release();
}

// Note: this is sensitive concurrent code
// currently it only works when one thread tries to read vispos and one thread
// writes to vispos
inline void 
GLPointCloud::update_glbuf_withsort(const AffineCompact3f &mvtrans,
                                    const AffineCompact3f &nmlmvtrans)
{
  //if (m_pc.is_stalepos())
  //  return;

  Matrix3XT<GLfloat> pos = m_pc.get_vispos();
  Matrix3XT<GLfloat> vispos = mvtrans * m_pc.get_vispos();
  sort_by_z(vispos, pos);
  update_glbuf(pos); // write to gl buffer

  m_pc.set_stalepos(true);
  (void) nmlmvtrans; // we dont map normals for point clouds
}
inline void 
GLPointCloud::update_glbuf_nosort()
{
  if (m_pc.is_stalepos())
    return;

  update_glbuf(m_pc.get_vispos()); // write to gl buffer
  m_pc.set_stalepos(true);
}

inline void 
GLPointCloud::update_shader(ShaderManager::ShaderType type,
                            ShaderManager &shaderman)
{
  if (this->m_prog)
  {
    this->m_prog->disableAttributeArray( "pos" );
  }

  if (type == ShaderManager::PARTICLE)
    this->m_prog = shaderman.get_particle_shader();
  else 
    this->m_prog = shaderman.get_additive_particle_shader();

  this->m_vao.bind();

  this->m_pos.bind();
  this->m_prog->enableAttributeArray( "pos" );
  this->m_prog->setAttributeBuffer( "pos", GL_FLOAT, 0, 3 );

  this->m_vao.release();

  this->m_ubo.bindToProg(this->m_prog->programId(), "Globals");
}
