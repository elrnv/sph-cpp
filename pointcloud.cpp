#include <QDebug>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLBuffer>
#include <assimp/scene.h>
#include <thread>
#include <limits>
#include <queue>
#include "pointcloud.h"

// PointCloud stuff
PointCloud::PointCloud(const aiMesh *mesh, Index matidx)
  : Primitive(matidx)
  , m_mindist(-1.f)
  , m_stalepos(true)
{
  // First copy all the vertices
  if (!mesh->HasPositions())
    return;

  aiVector3D *verts = mesh->mVertices;
  Size num_verts = mesh->mNumVertices;
  m_pos.resize(NoChange, num_verts);

  for (Size i = 0; i < num_verts; ++i)
    m_pos.col(i) << verts[i].x, verts[i].y, verts[i].z;

  prepare_vispos(); // copy data for visualizer
}

// Default constructor used for ghost and boundary particles
PointCloud::PointCloud(const Matrix3XT<Real> &pos)
  : Primitive(0)
  , m_mindist(-1.f)
  , m_pos(pos)
  , m_stalepos(true)
{
  prepare_vispos(); // copy data for visualizer
}

PointCloud::PointCloud(const PointCloud &orig)
  : Primitive(orig)
  , m_mindist(orig.m_mindist)
  , m_pos(orig.m_pos)
  , m_vispos(orig.m_vispos)
  , m_stalepos(true)
{ }

PointCloud::~PointCloud()
{ }

inline void
PointCloud::transform_in_place(const Affine3f &trans)
{
  m_pos = trans.template cast<Real>() * m_pos;
}

inline AlignedBox3f &
PointCloud::compute_bbox()
{
  m_bbox.setEmpty();
  Size num_verts = m_pos.cols();
  for (Size i = 0; i < num_verts; ++i)
    m_bbox.extend(Vector3d(m_pos.col(i)).cast<float>());
  return m_bbox;
}

Real
PointCloud::compute_mindist_brute()
{
  if (m_mindist != -1)
    return m_mindist;

  float dist2 = std::numeric_limits<float>::infinity();
  Size num_verts = m_pos.cols();
  
  for (Size i = 0; i < num_verts-1; ++i)
  {
    for (Size j = i+1; j < num_verts; ++j)
    {
      float curdist2 = (m_pos.col(i) - m_pos.col(j)).squaredNorm();
      if ( curdist2 < dist2 )
        dist2 = curdist2;
    }
  }
  return m_mindist = Real(std::sqrt(dist2));
}

Real
PointCloud::compute_mindist()
{
  if (m_mindist != -1)
    return m_mindist;

  float dist2 = std::numeric_limits<float>::infinity();
  
  Size num_verts = m_pos.cols();
  std::vector<Vector3f> S(num_verts);

  for (Size i = 0; i < num_verts; ++i)
    S[i] = Vector3f(m_pos.col(i).template cast<float>());

  std::sort(S.begin(), S.end(), 
      [](const Vector3f &p, const Vector3f &q) { return p[0] < q[0]; });

  std::deque<Vector3f> L;

  for ( Size i = 0; i < num_verts; ++i)
  {
    const Vector3f &q = S[i];

    while ( !L.empty() )
    {
      float dx = L.front()[0] - q[0];
      if ( dx*dx > dist2 )
        L.pop_front();
      else
        break;
    }
    for ( const Vector3f &p : L )
    {
      float curdist2 = (p - q).squaredNorm();
      if ( curdist2 < dist2 )
        dist2 = curdist2;
    }
    L.push_back(q);
  }

  return m_mindist = Real(std::sqrt(dist2));
}


// Visualization stuff
// Called from the dynamics thread
void
PointCloud::prepare_vispos()
{
  if (m_stalepos)
  {
    m_vispos = m_pos.template cast<float>();
    m_stalepos = false;
  }
}

// String representation
std::ostream&
operator<<(std::ostream& out, const PointCloud& pc)
{
  Size num_verts = pc.get_num_vertices();
  out << "pc({ ";
  for (Size i = 0; i < num_verts; ++i)
  {
    if (i > 0) out << ",\n     ";
    Vector3d v = pc.m_pos.col(i);
    out << v[0] << " " << v[1] << " " << v[2];
  }
  out << "}";
  return out;
}
