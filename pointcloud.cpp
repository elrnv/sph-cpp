#include <QDebug>
#include <QtGui/QOpenGLVertexArrayObject>
#include <QtGui/QOpenGLBuffer>
#include <assimp/scene.h>
#include <thread>
#include <limits>
#include <queue>
#include "pointcloud.h"
#include "dynamics.h"

// PointCloud stuff

template<typename REAL, typename SIZE>
PointCloudRS<REAL,SIZE>::PointCloudRS(const aiMesh *mesh)
  : m_mindist(-1.f)
{
  // First copy all the vertices
  if (!mesh->HasPositions())
    return;

  aiVector3D *verts = mesh->mVertices;
  SIZE num_verts = mesh->mNumVertices;
  m_pos.resize(NoChange, num_verts);

  for (SIZE i = 0; i < num_verts; ++i)
    m_pos.col(i) << verts[i].x, verts[i].y, verts[i].z;
}

template<typename REAL, typename SIZE>
PointCloudRS<REAL,SIZE>::~PointCloudRS()
{
}

template<typename REAL, typename SIZE>
void PointCloudRS<REAL,SIZE>::transform(const AffineCompact3f &trans)
{
  // convert positions to canonical homogenized coordinates
  Translation3f T(trans.translation());
  Matrix3XR<REAL> Tmtx; // translation matrix
  Tmtx.resize(NoChange, m_pos.cols());
  Tmtx.row(0).fill(T.x());
  Tmtx.row(1).fill(T.y());
  Tmtx.row(2).fill(T.z());
  m_pos = trans.linear().template cast<REAL>() * m_pos + Tmtx;
}

template<typename REAL, typename SIZE>
AlignedBox3f &PointCloudRS<REAL,SIZE>::compute_bbox()
{
  m_bbox.setEmpty();
  SIZE num_verts = m_pos.cols();
  for (SIZE i = 0; i < num_verts; ++i)
    m_bbox.extend(Vector3d(m_pos.col(i)).cast<float>());
  return m_bbox;
}

template<typename REAL, typename SIZE>
REAL PointCloudRS<REAL,SIZE>::compute_mindist_brute()
{
  if (m_mindist != -1)
    return m_mindist;

  float dist2 = std::numeric_limits<float>::infinity();
  SIZE num_verts = m_pos.cols();
  
  for (SIZE i = 0; i < num_verts-1; ++i)
  {
    for (SIZE j = i+1; j < num_verts; ++j)
    {
      float curdist2 = (m_pos.col(i) - m_pos.col(j)).squaredNorm();
      if ( curdist2 < dist2 )
        dist2 = curdist2;
    }
  }
  return m_mindist = REAL(std::sqrt(dist2));
}

template<typename REAL, typename SIZE>
REAL PointCloudRS<REAL,SIZE>::compute_mindist()
{
  if (m_mindist != -1)
    return m_mindist;

  float dist2 = std::numeric_limits<float>::infinity();
  
  SIZE num_verts = m_pos.cols();
  std::vector<Vector3f> S(num_verts);

  for (SIZE i = 0; i < num_verts; ++i)
    S[i] = Vector3f(m_pos.col(i).template cast<float>());

  std::sort(S.begin(), S.end(), 
      [](const Vector3f &p, const Vector3f &q) { return p[0] < q[0]; });

  std::deque<Vector3f> L;

  for ( SIZE i = 0; i < num_verts; ++i)
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

  return m_mindist = REAL(std::sqrt(dist2));
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

template class PointCloudRS<double, unsigned int>;

template std::ostream& operator<<(std::ostream& out, const PointCloud& pc);
