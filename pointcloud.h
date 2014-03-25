#ifndef POINTCLOUD_H
#define POINTCLOUD_H

#include <iostream>
#include <assimp/mesh.h>
#include "primitive.h"

// Partially resolve matrix template for convenience

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

  inline REAL get_radius() { return 0.5*compute_mindist(); }
  inline SIZE get_num_vertices() const { return m_pos.cols(); }
  inline Matrix3XR<REAL> &get_pos() { return m_pos; }
  inline bool is_pointcloud() const { return true; }

  inline void transform_in_place(const AffineCompact3f &trans);
  inline AlignedBox3f &compute_bbox();
  REAL compute_mindist();
  REAL compute_mindist_brute();

  friend std::ostream& operator<< <>(std::ostream& out, const PointCloudRS<REAL,SIZE>& pc);

protected:
  REAL m_mindist;
  Matrix3XR<REAL> m_pos;
}; // class PointCloudRS

// defaults
typedef PointCloudRS<double, unsigned int> PointCloud;

#endif // POINTCLOUD_H
