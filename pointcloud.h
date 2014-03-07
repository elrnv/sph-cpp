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

  Matrix3XR<REAL> &get_pos() { return m_pos; }
  inline SIZE get_num_vertices() const { return m_pos.cols(); }
  inline bool is_pointcloud() const { return true; }

  inline void transform(const AffineCompact3f &trans);

  friend std::ostream& operator<< <>(std::ostream& out, const PointCloudRS<REAL,SIZE>& pc);

protected:
  void compute_bbox();

  Matrix3XR<REAL> m_pos;
}; // class PointCloudRS

// defaults
typedef PointCloudRS<double, unsigned int> PointCloud;

#endif // POINTCLOUD_H
