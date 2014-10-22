#ifndef POINTCLOUD_H
#define POINTCLOUD_H

#include <iostream>
#include <assimp/mesh.h>
#include <boost/shared_ptr.hpp>
#include "primitive.h"
#include "dynparams.h"

// Partially resolve matrix template for convenience

class PointCloud;

std::ostream& operator<<(std::ostream& out, const PointCloud& pc);

// A cloud of points
class PointCloud : public Primitive 
{
public:
  explicit PointCloud(const aiMesh *pc);
  ~PointCloud();

  inline Real get_radius() { return 0.5*compute_mindist(); }
  virtual Real get_halo_radius() { return get_radius(); }
  inline Size get_num_vertices() const { return m_pos.cols(); }
  inline Real *pos_at(Size i)    { return m_pos.data() + i*3; }
  inline Matrix3XR<Real> &get_pos() { return m_pos; }
  inline void set_pos(const Matrix3XR<Real> &pos) { m_pos = pos; }
  inline bool is_pointcloud() const { return true; }

  void transform_in_place(const AffineCompact3f &trans);
  AlignedBox3f &compute_bbox();
  Real compute_mindist();
  Real compute_mindist_brute();

  friend std::ostream& operator<<(std::ostream& out, const PointCloud& pc);

protected:
  Real m_mindist;
  Matrix3XR<Real> m_pos;
}; // class PointCloud

typedef boost::shared_ptr< PointCloud > PointCloudPtr;

#endif // POINTCLOUD_H
