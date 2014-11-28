#ifndef POINTCLOUD_H
#define POINTCLOUD_H

#include <iostream>
#include <assimp/mesh.h>
#include <boost/shared_ptr.hpp>
#include "primitive.h"
#include "dynparams.h"
#include "eigen.h"

// Partially resolve matrix template for convenience

class PointCloud;

std::ostream& operator<<(std::ostream& out, const PointCloud& pc);

// A cloud of points
class PointCloud : public Primitive 
{
public:
  explicit PointCloud(const aiMesh *pc, Index matidx);
  explicit PointCloud(const Matrix3XR<Real> &pos);
  ~PointCloud();

  inline Real get_radius() { return 0.5*compute_mindist(); }
  virtual Real get_halo_radius() { return get_radius(); }
  inline Size get_num_vertices() const { return m_pos.cols(); }
  inline Real *pos_at(Size i)    { return m_pos.data() + i*3; }

  inline Matrix3XR<Real> &get_pos() { return m_pos; }
  inline const Matrix3XR<Real> &get_pos() const { return m_pos; }

  inline void set_pos(const Matrix3XR<Real> &pos) { m_pos = pos; }
  inline bool is_pointcloud() const { return true; }

  void transform_in_place(const AffineCompact3f &trans);
  AlignedBox3f &compute_bbox();
  Real compute_mindist();
  Real compute_mindist_brute();

  // manage position visualization data
  inline void prepare_vispos();
  template <typename T>
  inline Matrix3XR<T> &get_vispos() { return m_vispos.template cast<T>(); }
  inline bool is_stalepos() { return m_stalepos; }
  inline void set_stalepos(bool sp) { m_stalepos = sp; }

  friend std::ostream& operator<<(std::ostream& out, const PointCloud& pc);

protected:
  Real            m_mindist;      // minimum distance between a pair of points
  Matrix3XR<Real> m_pos;          // matrix of positions of all the points

  // the following data structures provide a mechanism for syncronizing computed
  // fluid data with its visual representative (like glpointcloud).
  Matrix3XR<Real>   m_vispos; // position data to be transferred to a visualizer
  std::atomic<bool> m_stalepos;  // true if data already collected
}; // class PointCloud

typedef boost::shared_ptr< PointCloud > PointCloudPtr;

#endif // POINTCLOUD_H
