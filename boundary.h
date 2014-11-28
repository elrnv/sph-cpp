#ifndef BOUNDARYPC_H
#define BOUNDARYPC_H

#include <vector>
#include <assimp/scene.h>
#include "types.h"
#include "dynparams.h"
#include "sphgrid.h"

// BoundaryPC (Boundary Point Cloud)
// Def'n: this boundary is a static cloud of points, which passively interacts
//        with other objects.

// Forward declaration
class PointCloud;

// A static cloud of points

class BoundaryPC
{
public:
  // dynamic point cloud from a regular updatable gl point cloud
  explicit BoundaryPC(const aiMesh *pc, Index matidx, DynParamsPtr params);
  explicit BoundaryPC(SPHGrid &grid, int particles_per_cell_length = 4);
  ~BoundaryPC();

  Matrix3XR<Real> generate_grid_box_pc(
      SPHGrid &grid, int particles_per_cell_length);

  inline Size get_num_vertices() const { return m_pc.get_num_vertices(); }
  inline Matrix3XR<Real> &get_pos()    { return m_pc.get_pos(); }
  inline const Matrix3XR<Real> &get_pos() const { return m_pc.get_pos(); }

  // kernel support radius
  inline Real get_kernel_radius()
  { 
    return m_params->kernel_inflation * m_pc.get_radius();
  }

  inline Real get_halo_radius()         { return get_kernel_radius(); }

  friend std::size_t hash_value( const BoundaryPC &bound ) 
  { 
    return hash_value(*(bound.m_params)); 
  }

protected:
  DynParamsPtr m_params;
  Real         m_kernel_radius;
  PointCloud   m_pc;   // The underlying dynamic cloud of points
}; // class BoundaryPC

typedef std::vector< BoundaryPC > BoundaryPCVec;
typedef boost::shared_ptr< BoundaryPC > BoundaryPCPtr;

#endif // BOUNDARYPC_H
