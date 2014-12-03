#ifndef BOUNDARYPC_H
#define BOUNDARYPC_H

#include <vector>
#include <assimp/scene.h>
#include "types.h"
#include "dynparams.h"
#include "pointcloud.h"

// BoundaryPC (Boundary Point Cloud)
// Def'n: this boundary is a static cloud of points, which passively interacts
//        with other objects.

// Forward declaration
class SPHGrid;

// A static cloud of points

class BoundaryPC
{
public:
  // dynamic point cloud from a regular updatable gl point cloud
  explicit BoundaryPC(const aiMesh *pc, Index matidx, RigidParamsPtr params);
  explicit BoundaryPC(SPHGrid &grid, int particles_per_cell_length = 4);
  ~BoundaryPC();

  Matrix3XT<Real> generate_grid_box_pc(
      SPHGrid &grid, int particles_per_cell_length);

  // interface for point cloud
  inline PointCloud &get_pc() { return m_pc; }
  void transform_in_place(const Affine3f &trans)
  {
    m_pc.transform_in_place(trans);
  }
  AlignedBox3f compute_bbox() { return m_pc.compute_bbox(); }
  inline Index get_material_idx() const { return m_pc.get_material_idx(); }
  inline void prepare_vispos() { m_pc.prepare_vispos(); }
  inline const Matrix3Xf &get_vispos() const { return m_pc.get_vispos(); }
  inline bool is_stalepos() { return m_pc.is_stalepos(); }
  inline void set_stalepos(bool sp) { m_pc.set_stalepos(sp); }

  inline Size get_num_vertices() const { return m_pc.get_num_vertices(); }
  inline Matrix3XT<Real> &get_pos()    { return m_pc.get_pos(); }
  inline const Matrix3XT<Real> &get_pos() const { return m_pc.get_pos(); }

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
  PointCloud     m_pc;   // The underlying dynamic cloud of points
  RigidParamsPtr m_params;
  Real           m_kernel_radius;
}; // class BoundaryPC

typedef std::vector< BoundaryPC > BoundaryPCVec;
typedef boost::shared_ptr< BoundaryPC > BoundaryPCPtr;

#endif // BOUNDARYPC_H
