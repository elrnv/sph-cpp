#include <ctime>
#include <algorithm>
#include <limits>
#include <fstream>
#include <sstream>
#include "fluid.h"
#include "pointcloud.h"

BoundaryPC::BoundaryPC(const aiMesh *pc, Index matidx, DynParamsPtr params)
  : m_pc(pc)
  , m_params(params)
  , m_kernel_radius(get_kernel_radius())
{ }

// A default transparent boundary resembing a box covering the given sph grid
BoundaryPC::BoundaryPC(SPHGrid &grid, int particles_per_cell_length = 4)
  : m_pc(generate_grid_box_pc(grid, particles_per_cell_length))
  , m_params(NULL)
  , m_kernel_radius(grid.get_cell_size())
{ }

// helper function to the constructor above
inline Matrix3XR<Real> 
BoundaryPC::generate_grid_box_pc(SPHGrid &grid, int particles_per_cell_length)
{
  float h = grid.get_cell_size();
  SPHGrid::Array3Index &gridsize = grid.get_grid_size();
  const Vector3f &bmin = grid.get_bmin();
  const Vector3f &bmax = grid.get_bmax();

  // place boundary particles slightly away from the actual boundary
  float pad = 0.5*grid.get_cell_size();

  Size nx = particles_per_cell_length*(gridsize[0]+2);
  Size ny = particles_per_cell_length*(gridsize[1]+2);
  Size nz = particles_per_cell_length*(gridsize[2]+2);
  float incx = (bmax[0] - bmin[0] + 2*pad)/nx;
  float incy = (bmax[1] - bmin[1] + 2*pad)/ny;
  float incz = (bmax[2] - bmin[2] + 2*pad)/nz;

  Size i,j,k; // indices
  Matrix3XR<Real> pos; // result
  
  auto f = [&](Size i, Size j, Size k)
  {
    pos << Vector3R<Real>(bmin[0] - pad + i*incx, 
                          bmin[1] - pad + j*incy,
                          bmin[2] - pad + k*incz);
    //fprintf(stderr, "v %f %f %f\n", pos[0], pos[1], pos[2]);
  };

  k = 0;
  for (i = 0; i < nx; ++i)
    for (j = 0; j < ny; ++j)
      f(i,j,k);
  k = nz;
  for (i = nx; i > 0; --i)
    for (j = ny; j > 0; --j)
      f(i,j,k);

  j = 0;
  for (i = nx; i > 0; --i)
    for (k = nz; k > 0; --k)
      f(i,j,k);
  j = ny;
  for (i = 0; i < nx; ++i)
    for (k = 0; k < nz; ++k)
      f(i,j,k);

  i = 0;
  for (j = 0; j < ny; ++j)
    for (k = nz; k > 0; --k)
      f(i,j,k);
  i = nx;
  for (j = ny; j > 0; --j)
    for (k = 0; k < nz; ++k)
      f(i,j,k);

  // two points not filled
  f(0, ny, nz);
  f(nx, 0, 0);

  return pos;
}


BoundaryPC::~BoundaryPC()
{ }
