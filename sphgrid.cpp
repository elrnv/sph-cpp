#include <algorithm>
#include "settings.h"
#include "glpointcloud.h"
#include "gltext.h"
//#include "quantityprocessor.h"
#include "fluid.h"
#include "dynamicsmanager.h"
#include "sphgrid.h"
#define BOOST_CHRONO_HEADER_ONLY
#include <boost/chrono.hpp>

// From [Solenthaler and Pajarola 2008], an alternative density
// (number_density * mass) is used
/*
class CBVolume : public CBQ<CBVolume>
{
public:
  inline void init_kernel(float h) { m_kern.init(h); }
  inline void init_particle(Particle &p) 
  { p.dinv = 0.0f; }
  inline void fluid(Particle &p, FluidParticle &near_p)
  { }
  inline void bound(Particle &p, Particle &near_p)
  {
    p.dinv += this->m_kern[ p.pos - near_p.pos ];
  }
  inline void finish_particle(Particle &p)
  {
    p.dinv = 1.0f/(p.dinv * this->m_kern.coef);
  }
private:
  CubicSplineKernel m_kern;
};
*/

// SPHGrid stuff


SPHGrid::SPHGrid(
    AlignedBox3f &box,
    DynamicsManager &dynman ) 
  : m_dynman(dynman)
  , m_h(2e-3f) // minimum possible cell size
  , m_hinv(500.0f)
  , m_bmin(box.corner(AlignedBox3f::BottomLeftFloor))// smallest boundary corner
  , m_bmax(box.corner(AlignedBox3f::TopRightCeil))   // largest boundary corner
{ }


SPHGrid::~SPHGrid() { }

void
SPHGrid::init()
{
  if (m_dynman.get_num_fluids() < 1)
    return;

  // set cell size to be the maximum of the kernel support over all fluids;
  m_h = std::max(float(m_dynman.get_max_radius()), m_h);

  m_hinv = 1.0f/m_h;
  glprintf_tr("cell size: %f\n", m_h);

  // determine the number of cells needed
  m_gridsize = {{
    static_cast<Index>((m_bmax[0] - m_bmin[0])*m_hinv),
    static_cast<Index>((m_bmax[1] - m_bmin[1])*m_hinv),
    static_cast<Index>((m_bmax[2] - m_bmin[2])*m_hinv) }};

  glprintf_tr("grid size: %dx%dx%d = %d\n", 
      m_gridsize[0], m_gridsize[1], m_gridsize[2],
      m_gridsize[0]*m_gridsize[1]*m_gridsize[2]);

  m_grid.resize( boost::extents[m_gridsize[0]][m_gridsize[1]][m_gridsize[2]] );
  
  m_dynman.populate_sph_grid(*this); // populate grid with particles

  //CBVolume proc;
  //proc.init_kernel(m_h);
  //m_bound_volume_proc = &proc;
  //compute_bound_quantity< CBVolume >();
  //m_bound_volume_proc = NULL;
}

void
SPHGrid::update_grid()
{
  clear_fluid_particles();
  // populate grid with fluid particles
  m_dynman.populate_sph_grid_with_fluids(*this);
}

inline void
SPHGrid::clear_fluid_particles()
{
  for (Size i = 0; i < m_gridsize[0]; ++i)
    for (Size j = 0; j < m_gridsize[1]; ++j)
      for (Size k = 0; k < m_gridsize[2]; ++k)
        m_grid[i][j][k].clear();
}


template<int PT>
inline void
SPHGrid::compute_density()
{
  FluidDataVecT<PT> &fldatavec = m_dynman.get_fluiddatas<PT>();
  if (!fldatavec.size())
    return;

#ifdef REPORT_DENSITY_VARIATION
  Real max_var[fldatavec.size()];
  Real avg_var[fldatavec.size()];
  
  int i = 0;
  for ( auto &fldata : fldatavec )
  {
    max_var[i] = avg_var[i] = 0.0f; // initialize
    fldata.m_max_var = &max_var[i];
    fldata.m_avg_var = &avg_var[i];
    i++;
  }
#endif

  compute_quantity<Density,PT>();

#ifdef REPORT_DENSITY_VARIATION
  FluidVec &fluids = m_dynman.get_fluids();
  i = 0;
  for (auto &fldata : fldatavec )
  {
    avg_var[i] = avg_var[i]/fluids[fldata.flidx].get_num_vertices();
    qDebug("Fluid %d:  max: %.0f, %.1f percent;    avg: %.0f, %.1f percent", i,
        max_var[i], 100.00f*max_var[i]/fluids[fldata.flidx].get_rest_density(), 
        avg_var[i], 100.00f*avg_var[i]/fluids[fldata.flidx].get_rest_density());
    i++;
  }
#endif
}

template<int PT>
inline void 
SPHGrid::compute_accel()
{
  FluidDataVecT<PT> &fldatavec = m_dynman.get_fluiddatas<PT>();
  if (!fldatavec.size())
    return;

  FluidVec &fluids = m_dynman.get_fluids();
  for ( auto &fldata : fldatavec )
    fluids[fldata.flidx].reset_accel();   
  // now we may assume all accelerations are zero

  compute_quantity<Accel,PT>();
}

// instance the quantity computation functions above for each type
#define INSTANCE_COMPUTE_FUNC_TEMPLATE(z, PT, func) \
  template void SPHGrid::func<PT>();
BOOST_PP_REPEAT(NUMFLUIDTYPES, INSTANCE_COMPUTE_FUNC_TEMPLATE, compute_density)
BOOST_PP_REPEAT(NUMFLUIDTYPES, INSTANCE_COMPUTE_FUNC_TEMPLATE, compute_accel)
