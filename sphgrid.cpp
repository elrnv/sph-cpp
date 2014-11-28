#include <ctime>
#include <algorithm>
#include <limits>
#include <fstream>
#include <iomanip>
#include "sph.h"
#include "glpointcloud.h"
#include "gltext.h"
//#include "quantityprocessor.h"
#include "fluid.h"
#include "settings.h"
#define BOOST_CHRONO_HEADER_ONLY
#include <boost/chrono.hpp>

//#define REPORT_DENSITY_VARIATION

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
  m_h = std::max(m_dynman.get_max_radius(), m_h);

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

inline void
SPHGrid::update_grid()
{
  clear_fluid_data();
  // populate grid with fluid particles
  m_dynman.populate_sph_grid_with_fluids(*this);
}

inline void
SPHGrid::clear_fluid_data()
{
  for (Size i = 0; i < m_gridsize[0]; ++i)
    for (Size j = 0; j < m_gridsize[1]; ++j)
      for (Size k = 0; k < m_gridsize[2]; ++k)
        m_grid[i][j][k].clear();
}

// variadic base case
template<int F, typename ParticleType>
inline void
SPHGrid::interact_with_neigh_cell( ParticleType &p, Cell &cell )
{ }

// variadic induction
template<int F, typename ParticleType, int NPT, int... NPTs>
inline void
SPHGrid::interact_with_neigh_cell( ParticleType &p, Cell &cell )
{
  auto &neigh_pvec = cell.template get_pvec<NPT>();
  for ( auto &neigh_p : neigh_pvec )
    p.neigh<F>(neigh_p); // process neighbouring particles

  interact_with_neigh_cell<F,ParticleType, NPTs...>(p,cell);
}

// variadic base case
template<int F>
inline void
SPHGrid::compute_quantity_in_cell( Size i,  Size j,  Size k, 
                                   Size nx, Size ny, Size nz )
{ }

// variadic induction
template<int F, int PT, int... PTs>
inline void
SPHGrid::compute_quantity_in_cell( Size i, Size j, Size k,
                                   Size nx, Size ny, Size nz )
{
  auto &pvec = m_grid[i][j][k].template get_pvec<PT>();
  if (pvec.empty())
    return;

  IndexRange xrange = range3(i,nx);
  IndexRange yrange = range3(j,ny);
  IndexRange zrange = range3(k,nz);
  // TODO: why do we need the GridView? just iterate directly in the grid
  GridView neigh_view = m_grid[ boost::indices[xrange][yrange][zrange] ];

  for ( auto &p : pvec )  // prepare data
    p.init<F>();

  Size xrange_size = xrange.finish() - xrange.start();
  Size yrange_size = yrange.finish() - yrange.start();
  Size zrange_size = zrange.finish() - zrange.start();
  for (Size near_i = 0; near_i < xrange_size; ++near_i)
  {
    for (Size near_j = 0; near_j < yrange_size; ++near_j)
    {
      for (Size near_k = 0; near_k < zrange_size; ++near_k)
      {
        Cell &cell = neigh_view[near_i][near_j][near_k];
        for ( auto &p : pvec )
          interact_with_neigh_cell<F,PT,decltype(p),ALL_PARTICLE_TYPES>(p,cell);
      }
    }
  }

  for ( auto &p : pvec )  // finalize data
    p.finish<F>();

  compute_quantity_in_cell<F,PTs...>(i,j,k,nx,ny,nz);
}

template<typename F, int... PT>
inline void
SPHGrid::compute_quantity()
{
//  clock_t s = clock();
  Size nx = m_gridsize[0];
  Size ny = m_gridsize[1];
  Size nz = m_gridsize[2];
  for (Size i = 0; i < nx; ++i)
  {
    for (Size j = 0; j < ny; ++j)
    {
      for (Size k = 0; k < nz; ++k)
      {
        compute_quantity_in_cell<F,PT...>(i,j,k,nx,ny,nz);
      } // for k
    } // for j
  } // for i

  //qDebug() << "average time" << (float(clock() - s) / CLOCKS_PER_SEC);
}

#if 0
template<int FT>
inline void FTiter<FT>::update_density(SPHGrid &g,float timestep)
{
  if (g.m_fluids[FT].size())
  {
    for (auto &fl : g.m_fluids[FT])
      fl->template cast<FT>()->m_fluid_density_update_proc.init( timestep ); 

    g.template compute_density_updateT<FT>();
  }
  FTiter<FT+1>::update_density(g,timestep); // recurse
}
#endif

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
  i = 0;
  for (auto &fldata : fldatavec )
  {
    avg_var[i] = avg_var[i]/fldata.fl.get_num_vertices();
    qDebug("Fluid %d:  max: %.0f, %.1f percent;    avg: %.0f, %.1f percent", i,
        max_var[i], 100.00f*max_var[i]/fldata.fl.get_rest_density(), 
        avg_var[i], 100.00f*avg_var[i]/fldata.fl.get_rest_density());
    i++;
  }
#endif
}

template<int PT, int PT2, int... PTs>
inline void
SPHGrid::compute_density()
{
  compute_density<PT>();       // do one
  compute_density<PT2, PTs>(); // recurse on the rest
}

template<int PT>
inline void 
SPHGrid::compute_accel()
{
  FluidDataVecT<PT> &fldatavec = m_dynman.get_fluiddatas<PT>();
  if (!fldatavec.size())
    return;

  for ( auto &fldata : fldatavec )
    fldata.fl.reset_accel();   
  // now we may assume all accelerations are zero

  compute_quantity<Accel,PT>();
}

template<int PT, int PT2, int... PTs>
inline void 
SPHGrid::compute_accel()
{
  compute_accel<PT>();
  compute_accel<PT2, PTs>();
}

#if 0
inline void SPHGrid::jacobi_pressure_solve(float dt, float factor)
{
  if (!m_fluids[ICS13].size())
    return;
  
  for (auto &fl : m_fluids[ICS13])
  {
    //fl->template cast<ICS13>()->m_fluid_prepare_jacobi_proc.init( dt );
    //fl->template cast<ICS13>()->m_fluid_jacobi_solve1_proc.init( dt );
    //fl->template cast<ICS13>()->m_fluid_jacobi_solve2_proc.init( dt, fl->get_avg_density(), fl->get_avg_pressure() );
    //fl->template cast<ICS13>()->reset_extern_accel();     // now may assume all accelerations are zero
  }

  //compute_fluid_quantity< CFPrepareJacobiT<ICS13>, ICS13 >();

  bool proceed = false;
  int iter = 0;
  do
  {
    //compute_fluid_quantity< CFJacobiSolveFirstT<ICS13>, ICS13 >();
    //compute_fluid_quantity< CFJacobiSolveSecondT<ICS13>, ICS13 >();
    proceed = false;
    for (auto &fl : m_fluids[ICS13])
    {
      proceed |= std::abs(fl->get_avg_density()/fl->get_num_vertices() -
          fl->get_rest_density()) > 1;
      qDebug() << "avg density  = " << fl->get_avg_density()/fl->get_num_vertices();
      qDebug() << "avg pressure = " << fl->get_avg_pressure()/fl->get_num_vertices();
      fl->get_avg_density() = 0.0f;
      fl->get_avg_pressure() = 0.0f;
    }
  } while(proceed);

  //compute_fluid_quantity< CFPressureAccelT<ICS13>, ICS13 >();

  for (auto &fl : m_fluids[ICS13])
    fl->get_vel() = fl->get_vel() + factor*dt*fl->get_extern_accel();
}
#endif

