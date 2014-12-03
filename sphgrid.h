#ifndef SPHGRID_H
#define SPHGRID_H

#include <vector>
#include <deque>
#include <boost/multi_array.hpp>
#include "types.h"
#include "pointcloud.h"
#include "boundary.h"
#include "fluiddata.h"
#include "particle.h"
#include "dynamicsmanager.h"

// Grid structure used to optimize computing particle properties using kernels
class SPHGrid
{
public:
  struct Cell
  {
    // The following macros define the data and their respective getters
    // of particle vectors for each particle type (PT) defined in types.h

    template<int PT>
    inline void push_particle(FluidDataT<PT> &fldata, FluidVec &fluids,
                              float color, Size vtxidx)
    { 
      Fluid &fl = fluids[fldata.flidx];
      Vector3T<Real> pos(fl.get_pos().col(vtxidx));
      Vector3T<Real> vel(fl.get_vel().col(vtxidx));
      ParticleVecT<PT> &pvec = get_pvec<PT>();
      pvec.push_back(
          ParticleT<PT>( pos, vel, 
            fl.accel_at(vtxidx), 
            fl.extern_accel_at(vtxidx), 
            fl.dinv_at(vtxidx), 
            color, fldata));
    }

    void push_particle(BoundaryPC &bnd, Size vtxidx);

    inline void clear() 
    {
#define CLEAR_PVEC_TEMPLATE(z, PT, _) m_pvec_##PT.clear();
      BOOST_PP_REPEAT(NUMTYPES, CLEAR_PVEC_TEMPLATE, _)
    }
    
  template<int PT>
  ParticleVecT<PT> & get_pvec();

  template<int PT>
  const ParticleVecT<PT> & get_pvec() const;

    // data members
#define PARTICLE_VEC_MEMBER(z, PT, _) \
    ParticleVecT<PT> m_pvec_##PT;

    BOOST_PP_REPEAT(NUMTYPES, PARTICLE_VEC_MEMBER, _)
  };

  typedef boost::multi_array< Cell, 3 > Array3;
  typedef typename Array3::index        Index;
  typedef boost::array<Index, 3>        Array3Index;
  typedef typename Array3::index_range  IndexRange;

  typedef typename Array3::template array_view<3>::type GridView;

  // Constructors/Destructor
  SPHGrid(AlignedBox3f &box, DynamicsManager &dynman);
  ~SPHGrid();
  
  void init();
  void update_grid();
  void clear_fluid_particles();

  inline Index clamp(Index d, Index min, Index max)
  {
    return std::max(std::min(d, max), min);
  }

  inline Array3Index get_voxel_index(const Vector3T<Real> &pos)
  {
    return {{ 
      clamp(static_cast<Index>(m_hinv*(pos[0]-m_bmin[0])), 0, m_gridsize[0]-1),
      clamp(static_cast<Index>(m_hinv*(pos[1]-m_bmin[1])), 0, m_gridsize[1]-1),
      clamp(static_cast<Index>(m_hinv*(pos[2]-m_bmin[2])), 0, m_gridsize[2]-1) }};
  }

  inline Cell &get_cell_at(const Array3Index &idx) { return m_grid(idx); }

  inline float get_cell_size() { return m_h; }
  inline Array3Index get_grid_size() { return m_gridsize; }
  inline const Vector3f &get_bmin() const { return m_bmin; }
  inline const Vector3f &get_bmax() const { return m_bmax; }

  // Low level quantity processing functions
  template<int F, int... PTs>
  void compute_quantity();

  template<int PT>
  inline void reset_accel()
  {
    FluidDataVecT<PT> &fldatavec = m_dynman.get_fluiddatas<PT>();
    FluidVec &fluids = m_dynman.get_fluids();
    for ( auto &fldata : fldatavec )
      fluids[fldata.flidx].reset_accel();   
  }

  template<int PT, int PT2, int... PTs>
  inline void reset_accel()
  {
    reset_accel<PT>();
    reset_accel<PT2,PTs...>();
  }

  //template<int PT, int PT2, int... PTs>
  //void compute_accel();

  //template<int PT, int PT2, int... PTs>
  //void compute_density();

  //void jacobi_pressure_solve(float dt,float factor);

private: // member functions

  // two helper functions for compute_quantity
  template<int F> // base case
  void compute_quantity_in_cell( Size i,  Size j,  Size k,
                                 Size nx, Size ny, Size nz );
  template<int F, int PT, int... PTs>
  void compute_quantity_in_cell( Size i,  Size j,  Size k,
                                 Size nx, Size ny, Size nz );

  template<int F, typename ParticleType> // base case
  void interact_with_neigh_cell( ParticleType &p, Cell &cell );
  template<int F, typename ParticleType, int NPT, int... NPTs>
  void interact_with_neigh_cell( ParticleType &p, Cell &cell );

  // utility function used in the constructor to get a range of two elements
  // centered at x (or 2 if x is on the boundary)
  inline IndexRange range3(Size x, Size hi)
  {
    return IndexRange(x == 0 ? 0 : x-1, x == hi-1 ? hi : x+2 );
  }

private: // member variables
  // array of simulated interacting fluids
  DynamicsManager &m_dynman;

  // array of cells containing xyzp (position and density) for each vertex
  Array3      m_grid;
  Array3Index m_gridsize;

  float m_h;    // cell size
  float m_hinv; // 1 / cell_size

  // the following two values define the axis aligned bounding box
  Vector3f m_bmin; // min boundary corner
  Vector3f m_bmax; // max boundary corner

  //CBVolume *m_bound_volume_proc;

}; // class SPHGridRS

// particle getters
#define PARTICLE_VEC_GETTER(z, PT, _) \
  template<> \
  inline ParticleVecT<PT> & \
  SPHGrid::Cell::get_pvec<PT>() { return m_pvec_##PT; }
#define PARTICLE_VEC_CONST_GETTER(z, PT, _) \
  template<> \
  inline const ParticleVecT<PT> & \
  SPHGrid::Cell::get_pvec<PT>() const { return m_pvec_##PT; }

BOOST_PP_REPEAT(NUMTYPES, PARTICLE_VEC_GETTER, _)
BOOST_PP_REPEAT(NUMTYPES, PARTICLE_VEC_CONST_GETTER, _)

#undef PARTICLE_VEC_GETTER
#undef PARTICLE_VEC_CONST_GETTER

inline void 
SPHGrid::Cell::push_particle(BoundaryPC &bnd, Size vtxidx)
{
  const Vector3T<Real> &pos = bnd.get_pos().col(vtxidx);
  get_pvec<STATIC>().push_back(ParticleT<STATIC>(pos,
        bnd.get_kernel_radius()));
}

// Variadic template definitions

// variadic base case
template<int F, typename ParticleType>
inline void
SPHGrid::interact_with_neigh_cell( ParticleType &p, Cell &cell )
{ (void) p; (void) cell; }

// variadic induction
template<int F, typename ParticleType, int NPT, int... NPTs>
inline void
SPHGrid::interact_with_neigh_cell( ParticleType &p, Cell &cell )
{
  auto &neigh_pvec = cell.template get_pvec<NPT>();
  for ( auto &neigh_p : neigh_pvec )
    p.template neigh<F>(neigh_p); // process neighbouring particles

  interact_with_neigh_cell<F,ParticleType, NPTs...>(p,cell);
}

// variadic base case
template<int F>
inline void
SPHGrid::compute_quantity_in_cell( Size i,  Size j,  Size k, 
                                   Size nx, Size ny, Size nz )
{ (void) i; (void) j; (void) k; (void) nx; (void) ny; (void) nz; }

// variadic induction
template<int F, int PT, int... PTs>
inline void
SPHGrid::compute_quantity_in_cell( Size i, Size j, Size k,
                                   Size nx, Size ny, Size nz )
{
  auto &pvec = m_grid[i][j][k].template get_pvec<PT>();
  if (!pvec.empty())
  {
    IndexRange xrange = range3(i,nx);
    IndexRange yrange = range3(j,ny);
    IndexRange zrange = range3(k,nz);
    // TODO: why do we need the GridView? just iterate directly in the grid
    GridView neigh_view = m_grid[ boost::indices[xrange][yrange][zrange] ];

    for ( auto &p : pvec )  // prepare data
      p.template init<F>();

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
            interact_with_neigh_cell<F,decltype(p),ALL_PARTICLE_TYPES>(p,cell);
        }
      }
    }

    for ( auto &p : pvec )  // finalize data
      p.template finish<F>();
  }

  compute_quantity_in_cell<F,PTs...>(i,j,k,nx,ny,nz);
}

template<int F, int... PTs>
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
        compute_quantity_in_cell<F,PTs...>(i,j,k,nx,ny,nz);
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

#if 0
template<int PT, int PT2, int... PTs>
inline void
SPHGrid::compute_density()
{
  compute_density<PT>();       // do one
  compute_density<PT2, PTs...>(); // recurse on the rest
}

template<int PT, int PT2, int... PTs>
inline void 
SPHGrid::compute_accel()
{
  compute_accel<PT>();
  compute_accel<PT2, PTs...>();
}
#endif

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


#endif // SPHGRID_H
