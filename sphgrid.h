#ifndef SPHGRID_H
#define SPHGRID_H

#include <vector>
#include <deque>
#include <boost/multi_array.hpp>
#include "types.h"
#include "pointcloud.h"
#include "fluid.h"
#include "dynamicsmanager.h"
#include "particle.h"

//class CBVolume;

// Grid structure used to optimize computing particle properties using kernels
class SPHGrid
{
public:

  // grid-cell data definitions
  template<int PT>
  using ParticleVecT = std::vector< ParticleT<PT> >;

  struct Cell
  {
    // The following macros define the data and their respective getters
    // of particle vectors for each particle type (PT) defined in types.h

    template<int PT>
    inline void push_particle(FluidDataT<PT> &fldata, float color, Size vtxidx)
    { 
      Vector3R<Real> pos(fldata.fl.get_pos().col(vtxidx));
      Vector3R<Real> vel(fldata.fl.get_vel().col(vtxidx));
      get_pvec<PT>().push_back(
          FluidParticleT<PT>( pos, vel, 
            fldata.fl.accel_at(vtxidx), 
            fldata.fl.extern_accel_at(vtxidx), 
            fldata.fl.dinv_at(vtxidx), 
            color, fldata));
    }

    inline void push_particle(BoundaryPC &bnd, Size vtxidx)
    {
      Vector3R<Real> pos(bnd->get_pos().col(vtxidx));
      get_pvec<STATIC>().push_back(ParticleT<STATIC>(pos,
            bnd.get_kernel_radius()));
    }

    inline void clear() 
    {
#define CLEAR_PVEC_TEMPLATE(z, PT, _) m_pvec_##PT.clear();
      BOOST_PP_REPEAT(NUMTYPES, CLEAR_PVEC_TEMPLATE, _)
    }
    
    // particle getters
#define PARTICLE_VEC_GETTER(z, PT, _) \
    inline DynamicParticleVecT<PT> & \
    get_pvec<PT>() { return m_pvec_##PT; }
#define PARTICLE_VEC_CONST_GETTER(z, PT, _) \
    inline const DynamicParticleVecT<PT> & \
    get_pvec<PT>() const { return m_pvec_##PT; }

    BOOST_PP_REPEAT(NUMTYPES, PARTICLE_VEC_GETTER, _)
    BOOST_PP_REPEAT(NUMTYPES, PARTICLE_VEC_CONST_GETTER, _)

#undef PARTICLE_VEC_GETTER
#undef PARTICLE_VEC_CONST_GETTER

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
  SPHGrid(const Vector3f &bmin, const Vector3f &bmax, DynamicsManager &dynman);
  ~SPHGrid();
  
  void init();
  void update_grid();
  void clear_fluid_data();

  inline Index clamp(Index d, Index min, Index max)
  {
    return std::max(std::min(d, max), min);
  }

  inline Array3Index get_voxel_index(const Vector3R<Real> &pos)
  {
    return {{ 
      clamp(static_cast<Index>(m_hinv*(pos[0]-m_bmin[0])), 0, m_gridsize[0]-1),
      clamp(static_cast<Index>(m_hinv*(pos[1]-m_bmin[1])), 0, m_gridsize[1]-1),
      clamp(static_cast<Index>(m_hinv*(pos[2]-m_bmin[2])), 0, m_gridsize[2]-1) }};
  }

  inline Cell &get_cell_at(const Array3Index &idx) { return m_grid(idx); }

  /*
  template <typename ProcessPairFunc>
  typename std::enable_if< std::is_same< ProcessPairFunc, CBVolume >::value,
           CBVolume >::type &get_proc()
           {
             return *m_bound_volume_proc;
           }

  template <typename ProcessPairFunc, typename ParticleType>
  typename std::enable_if< std::is_same< ParticleType, Particle >::value,
           ProcessPairFunc & >::type determine_proc(const ParticleType &p)
           {
             Q_UNUSED(p); return get_proc<ProcessPairFunc>();
           }
*/

  inline float get_cell_size() { return m_h; }
  inline Array3Index get_grid_size() { return m_gridsize; }
  inline const Vector3f &get_bmin() const { return m_bmin; }
  inline const Vector3f &get_bmax() const { return m_bmax; }

//  template <typename ProcessPairFunc, typename ParticleType>
//  typename std::enable_if< std::is_base_of< FluidParticle, ParticleType >::value,
//           ProcessPairFunc & >::type determine_proc(const FluidParticle &p)
//           {
//             return get_fluid<extract_fluid_type<ParticleType>::type>(p.id)
//               ->template get_proc<ProcessPairFunc>();
//           }

  // Low level quantity processing functions
  template<int F, int... PTs>
  void compute_quantity();

  //void jacobi_pressure_solve(float dt,float factor);

private: // member functions

  // two helper functions for compute_quantity
  template<int F, int... PTs>
  void compute_quantity_in_cell( Size i,  Size j,  Size k,
                                 Size nx, Size ny, Size nz );

  template<int F, typename ParticleType, int... NPTs>
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

#endif // SPHGRID_H
