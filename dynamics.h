#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <vector>
#include <deque>
#include <boost/multi_array.hpp>
#include "types.h"
#include "kernel.h"
#include "pointcloud.h"
#include "fluid.h"
#include "fluidmanager.h"
#include "glpointcloud.h"
#include "particle.h"

#define M_G 9.81f

// Compute routine used to compute SPH quantities

//class CBVolume;

// forward declarations

class UniformGrid;

// template recursion mechanism to iterate through each fluid type
template<int FT>
struct FTiter
{
  // high level processing (ones that should be called by the integrator)
  static void init_fluid_processors(UniformGrid &);
  static void update_density(UniformGrid &, float);
  static void compute_density(UniformGrid &);
  static void compute_accel(UniformGrid &);
  static void compute_pressure(UniformGrid &);
};

 // base case
template <>
struct FTiter<NUMTYPES>
{
  inline static void init_fluid_processors(UniformGrid &) { }
  inline static void update_density(UniformGrid &,float) { }
  inline static void compute_density(UniformGrid &) { }
  inline static void compute_accel(UniformGrid &) { }
  inline static void compute_pressure(UniformGrid &) { }
};

// Grid structure used to optimize computing particle properties using kernels

class UniformGrid
{
public:

  // grid-cell data definitions
  typedef std::vector< FluidParticle > DynamicParticles;
  typedef std::vector< Particle > StaticParticles;

  struct Cell
  {
    Index fluididx; // dynamic fluid particles
    Index boundidx; // static boundary particles
    //Index rigididx; // dynamic rigid body particles

    Index get_idx()
    {
      return ParticleType

    }

    template <typename ParticleType>
    typename std::enable_if<std::is_same< ParticleType, Particle >::value,
             std::vector< ParticleType > >::type &get_vec() { return boundvec; }

    template <typename ParticleType>
    typename std::enable_if<std::is_base_of< FluidParticle, ParticleType >::value,
             std::vector< FluidParticle > >::type &get_vec() { return fluidvec; }
  };

  typedef boost::multi_array< Cell, 3 > Array3;
  typedef typename Array3::index        Index;
  typedef boost::array<Index, 3>        Array3Index;
  typedef typename Array3::index_range  IndexRange;

  typedef typename Array3::template array_view<3>::type GridView;

  // Define a set of different types of fluids
  typedef std::deque< FluidPtr > FluidVec;

  // Constructors/Destructor
  UniformGrid(const Vector3f &bmin, const Vector3f &bmax);
  ~UniformGrid();
  
  void add_fluid(GLPointCloud &glpc);

  template<int FT>
  inline FluidT<FT> *
  get_fluid(unsigned int id) 
  { 
    return m_fluids[FT][id]->template cast<FT>();
  }

  void init();
  void update_grid();
  void populate_fluid_data();
  void populate_bound_data();
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

  template <typename ProcessPairFunc, typename ParticleType>
  typename std::enable_if< std::is_base_of< FluidParticle, ParticleType >::value,
           ProcessPairFunc & >::type determine_proc(const FluidParticle &p)
           {
             return get_fluid<extract_fluid_type<ParticleType>::type>(p.id)
               ->template get_proc<ProcessPairFunc>();
           }

  // Low level quantity processing functions
  template<typename ProcessPairFunc, typename ParticleType>
  void compute_quantity();

  template<typename Func>
  inline void compute_bound_quantity()
  { compute_quantity<Func, Particle >(); }

  template<typename Func, int FT>
  inline void compute_fluid_quantity()
  { compute_quantity<Func, FluidParticleT<FT> >(); }

  template<int FT> void compute_accelT();
  template<int FT> void compute_surface_normalT();
  template<int FT> void compute_surface_tensionT();

  template<int FT> void compute_densityT();

  template<int FT> void compute_density_updateT();

  void jacobi_pressure_solve(float dt,float factor);

  // run dynamic simulation
  void run();

  // execute one substep
  bool step(float dt, bool first_step, float *substep_t = NULL);

  // request stop which will be checked by the owner thread
  void request_stop() { m_stop_requested = true; }
  void toggle_pause() 
  { 
    std::unique_lock<std::mutex> locker(m_pause_lock);
    m_pause = !m_pause;
    if (!m_pause)
      m_pause_cv.notify_all();
  }
  void un_pause() 
  { 
    std::unique_lock<std::mutex> locker(m_pause_lock);
    if (m_pause)
    {
      m_pause = false;
      m_pause_cv.notify_all();
    }
  }

  bool check_and_write_hash();

  friend std::size_t hash_value( const UniformGrid &ug )
  {
    std::size_t seed = 0;
    for (int j = 0; j < NUMTYPES; ++j)
      for (auto &fl : ug.m_fluids[j])
        boost::hash_combine(seed, hash_value(*fl));
    return seed;
  }

  template <int FT>
  friend struct FTiter;

private: // member functions

  // utility function used in the constructor to get a range of two elements
  // centered at x (or 2 if x is on the boundary)
  inline IndexRange range3(Size x, Size hi)
  {
    return IndexRange(x == 0 ? 0 : x-1, x == hi-1 ? hi : x+2 );
  }

  // wrapper for fluid getter using the fluid manager

  template< int FT >
  inline get_fluids<FT> { return m_fluidman.get_fluids<FT>(); }

private: // member variables
  // array of simulated interacting fluids
  FluidManager m_fluidman;

  Size m_num_fluids;

  // array of cells containing xyzp (position and density) for each vertex
  Array3      m_grid;
  Array3Index m_gridsize;

  float m_h;    // cell size
  float m_hinv; // 1 / cell_size

  // the following two values define the axis aligned bounding box
  Vector3f m_bmin; // min boundary corner
  Vector3f m_bmax; // max boundary corner

  //CBVolume *m_bound_volume_proc;

  std::atomic<bool> m_stop_requested;

  std::mutex m_pause_lock;
  std::condition_variable m_pause_cv;
  std::atomic<bool> m_pause;

}; // class UniformGridRS

#endif // DYNAMICS_H
