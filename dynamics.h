#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <deque>
#include <vector>
#include <boost/multi_array.hpp>
#include "particle.h"

#define M_G 9.81f

// Compute routine used to compute SPH quantities

template<typename REAL,typename SIZE>
class CBVolumeRS;
// forward declarations
template<typename REAL, typename SIZE>
class FluidRS;
template<typename REAL, typename SIZE, int FT>
class FluidRST;

template<typename REAL, typename SIZE>
using FluidPtrRS = boost::shared_ptr< FluidRS<REAL,SIZE> >;
template<typename REAL, typename SIZE, int FT>
using FluidPtrRST = boost::shared_ptr< FluidRST<REAL,SIZE,FT> >;

template<typename REAL, typename SIZE>
class GLPointCloudRS;

template <typename REAL, typename SIZE>
class UniformGridRS;

// template recursion mechanism to iterate through each fluid type
template<typename REAL, typename SIZE, int FT>
struct FTiter
{
  // high level processing (ones that should be called by the integrator)
  inline static void init_fluid_processors(UniformGridRS<REAL,SIZE> &);
  inline static void populate_fluid_data(UniformGridRS<REAL,SIZE> &g);
  inline static void update_density(UniformGridRS<REAL,SIZE> &, float);
  inline static void compute_density(UniformGridRS<REAL,SIZE> &);
  inline static void compute_accel(UniformGridRS<REAL,SIZE> &);
};

template<typename REAL, typename SIZE> // base case
struct FTiter<REAL, SIZE, NUMTYPES>
{
  inline static void init_fluid_processors(UniformGridRS<REAL,SIZE> &) { }
  inline static void populate_fluid_data(UniformGridRS<REAL,SIZE> &) { }
  inline static void update_density(UniformGridRS<REAL,SIZE> &,float) { }
  inline static void compute_density(UniformGridRS<REAL,SIZE> &) { }
  inline static void compute_accel(UniformGridRS<REAL,SIZE> &) { }
};

// Grid structure used to optimize computing particle properties using kernels
template <typename REAL, typename SIZE>
class UniformGridRS
{
public:

  // grid-cell data definitions
  typedef std::deque< FluidParticleRS<REAL,SIZE> > DynamicParticles;
  typedef std::deque< ParticleR<REAL> > StaticParticles;

  struct Cell
  {
    DynamicParticles fluidvec; // dynamic fluid particles
    StaticParticles boundvec; // static boundary particles
    //DynamicParticles rigidvec; // dynamic rigid body object data

    template <typename ParticleType>
    typename std::enable_if<std::is_same< ParticleType, ParticleR<REAL> >::value,
             StaticParticles &>::type get_vec() { return boundvec; }

    template <typename ParticleType>
    typename std::enable_if<std::is_same< FluidParticleRS<REAL,SIZE>, ParticleType >::value,
             DynamicParticles &>::type get_vec() { return fluidvec; }
  };

  typedef boost::multi_array< Cell, 3 > Array3;
  typedef typename Array3::index        Index;
  typedef boost::array<Index, 3>        Array3Index;
  typedef typename Array3::index_range  IndexRange;

  typedef typename Array3::template array_view<3>::type GridView;

  // Define a set of different types of fluids
  typedef std::vector< FluidPtrRS<REAL,SIZE> > FluidVec;

  // Constructors/Destructor
  UniformGridRS(const Vector3f &bmin, const Vector3f &bmax);
  ~UniformGridRS();
  
  inline unsigned int get_num_fluids() { return m_num_fluids; }
  inline void add_fluid(GLPointCloudRS<REAL,SIZE> &glpc) 
  {
    if (glpc.is_dynamic())
    {
      FluidPtrRS<REAL,SIZE> fl = boost::static_pointer_cast< FluidRS<REAL,SIZE> >(glpc.m_pc);
      m_fluids[fl->get_type()].push_back(fl);
      m_num_fluids += 1;
    }
  }

  void init();
  inline void update_grid();
  inline void populate_bound_data();
  inline void clear_fluid_data();

  inline Index clamp(Index d, Index min, Index max)
  {
    return std::max(std::min(d, max), min);
  }

  inline Array3Index get_voxel_index(const Vector3R<REAL> &pos)
  {
    return {{ 
      clamp(static_cast<Index>(m_hinv*(pos[0]-m_bmin[0])), 0, m_gridsize[0]-1),
      clamp(static_cast<Index>(m_hinv*(pos[1]-m_bmin[1])), 0, m_gridsize[1]-1),
      clamp(static_cast<Index>(m_hinv*(pos[2]-m_bmin[2])), 0, m_gridsize[2]-1) }};
  }

  template <typename ProcessPairFunc>
  typename std::enable_if< std::is_same< ProcessPairFunc, CBVolumeRS<REAL,SIZE> >::value,
           CBVolumeRS<REAL,SIZE> >::type &get_proc()
           {
             return *m_bound_volume_proc;
           }

  template <typename ProcessPairFunc, typename ParticleType, int FT>
  typename std::enable_if< std::is_same< ParticleType, ParticleR<REAL> >::value,
           ProcessPairFunc & >::type determine_proc(const ParticleType &p)
           {
             Q_UNUSED(p); return get_proc<ProcessPairFunc>();
           }

  template <typename ProcessPairFunc, typename ParticleType, int FT>
  typename std::enable_if< std::is_same< FluidParticleRS<REAL,SIZE>, ParticleType >::value,
           ProcessPairFunc & >::type determine_proc(const FluidParticleRS<REAL,SIZE> &p)
           {
             return static_cast<FluidRST<REAL,SIZE,FT>*>(p.fl)->template get_proc<ProcessPairFunc>();
           }

  // Low level quantity processing functions
  template<typename ProcessPairFunc, typename ParticleType, int FT>
  inline void compute_quantity();

  template<typename Func>
  inline void compute_bound_quantity()
  { compute_quantity<Func, ParticleR<REAL>, NOTFLUID >(); }

  template<typename Func, int FT>
  inline void compute_fluid_quantity()
  { compute_quantity<Func, FluidParticleRS<REAL, SIZE>, FT >(); }

  template<int FT>
  inline void compute_accelT();
  
  template<int FT>
  inline void compute_densityT();

  template<int FT>
  inline void compute_density_updateT();

#if 0
  void compute_initial_density();
#endif


  // run dynamic simulation
  void run();

  // request stop which will be checked by the owner thread
  void request_stop() { m_stop_requested = true; }

  friend FTiter<REAL,SIZE,NOTFLUID>;
  friend FTiter<REAL,SIZE,DEFAULT>;
  friend FTiter<REAL,SIZE,MCG03>;
  friend FTiter<REAL,SIZE,BT07>;
  friend FTiter<REAL,SIZE,AIAST12>;

private: // member functions

  // utility function used in the constructor to get a range of two elements
  // centered at x (or 2 if x is on the boundary)
  inline IndexRange range3(SIZE x, SIZE hi)
  {
    return IndexRange(x == 0 ? 0 : x-1, x == hi-1 ? hi : x+2 );
  }

private: // member variables
  // array of simulated interacting fluids
  FluidVec m_fluids[NUMTYPES];

  unsigned int m_num_fluids;

  // array of cells containing xyzp (position and density) for each vertex
  Array3      m_grid;
  Array3Index m_gridsize;

  float m_h;    // cell size
  float m_hinv; // 1 / cell_size

  // the following two values define the axis aligned bounding box
  Vector3f m_bmin; // min boundary corner
  Vector3f m_bmax; // max boundary corner

  CBVolumeRS<REAL,SIZE> *m_bound_volume_proc;

  std::atomic<bool> m_stop_requested;

}; // class UniformGridRS

typedef UniformGridRS<double, unsigned int> UniformGrid;

#endif // DYNAMICS_H
