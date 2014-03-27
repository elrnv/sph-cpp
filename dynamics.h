#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <vector>
#include <boost/multi_array.hpp>
#include "kernel.h"
#include "pointcloud.h"
#include "particle.h"
#include "quantityprocessor.h"

#define M_G 9.81f

// Compute routines used to compute SPH quantities

// From [Solenthaler and Pajarola 2008], an alternative density
// (number_density * mass) is used
template<typename REAL>
class CBVolumeR : public CBQPoly6R<REAL, CBVolumeR<REAL> >
{
public:
  inline void init_particle(ParticleR<REAL> &p) 
  { p.dinv = 0.0f; }
  inline void fluid(ParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  { }
  inline void bound(ParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    p.dinv += this->m_kern[ p.pos - near_p.pos ];
  }
  inline void finish_particle(ParticleR<REAL> &p)
  {
    p.dinv = 1.0f/(p.dinv * this->m_kern.coef);
  }
};

// forward declarations
template<typename REAL, typename SIZE>
class FluidRS;

template<typename REAL, typename SIZE>
class GLPointCloudRS;

// Grid structure used to optimize computing particle properties using kernels
template <typename REAL, typename SIZE>
class UniformGridRS
{
public:

  // grid-cell data definitions
  typedef std::vector< FluidParticleR<REAL> > DynamicParticles;
  typedef std::vector< ParticleR<REAL> > StaticParticles;

  struct Cell
  {
    DynamicParticles fluidvec; // dynamic fluid particles
    StaticParticles boundvec; // static boundary particles
    //DynamicParticles rigidvec; // dynamic rigid body object data

    template <typename ParticleType>
    typename std::enable_if<std::is_same< ParticleType, ParticleR<REAL> >::value,
             std::vector< ParticleType > >::type &get_vec() { return boundvec; }

    template <typename ParticleType>
    typename std::enable_if<std::is_same< ParticleType, FluidParticleR<REAL> >::value,
             std::vector< ParticleType > >::type &get_vec() { return fluidvec; }
  };

  typedef boost::multi_array< Cell, 3 > Array3;
  typedef typename Array3::index        Index;
  typedef boost::array<Index, 3>        Array3Index;
  typedef typename Array3::index_range  IndexRange;

  typedef typename Array3::template array_view<3>::type GridView;

  typedef std::vector< FluidRS<REAL, SIZE> * > FluidVec;

  // Constructors/Destructor
  UniformGridRS(const Vector3f &bmin, const Vector3f &bmax);
  ~UniformGridRS();
  
  void add_fluid(FluidRS<REAL,SIZE> *fl);
  void init();
  inline void update_grid();
  inline void populate_fluid_data();
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

#if 0
  void compute_initial_density();
#endif

  void compute_accel();
  void compute_pressure();
  void compute_density();
  void update_density(float timestep);

#if 0
  template<typename ProcessPairFunc>
  inline void compute_bound_quantity(ProcessPairFunc process);

  template<typename ProcessPairFunc>
  inline void compute_fluid_quantity(ProcessPairFunc process);
#endif

  template <typename ProcessPairFunc>
  typename std::enable_if< std::is_same< ProcessPairFunc, CBVolumeR<REAL> >::value,
           CBVolumeR<REAL> >::type &get_proc()
           {
             return *m_bound_volume_proc;
           }

  template <typename ProcessPairFunc>
  ProcessPairFunc &determine_proc(const ParticleR<REAL> &p)
           {
             Q_UNUSED(p); return get_proc<ProcessPairFunc>();
           }

  template <typename ProcessPairFunc>
  ProcessPairFunc &determine_proc(const FluidParticleR<REAL>  &p)
           {
             return m_fluids[p.id]->template get_proc<ProcessPairFunc>();
           }

  template<typename ProcessPairFunc, typename ParticleType>
  inline void compute_quantity();

  template<typename Func>
  inline void compute_bound_quantity()
  { compute_quantity<Func, ParticleR<REAL> >(); }

  template<typename Func>
  inline void compute_fluid_quantity()
  { compute_quantity<Func, FluidParticleR<REAL> >(); }

  // run dynamic simulation
  void run();

  // request stop which will be checked by the owner thread
  void request_stop() { m_stop_requested = true; }


private: // member functions
  // utility function used in the constructor to get a range of two elements
  // centered at x (or 2 if x is on the boundary)
  inline IndexRange range3(SIZE x, SIZE hi)
  {
    return IndexRange(x == 0 ? 0 : x-1, x == hi-1 ? hi : x+2 );
  }

private: // member variables
  FluidVec m_fluids; // array of simulated interacting fluids

  // array of cells containing xyzp (position and density) for each vertex
  Array3      m_grid;
  Array3Index m_gridsize;

  float m_h;    // cell size
  float m_hinv; // 1 / cell_size

  // the following two values define the axis aligned bounding box
  Vector3f m_bmin; // min boundary corner
  Vector3f m_bmax; // max boundary corner

  CBVolumeR<REAL> *m_bound_volume_proc;

  std::atomic<bool> m_stop_requested;
}; // class UniformGridRS

typedef UniformGridRS<double, unsigned int> UniformGrid;

#endif // DYNAMICS_H
