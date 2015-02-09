#ifndef DYNAMICSMANAGER_H
#define DYNAMICSMANAGER_H

#include <ctime>
#include <algorithm>
#include <fstream>
#include <iomanip> // for setfill and setw

#include "types.h"
#include "settings.h"
#include "gltext.h"
#include "boundary.h"
#include "fluid.h"
#include "fluiddata.h"
#include "materialmanager.h"
#include "sphgrid.h"
#define BOOST_CHRONO_HEADER_ONLY
#include <boost/chrono.hpp>

typedef boost::chrono::process_real_cpu_clock real_clock;
typedef boost::chrono::nanoseconds ns_t;
typedef real_clock::time_point real_t;

// Helper template using two arguments from the data arguemnt
// where data is of the form (fluid, for loop body)

#define FLUID_TYPE_CASE_FUNC_TEMPLATE(z, PT, data) \
  case PT: \
    BOOST_PP_TUPLE_ELEM(2, 0, data)<PT>BOOST_PP_TUPLE_ELEM(2, 1, data); \
    break;

#define FLUID_TYPE_CASE_FUNC_CALL(func, params) \
  BOOST_PP_REPEAT(NUMFLUIDSPHTYPES, FLUID_TYPE_CASE_FUNC_TEMPLATE, (func, params))

class DynamicsManager
{
public:
  DynamicsManager() 
    : m_bbox(Vector3f(0.0f,0.0f,0.0f), Vector3f(0.0f,0.0f,0.0f))
    , m_grid(*this)
    , m_stop_requested(false)
    , m_pause(false)
  { }
  ~DynamicsManager() { }


  /// Manager interface

  // normalization routines used to fit the dynamic objects within a unit cube
  void normalize_models()
  {
    normalize_models(Vector3f(0.0f, 0.0f, 0.0f),Vector3f(0.0f, 0.0f, 0.0f));
  }
  void normalize_models(
      const Vector2f &ext_x,
      const Vector2f &ext_y,
      const Vector2f &ext_z)
  {
    normalize_models(Vector3f(ext_x[0], ext_y[0], ext_z[0]),
                     Vector3f(ext_x[1], ext_y[1], ext_z[1]));
  }

  void normalize_models(const Vector3f &ext_blf, const Vector3f &ext_trc)
  {
    compute_bbox();
    Vector3f blf( m_bbox.corner(Eigen::AlignedBox3f::BottomLeftFloor) );
    Vector3f trc( m_bbox.corner(Eigen::AlignedBox3f::TopRightCeil) );
    m_bbox.extend( blf - ext_blf );
    m_bbox.extend( trc + ext_trc );

    Vector3f sizevec = m_bbox.sizes();
    transform_models(Affine3f::Identity() * Scaling(2.0f/sizevec.maxCoeff()));
    Vector3f box_center = m_bbox.center();
    transform_models(Affine3f::Identity() * Translation3f(-box_center));
  }

  void transform_models(const Affine3f &trans)
  {
    for ( auto &fl : m_fluids )
      fl.transform_in_place(trans);
    for ( auto &bound : m_boundaries )
      bound.transform_in_place(trans);
  }

  AlignedBox3f &compute_bbox()
  {
    m_bbox.setEmpty();
    for ( auto &fl : m_fluids )
      m_bbox.extend(fl.compute_bbox());
    for ( auto &bound : m_boundaries )
      m_bbox.extend(bound.compute_bbox());
    return m_bbox;
  }

  inline void set_bbox(const AlignedBox3f &bbox) { m_bbox = bbox; }

  void add_dynamic_object(const aiMesh *mesh, Index mat_idx,
                          MaterialManager &matman, DynParamsPtr dynparams)
  {
    if ( mesh->mPrimitiveTypes & aiPrimitiveType_POINT )
    {
      if ( dynparams->type == DynParams::FLUID )
      {
        FluidParamsPtr fluidparams =
          boost::static_pointer_cast<FluidParams>(dynparams);
        switch (fluidparams->fluid_type)
        {
          FLUID_TYPE_CASE_FUNC_CALL(add_fluid,(mesh,mat_idx,matman,fluidparams))
          default: break;
        }
      }
    }
    // else ignore. TODO: extend this to include other dynamic objects
  }

  void add_static_object(const aiMesh *mesh, Index mat_idx)
  {
    // TODO: implement this for arbitrary boundaries
    (void) mesh; (void) mat_idx; // meanwhile suppress warnings
  }

  // Unit box boundary particles (transparent)
  // PRE: assume grid has been initialized
  void add_default_boundary(Index matidx)
  {
    m_boundaries.push_back(BoundaryPC(m_grid, matidx));
  }

  template<int PT>
  void add_fluid(const aiMesh *mesh, Index matidx, 
                 MaterialManager &matman, FluidParamsPtr fparams)
  {
    m_fluids.push_back(Fluid(mesh, matidx, matman, fparams));
  }

  inline void init_fluids(const AlignedBox3f &box)
  {
    for ( auto &fl : m_fluids )
      fl.init(box);
  }

  inline void init_sphgrid(const AlignedBox3f &box)
  {
    m_grid.init(box);
  }

  // PRE: assume that fluids were already initialized (init_fluids was called)
  inline void generate_fluiddatas()
  {
    generate_fluiddatas<ALL_FLUID_SPH_TYPES>();
  }

  inline Size get_num_fluids() { return m_fluids.size(); }
  inline Size get_num_boundaries() { return m_boundaries.size(); }

  inline void clear_fluid_sph_particles()
  {
    m_grid.clear_fluid_particles();
  }

  inline void clear_rigid_sph_particles()
  {
    m_grid.clear_rigid_particles();
  }

  inline void clear()
  {
    clear_fluid_sph_particles();
    clear_rigid_sph_particles();
    m_fluids.clear();
    m_boundaries.clear();
    clear_fluiddatas<ALL_FLUID_SPH_TYPES>();
  }

  template <int PT>
  inline void clear_fluiddatas()
  {
    get_fluiddatas<PT>().clear();
  }

  template <int PT, int PT2, int... PTs>
  inline void clear_fluiddatas()
  {
    clear_fluiddatas<PT>();
    clear_fluiddatas<PT2, PTs...>();
  }

  inline Real get_max_radius()
  {
    Real h = 0.0;
    for ( auto &fl : m_fluids )
    {
      if (h < fl.get_kernel_radius())
        h = fl.get_kernel_radius();
    }
    return h;
  }

  // prepare fluid data for visualization. This involves making a copy of
  // visualized data within each fluid, so that a visualizer can pick up a fresh
  // copy concurrently
  inline void prepare_vis_data()
  {
    for ( auto &fl : m_fluids )
      fl.prepare_vispos();
  }

  inline void clear_cache()
  {
    for ( auto &fl : m_fluids )
      fl.clear_cache();
  }
  inline void load_cached(unsigned int frame)
  {
    for ( auto &fl : m_fluids )
      fl.load_cached(frame);
  }

  // load all the cached frames.
  // return true if all the frames have been saved for all the fluids in the
  // scene
  inline bool load_saved_cache()
  {
    bool all_cached = true;
    for ( auto &fl : m_fluids )
      all_cached &= fl.load_saved_cache();
    return all_cached;
  }

  // for a frame to be cached, it must be cached for every fluid
  inline bool is_cached(unsigned int frame)
  {
    bool all_cached = true;
    for ( auto &fl : m_fluids )
      all_cached &= fl.is_cached(frame);
    return all_cached;
  }

  // store given frame for each fluid
  inline void cache(unsigned int frame)
  {
    for ( auto &fl : m_fluids )
      fl.cache(frame);
  }

  inline bool check_and_write_hash() 
  {
    if (global::dynset.savedir.empty())
      return false;

    // open save hash file
    std::fstream fs(global::dynset.hashfile, std::fstream::in );

    std::size_t savedhash = 0;

    if (fs.is_open())
    {
      fs >> std::hex >> savedhash;
      fs.close();
    }

    fs.open(global::dynset.hashfile, std::fstream::out);

    std::size_t curhash = 0;
    boost::hash_combine(curhash, hash_value(global::dynset));
    boost::hash_combine(curhash, hash_value(global::sceneset));
    boost::hash_combine(curhash, hash_value(*this));

    fs.seekp(0);
    fs << std::setfill('0') << std::setw(sizeof(std::size_t) >> 1) << std::hex << curhash;
    fs.close();

    if (savedhash == curhash)
      return true;

    return false;
  }

  template<int PT>
  inline void reset_accel()
  {
    FluidDataVecT<PT> &fldatavec = get_fluiddatas<PT>();
    FluidVec &fluids = get_fluids();
    for ( auto &fldata : fldatavec )
      fluids[fldata.flidx].reset_accel();   
  }

  template<int PT, int PT2, int... PTs>
  inline void reset_accel()
  {
    reset_accel<PT>();
    reset_accel<PT2,PTs...>();
  }


  void run();

  // one leapfrom integrator step
  inline bool step(float dt, bool first_step, float *substep_t=NULL);

  // request stop which will be checked by the owner thread
  void request_stop() { m_stop_requested = true; }
  void unrequest_stop() { m_stop_requested = false; }
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

  friend std::size_t hash_value( const DynamicsManager &dynman )
  {
    std::size_t seed = 0;
    for ( auto &fl : dynman.m_fluids )
      boost::hash_combine(seed, hash_value(fl));
    return seed;
  }

  inline void glprint_fluids(const MaterialManager &matman)
  {
    for ( auto &fl : m_fluids )
    {
      const Vector3f &color = matman[fl.get_material_idx()].kd();
      glprintf_trcv(color, " Fluid ---");
      SPHParticleType pt = fl.get_type();
      glprintf_trcv(color, SPHParticleTypeString[pt]);
      glprintf_trcv(color, "--- \n");
      glprintf_trcv(color, "# of particles: %d  \n", fl.get_num_vertices());
      glprintf_trcv(color, "rest density: %.2f  \n", fl.get_rest_density());
      glprintf_trcv(color, "particle mass: %.2f  \n", fl.get_mass());
      glprintf_trcv(color, "kernel radius: %.2f  \n", fl.get_kernel_radius());
      glprintf_trcv(color, "viscosity: %.2f  \n", fl.get_viscosity());
      glprintf_trcv(color, "surface tension: %.2f  \n", fl.get_surface_tension());
      glprintf_trcv(color, "sound speed: %.2e  \n", std::sqrt(fl.get_sound_speed2()));
    }
  }

  void populate_sph_grid_with_boundaries();
  void populate_sph_grid_with_fluids();

  void populate_sph_grid();

// The following macros generate the fluid vector members and their respective
// getters. 

  template<int PT>
  FluidDataVecT<PT> &get_fluiddatas();
  template<int PT>
  const FluidDataVecT<PT> &get_fluiddatas() const;

  BoundaryPCVec &get_boundaries() { return m_boundaries; }
  FluidVec      &get_fluids() { return m_fluids; }

private: // routine members

  // helper functions to populate_sph_grid
  template <int PT> // base case
  inline void push_fluiddatas_to_sph_grid(Real &color);

  template <int PT, int PT2, int... PTs>
  inline void push_fluiddatas_to_sph_grid(Real &color);

  // helper functions for non templated generate_fluiddatas
  template<int PT>
  inline void generate_fluiddatas()
  {
    // This is one way to do it and is inefficient if there are a lot of fluids.
    // But it's only done at initialization so who cares
    Size num_fluids = m_fluids.size();
    for ( Size i = 0; i < num_fluids; ++i )
    {
      if ( m_fluids[i].get_type() == PT )
      {
        get_fluiddatas<PT>().push_back(FluidDataT<PT>(i, m_fluids));
        get_fluiddatas<PT>().back().init_kernel(m_fluids[i].get_kernel_radius());
      }
    }
  }

  template<int PT, int PT2, int... PTs> // recursive definition
  inline void generate_fluiddatas();

private: // data members
#define FLUIDVEC_MEMBER(z, PT, _) \
  FluidDataVecT<PT> m_fluiddatas_##PT;

  BOOST_PP_REPEAT(NUMFLUIDSPHTYPES, FLUIDVEC_MEMBER, _)

  BoundaryPCVec m_boundaries;

  FluidVec m_fluids; // fluid bases referenced by FluidDataT types

  AlignedBox3f m_bbox; // bounding box of all stored models

  SPHGrid      m_grid; // sph grid

  // dynamics thread concurrency primitives
  std::atomic<bool> m_stop_requested;

  std::mutex m_pause_lock;
  std::condition_variable m_pause_cv;
  std::atomic<bool> m_pause;
};


#define FLUIDVEC_GETTER(z, PT, _) \
  template<> \
  inline FluidDataVecT<PT> &DynamicsManager::get_fluiddatas<PT>() { return m_fluiddatas_##PT; }
#define FLUIDVEC_CONST_GETTER(z, PT, _) \
  template<> \
  inline const FluidDataVecT<PT> &DynamicsManager::get_fluiddatas<PT>() const \
  { return m_fluiddatas_##PT; }

BOOST_PP_REPEAT(NUMFLUIDSPHTYPES, FLUIDVEC_GETTER, _)
BOOST_PP_REPEAT(NUMFLUIDSPHTYPES, FLUIDVEC_CONST_GETTER, _)

#undef FLUIDVEC_GETTER
#undef FLUIDVEC_CONST_GETTER


// Variadic template definitions, uninteresting stuff
template <int PT, int PT2, int... PTs>
inline void
DynamicsManager::push_fluiddatas_to_sph_grid(Real &color)
{
  push_fluiddatas_to_sph_grid<PT>(color);
  push_fluiddatas_to_sph_grid<PT2, PTs...>(color);
}

template<int PT, int PT2, int... PTs>
inline void 
DynamicsManager::generate_fluiddatas()
{
  generate_fluiddatas<PT>();
  generate_fluiddatas<PT2,PTs...>();
}

#endif // DYNAMICSMANAGER_H
