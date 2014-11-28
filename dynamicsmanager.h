#ifndef DYNAMICSMANAGER_H
#define DYNAMICSMANAGER_H

#include <ctime>
#include <algorithm>
#include <limits>
#include <fstream>
#include <iomanip>

#include "types.h"
#include "gltext.h"
#include "boundary.h"
#include "fluid.h"
#include "fluidimpl.h"
#include "settings.h"
#define BOOST_CHRONO_HEADER_ONLY
#include <boost/chrono.hpp>

typedef boost::chrono::process_real_cpu_clock real_clock;
typedef boost::chrono::nanoseconds ns_t;
typedef real_clock::time_point real_t;

// Helper template using two arguments from the data arguemnt
// where data is of the form (fluid, for loop body)

#define FLUID_TYPE_CASE_FUNC_TEMPLATE(z, PT, data) \
  case PT: \
    BOOST_PP_TUPLE_ELEM(2, 0, data)##BOOST_PP_TUPLE_ELEM(2, 1, data); \
    break;

#define FLUID_TYPE_CASE_FUNC_CALL(func, params) \
  BOOST_PP_REPEAT(NUMFLUIDTYPES, FLUID_TYPE_CASE_FUNC_TEMPLATE, (func, params))

class DynamicsManager
{
public:
  DynamicsManager() 
    : m_stop_requested(false)
    , m_pause(false)
  { }
  ~DynamicsManager() { }

  // templated "typedef"
  template<int PT>
  using FluidDataVecT = std::vector< FluidDataT<PT> >;

  /// Manager interface

  // normalization routines used to fit the dynamic objects within a unit cube
  void normalize_models()
  {
    normalize_model(Vector3f(0.0f, 0.0f, 0.0f),Vector3f(0.0f, 0.0f, 0.0f));
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
    transform_models(Scaling(2.0f/sizevec.maxCoeff()));
    Vector3f box_center = m_bbox.center();
    transform_models(Translation(-box_center));
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
                          DynParamsPtr dynparams)
  {
    if ( mesh->mPrimitiveTypes & aiPrimitiveType_POINT )
    {
      if ( dynparams->type == FLUID )
      {
        FluidParamsPtr fluidparams =
          boost::static_pointer_cast<FluidParams>(dynparams);
        switch (fluidparams->fluid_type)
        {
          FLUID_TYPE_CASE_FUNC_CALL(add_fluid, (mesh, mat_idx, fluidparams))
          default: break;
        }
      }
    }
    // else ignore. TODO: extend this to include other dynamic objects
  }

  void add_static_object(const aiMesh *mesh, Index mat_idx)
  {
    // implement this for arbitrary boundaries
  }

  template<int PT>
  void add_fluid(const aiMesh *mesh, Index mat_idx, FluidParamsPtr fparams)
  {
    m_fluids.push_back(Fluid(mesh, mat_idx, fparams));
    get_fluiddatas<PT>().push_back(FluidDataT<PT>(m_fluids.back()));
  }

  void init_fluids(SPHGrid &simulation_grid)
  {
    for ( auto &fl : m_fluids )
      fl.init(simulation_grid);
  }

  inline Size get_num_fluids() { return m_fluids.size(); }
  inline Size get_num_boundaries() { return m_boundaries.size(); }

  inline void clear()
  {
    m_fluids.clear();
    m_boundaries.clear();
    clear_fluiddatas<ALL_FLUID_PARTICLE_TYPES>();
  }

  template <int PT>
  inline void clear_fluiddatas()
  {
    get_fluiddatas<PT>().clear();
  }

  template <int PT, int PT2, int... PTs>
  inline void clear_fluiddatas()
  {
    clear_fluiddata<PT>();
    clear_fluiddata<PT2, PTs>();
  }

  inline Real get_max_radius()
  {
    Real h = 0.0;
    for ( auto &fl : m_fluids )
    {
      if (h < fl->get_kernel_radius())
        h = fl->get_kernel_radius();
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

  void run(SPHGrid &grid)
  {
    if (get_num_fluids() < 1)
      return;

    // timestep
    float dt = 1.0f/(global::dynset.fps * global::dynset.substeps);

    glprintf_tr("step: %.2es\n", dt);
    float rdt = dt;

    bool all_cached = check_and_write_hash();

    if (all_cached) // try to load cached frames into the fluid objects
      all_cached &= load_saved_cache();
    
    if (!all_cached)
    {
      // Initialize the fluid for init_steps steps before simulating
      // temporarily disable gravity
      Vector3f grav = global::dynset.gravity;
      global::dynset.gravity = Vector3f(0.0,0.0,0.0);
      for (unsigned int iter = 0; iter < global::dynset.init_steps; ++iter)
      { // for each simulation substep
        if (!step(dt, iter == 0))
          break;
      } // for each substep

      global::dynset.gravity = grav; // restore gravity

      if (m_stop_requested)
        return;

      cache(0);
    }

    real_t start_time;
    float file_read_t = 0.0f;
    float frame_t = 0.0f;
    float substep_t = 0.0f;
    unsigned int file_reads = 0;
    unsigned int frame_count = 0;
    unsigned int frame = 1;

    for ( ; ; ++frame, ++frame_count ) // for each frame
    {
      start_time = real_clock::now();
      clock_t s = clock();

      // check if we have the next frame cached for ALL of the fluids
      cached = is_cached(frame);

      if (cached)
      {
        load_cached(frame); // load cached frame
        prepare_vis_data(); // notify gl we have new positions

        file_reads += 1;
        file_read_t += float(clock() - s) / CLOCKS_PER_SEC;

        ns_t elapsed = (real_clock::now() - start_time);
        ns_t perframe(boost::int_least64_t(1.0e9/global::dynset.fps));
        if (elapsed < perframe)
          std::this_thread::sleep_for(
              std::chrono::nanoseconds((perframe - elapsed).count()));
      }
      else // if not all cached, compute substeps
      {
        substep_t = 0.0f;
        for (unsigned int iter = 0; iter < global::dynset.substeps; ++iter)
        { // for each simulation substep
          bool first_step = iter == 0 && frame == 1 && global::dynset.init_steps == 0;
          if (!step(dt, first_step, &substep_t))
            break;
        } // for each substep

        cache(frame);
      }

      frame_t += float(clock() - s) / CLOCKS_PER_SEC;

      glprintf_tl("\ravg time per:  substep = %.2es,   frame = %.2es,   read = %.2es",
          substep_t / float(global::dynset.substeps),
          frame_t / float(frame_count),
          file_read_t / float(file_reads));

      {
        std::unique_lock<std::mutex> locker(m_pause_lock);

        while (m_pause) // avoid spurious wakeup
          m_pause_cv.wait(locker);
      }

      if (m_stop_requested) break;

      if (frame >= global::dynset.frames)
      { // restart simulation
        // for timing
        file_reads = 0;
        frame_count = 0;
        frame_t = 0.0f;
        file_read_t = 0.0f;

        // actual frame
        frame = 0;
      }
    }

  }

  // one leapfrom integrator step
  inline bool step(float dt, bool first_step, SPHGrid &grid, float *substep_t = NULL)
  {
#if 0
    ///////// testing adaptive time step
    float fdt = std::numeric_limits<float>::infinity();
    float cvdt = std::numeric_limits<float>::infinity();
    for (int j = 0; j < NUMFLUIDTYPES; ++j)
      for (auto &fl : m_fluids[j])
      {
        for (Size i = 0; i < fl->get_num_vertices(); ++i)
        {
          fdt = std::min(fdt,
              float(fl->get_radius() /
                (fl->get_mass()*Vector3R<Real>(fl->get_extern_accel().col(i)).norm())));
          // check for nan and infs
          //if (!std::isfinite(fl->vel_at(i)[0]) ||
          //    !std::isfinite(fl->vel_at(i)[1]) || 
          //    !std::isfinite(fl->vel_at(i)[2]))
          //  qDebug() << "Found NaN at " << i;
        }
        cvdt = std::min(cvdt,
            float(fl->get_radius()/(fl->get_sound_speed2()*(1+0.6*fl->get_viscosity()))));
      }

    float new_rdt = std::min(0.25*fdt, 0.4*cvdt);
    if (new_rdt != rdt)
    {
      glprintf_tr("\rrecommended step: %.2es", new_rdt);
      rdt = new_rdt;
      //dt = rdt;
    }
    //////// end of adaptive timestep test
#endif

    grid.update_grid(); // prepare grid for simulation step

    clock_t prev_t = clock();

    //if (iter)
    // PTiter<DEFAULT>::update_density(*this, dt);
    //else
    grid.compute_density();


    //PTiter<DEFAULT>::compute_accel(*this); // update m_accel

    //if (m_fluids[MCG03].size()) // extra steps to compute surface tension for MCG03
    //{
    //  compute_surface_normalT<MCG03>();
    //  compute_surface_tensionT<MCG03>();
    //}

    if (substep_t)
      *substep_t += float(clock() - prev_t) / CLOCKS_PER_SEC;

    float factor = first_step ? 0.5f : 1.0f; // leap-frog method has a different first step

    for ( auto &fl : m_fluids )
      fl->get_vel() = fl->get_vel() + factor*dt*fl->get_accel();

    //jacobi_pressure_solve(dt,factor); // only relevant for IISPH fluids

    for ( auto &fl : m_fluids )
    {
      fl->get_pos() = (fl->get_pos() + dt*fl->get_vel()).eval();
      if (fl.get_type() == MCG03)
        fl->resolve_collisions();
      else
        fl->clamp(0.5*grid.get_cell_size()-0.01,0.01f);

      // prepare velocities for acceleration computation in next step
      fl->get_vel() = fl->get_vel() + 0.5*dt*fl->get_accel();
    }

    if (m_stop_requested) return false;

    // TODO: put this step after each frame instead of after each substep
    for ( auto &fl : m_fluids )
      fl->prepare_vispos(); // copy positions for visualization

    return true;
  }

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

  friend std::size_t hash_value( const DynamicsManager &dynman )
  {
    std::size_t seed = 0;
    for ( auto &fl : dynman.m_fluids )
      boost::hash_combine(seed, hash_value(*fl));
    return seed;
  }

  inline void glprint_fluids(const MaterialManager &matman)
  {
    for ( auto &fl : m_fluids )
    {
      const Vector3f &color = matman.get_material(fl.get_material_idx()).kd();
      glprintf_trcv(, " Fluid ---");
      ParticleType pt = fl.get_type();
      switch (pt)
      {
        case 1: glprintf_trcv(color, "MCG03"); break;
        case 2: glprintf_trcv(color, "BT07"); break;
        case 3: glprintf_trcv(color, "ICS13"); break;
        default: break;
      }
      glprintf_trcv(color, "--- \n");
      glprintf_trcv(color, "# of particles: %d  \n", fl.get_num_vertices());
      glprintf_trcv(color, "rest density: %.2f  \n", fl.get_rest_density());
      glprintf_trcv(color, "particle mass: %.2f  \n", fl.get_mass());
      glprintf_trcv(color, "particle radius: %.2f  \n", fl.get_radius());
      glprintf_trcv(color, "viscosity: %.2f  \n", fl.get_viscosity());
      glprintf_trcv(color, "surface tension: %.2f  \n", fl.get_surface_tension());
      glprintf_trcv(color, "sound speed: %.2e  \n", std::sqrt(fl.get_sound_speed2()));
    }
  }

  inline void populate_sph_grid_with_boundaries(SPHGrid &grid)
  {
    for ( auto &bound : m_boundaries )
    {
      Size num_vtx = bound.get_num_vertices();
      for ( Size i = 0; i < num_vtx; ++i )
      {
        Vector3R<Real> pos(bound.get_pos().col(i));
        typename SPHGrid::Array3Index idx = grid.get_voxel_index(pos);
        grid.get_cell_at(idx).push_particle(bound, i);
      }
    }
  }

  inline void populate_sph_grid_with_fluids(SPHGrid &grid)
  {
    Real color = 1.0;
    push_fluiddatas_to_sph_grid<ALL_FLUID_PARTICLE_TYPES>(grid, color);
  }

  inline void populate_sph_grid(SPHGrid &grid)
  { 
    // fluids
    populate_sph_grid_with_fluids(grid);

    // statics
    populate_sph_grid_with_boundaries(grid);
  }

// The following macros generate the fluid vector members and their respective
// getters. 
#define FLUIDVEC_GETTER(z, PT, _) \
  inline FluidDataVecT<PT> &get_fluiddatas<PT>() { return m_fluiddatas_##PT; }
#define FLUIDVEC_CONST_GETTER(z, PT, _) \
  inline const FluidDataVecT<PT> &get_fluiddatas<PT>() const \
  { return m_fluiddatas_##PT; }

  BOOST_PP_REPEAT(NUMFLUIDTYPES, FLUIDVEC_GETTER, _)
  BOOST_PP_REPEAT(NUMFLUIDTYPES, FLUIDVEC_CONST_GETTER, _)

#undef FLUIDVEC_GETTER
#undef FLUIDVEC_CONST_GETTER

  BoundaryPCVec &get_bounds() { return m_boundaries; }
  FluidVec      &get_fluids() { return m_fluids; }

private: // routine members

  // helper functions to populate_sph_grid
  template <int PT> // base case
  inline void push_fluiddatas_to_sph_grid(SPHGrid &grid, Real &color)
  {
    FluidDataVecT<PT> fldatavec = get_fluiddatas<PT>();
    for ( auto &fldata : fldatavec )
    {
      Size num_vtx = fldata.fl.get_num_vertices();
      for ( Size i = 0; i < num_vtx; ++i )
      {
        Vector3R<Real> pos(fldata.fl.get_pos().col(i));
        typename SPHGrid::Array3Index idx = grid.get_voxel_index(pos);
        grid.get_cell_at(idx).push_particle<PT>(fldata, color, i);
      }
      color += 1.0f;
    }
  }

  template <int PT, int PT2, int... PTs>
  inline void push_fluiddatas_to_sph_grid(SPHGrid &grid, Real &color)
  {
    push_fluiddatas_to_sph_grid<PT>(grid, color);
    push_fluiddatas_to_sph_grid<PT2, PTs...>(grid, color);
  }

private: // data members
#define FLUIDVEC_MEMEBER(z, PT, _) \
  FluidDataVecT<PT> m_fluiddatas_##PT;

  BOOST_PP_REPEAT(NUMFLUIDTYPES, FLUIDVEC_MEMBER, _)

  BoundaryPCVec m_boundaries;

  FluidVec m_fluids; // fluid bases referenced by FluidDataT types

  AlignedBox3f m_bbox; // bounding box of all stored models

  // dynamics thread concurrency primitives
  std::atomic<bool> m_stop_requested;

  std::mutex m_pause_lock;
  std::condition_variable m_pause_cv;
  std::atomic<bool> m_pause;
};

#endif // DYNAMICSMANAGER_H
