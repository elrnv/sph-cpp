#ifndef DYNAMICSMANAGER_H
#define DYNAMICSMANAGER_H

#include "types.h"
#include "fluidimpl.h"

// Helper template using two arguments from the data arguemnt
// where data is of the form (fluid, for loop body)

# if 1
#define FOREACH_FLUIDDATA_TEMPLATE(z, PT, data) \
  for ( BOOST_PP_TUPLE_ELEM(2, 0, data) : m_fluiddatas_##PT ) \
  { BOOST_PP_TUPLE_ELEM(2, 1, data); }

#define FOREACH_FLUIDDATA(fl, statement) \
  do { \
    BOOST_PP_REPEAT(NUMFLUIDTYPES, FOREACH_FLUIDDATA_TEMPLATE, (fl, statement)) \
  } while (false)


#define FLUID_TYPE_CASE_FUNC_TEMPLATE(z, PT, data) \
  case PT: \
    BOOST_PP_TUPLE_ELEM(2, 0, data)##BOOST_PP_TUPLE_ELEM(2, 1, data); \
    break;

#define FLUID_TYPE_CASE_FUNC_CALL(func, params) \
  BOOST_PP_REPEAT(NUMFLUIDTYPES, FLUID_TYPE_CASE_FUNC_TEMPLATE, (func, params))

#endif

class DynamicsManager
{
public:
  DynamicsManager() { }
  ~DynamicsManager() { }

  // templated "typedef"
  template<int PT>
  using FluidDataVecT = std::vector< FluidDataT<PT> >;

  /// Manager interface
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

  inline Size get_numfluids() { return m_fluids.size(); }

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

  inline void glprint_fluids()
  {
    for ( auto &fl : m_fluids )
    {
      glprintf_trcv(fl->get_color(), " Fluid ---");
      ParticleType pt = fl.get_type();
      switch (pt)
      {
        case 1: glprintf_trcv(fl->get_color(), "MCG03"); break;
        case 2: glprintf_trcv(fl->get_color(), "BT07"); break;
        case 3: glprintf_trcv(fl->get_color(), "ICS13"); break;
        default: break;
      }
      glprintf_trcv(fl->get_color(), "--- \n");
      glprintf_trcv(fl->get_color(), "# of particles: %d  \n", fl->get_num_vertices());
      glprintf_trcv(fl->get_color(), "rest density: %.2f  \n", fl->get_rest_density());
      glprintf_trcv(fl->get_color(), "particle mass: %.2f  \n", fl->get_mass());
      glprintf_trcv(fl->get_color(), "particle radius: %.2f  \n", fl->get_radius());
      glprintf_trcv(fl->get_color(), "viscosity: %.2f  \n", fl->get_viscosity());
      glprintf_trcv(fl->get_color(), "surface tension: %.2f  \n", fl->get_surface_tension());
      glprintf_trcv(fl->get_color(), "sound speed: %.2e  \n", std::sqrt(fl->get_sound_speed2()));
    }
  }

  inline void populate_sph_grid(SPHGrid &grid)
  { 
    // fluids
    Real color = 1.0f;
    FOREACH_FLUIDDATA(auto &fldata, push_fluiddata_to_sph_grid<PT>(fldata,grid,color));

    // statics
    push_bounddatas_to_sph_grid(grid);
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
  FluidVec &get_fluids() { return m_fluids; }

private: // routine members
  // helper function to populate_sph_grid
  template <int PT>
  inline void push_fluiddata_to_sph_grid<PT>(
      FluidDataT<PT> &fldata, SPHGrid &grid, Real &color)
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

  inline void push_bounddatas_to_sph_grid<PT>(SPHGrid &grid)
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

private: // data members
#define FLUIDVEC_MEMEBER(z, PT, _) \
  FluidDataVecT<PT> m_fluiddatas_##PT;

  BOOST_PP_REPEAT(NUMFLUIDTYPES, FLUIDVEC_MEMBER, _)

  BoundaryPCVec m_boundaries;

  FluidVec m_fluids; // fluid bases referenced by FluidDataT types
};

#endif // DYNAMICSMANAGER_H
