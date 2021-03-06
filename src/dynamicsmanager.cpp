#include <thread>
#include "sphgrid.h"
#include "dynamicsmanager.h"

void 
DynamicsManager::run()
{
  if (get_num_fluids() < 1)
    return;

  // timestep
  float dt = 1.0f/(global::dynset.fps * global::dynset.substeps);

  glprintf_tr("step: %.2es\n", dt);
  //float rdt = dt;

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

  prepare_vis_data();

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
    bool cached = is_cached(frame);

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

    prepare_vis_data();

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

inline bool
DynamicsManager::step(float dt, bool first_step, float *substep_t)
{
#if 0
  ///////// testing adaptive time step
  float fdt = std::numeric_limits<float>::infinity();
  float cvdt = std::numeric_limits<float>::infinity();
  for (int j = 0; j < NUMFLUIDSPHTYPES; ++j)
    for (auto &fl : m_fluids[j])
    {
      for (Size i = 0; i < fl.get_num_vertices(); ++i)
      {
        fdt = std::min(fdt,
            float(fl.get_radius() /
              (fl->get_mass()*Vector3T<Real>(fl->get_extern_accel().col(i)).norm())));
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

  clear_fluid_sph_particles(); // prepare grid for simulation step
  populate_sph_grid_with_fluids();
  //m_grid.update_particle_cell_positions();

  clock_t prev_t = clock();

  //if (iter)
  // PTiter<DEFAULT>::update_density(*this, dt);
  //else
  //m_grid.compute_quantity<Density, ALL_FLUID_SPH_TYPES>();
  m_grid.compute_density<ALL_FLUID_SPH_TYPES>();
  reset_accel<ALL_FLUID_SPH_TYPES>();
  m_grid.compute_quantity<Accel,ALL_FLUID_SPH_TYPES>(); // update m_accel

  //if (m_fluids[MCG03].size()) // extra steps to compute surface tension for MCG03
  //{
  //  compute_surface_normalT<MCG03>();
  //  compute_surface_tensionT<MCG03>();
  //}

  if (substep_t)
    *substep_t += float(clock() - prev_t) / CLOCKS_PER_SEC;

  float factor = first_step ? 0.5f : 1.0f; // leap-frog method has a different first step

  for ( auto &fl : m_fluids )
    fl.get_vel() = fl.get_vel() + factor*dt*fl.get_accel();

  //jacobi_pressure_solve(dt,factor); // only relevant for IISPH fluids

  for ( auto &fl : m_fluids )
  {
    fl.get_pos() = (fl.get_pos() + dt*fl.get_vel()).eval();
    if (fl.get_type() == MCG03)
      fl.resolve_collisions();
    //else
    //  fl.clamp(0.5*m_grid.get_cell_size()-0.01,0.01f);

    // prepare velocities for acceleration computation in next step
    fl.get_vel() = fl.get_vel() + 0.5*dt*fl.get_accel();
  }

  if (m_stop_requested) return false;

  prepare_vis_data();

  return true;
}

void 
DynamicsManager::populate_sph_grid_with_boundaries()
{
  for ( auto &bound : m_boundaries )
  {
    Size num_vtx = bound.get_num_vertices();
    for ( Size i = 0; i < num_vtx; ++i )
    {
      Vector3T<Real> pos(bound.get_pos().col(i));
      typename SPHGrid::Array3Index idx = m_grid.get_voxel_index(pos);
      m_grid.get_cell_at(idx).push_particle(bound, i);
    }
  }

  // compute the needed boundary data right away
  m_grid.compute_quantity<Volume,STATIC>();
}

void 
DynamicsManager::populate_sph_grid_with_fluids()
{
  Real color = 1.0;
  push_fluiddatas_to_sph_grid<ALL_FLUID_SPH_TYPES>(color);
}

void 
DynamicsManager::populate_sph_grid()
{ 
  // fluids
  populate_sph_grid_with_fluids();

  // statics
  populate_sph_grid_with_boundaries();
}

template <int PT> // base case
inline void
DynamicsManager::push_fluiddatas_to_sph_grid(Real &color)
{
  FluidDataVecT<PT> &fldatavec = get_fluiddatas<PT>();
  for ( auto &fldata : fldatavec )
  {
    Size num_vtx = m_fluids[fldata.flidx].get_num_vertices();
    for ( Size i = 0; i < num_vtx; ++i )
    {
      Vector3T<Real> pos(m_fluids[fldata.flidx].get_pos().col(i));
      typename SPHGrid::Array3Index idx = m_grid.get_voxel_index(pos);
      m_grid.get_cell_at(idx).push_particle<PT>(fldata, m_fluids, color, i);
    }
    color += 1.0f;
  }
}

// instance the function above for each fluid type
#define INSTANCE_PUSH_FLUIDDATAS_TEMPLATE(z, PT, _) \
  template void DynamicsManager::push_fluiddatas_to_sph_grid<PT>(Real &);
BOOST_PP_REPEAT(NUMFLUIDSPHTYPES, INSTANCE_PUSH_FLUIDDATAS_TEMPLATE, _)
