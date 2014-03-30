#include <ctime>
#include <algorithm>
#include <limits>
#include "dynamics.h"
#include "glpointcloud.h"
#include "gltext.h"
#include "fluid.h"
#include "settings.h"

// UniformGrid stuff

template<typename REAL, typename SIZE>
UniformGridRS<REAL,SIZE>::UniformGridRS(
    const Vector3f &bmin,
    const Vector3f &bmax ) 
  : m_num_fluids(0)
  , m_h(2e-3f) // minimum possible cell size
  , m_hinv(500.0f)
  , m_bmin(bmin) // smallest boundary corner
  , m_bmax(bmax) // largest boundary corner
  , m_stop_requested(false)
{ }

template<typename REAL, typename SIZE>
UniformGridRS<REAL,SIZE>::~UniformGridRS() { }
  
template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::init()
{
  FTiter<REAL,SIZE,DEFAULT>::init_fluid_processors(*this); // populate grid with particles
  // set cell size to be the maximum of the kernel support over all fluids;
  for (int j = 0; j < NUMTYPES; ++j)
    for (auto &fl : m_fluids[j])
    {
      if (m_h < fl->get_kernel_radius())
        m_h = fl->get_kernel_radius();
    }

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
  
  FTiter<REAL,SIZE,DEFAULT>::populate_fluid_data(*this); // populate grid with particles
  populate_bound_data();

  int count = 0;
  for (int j = 0; j < NUMTYPES; ++j)
  {
    for (auto &fl : m_fluids[j])
    {
      glprintf_trcv(fl->get_color(), "--- Fluid %d ---\n", ++count);
      glprintf_trcv(fl->get_color(), "# of particles: %d  \n", fl->get_num_vertices());
      glprintf_trcv(fl->get_color(), "rest density: %.2f  \n", fl->get_rest_density());
      glprintf_trcv(fl->get_color(), "particle mass: %.2f  \n", fl->get_mass());
      glprintf_trcv(fl->get_color(), "viscosity: %.2f  \n", fl->get_viscosity());
      glprintf_trcv(fl->get_color(), "surface tension: %.2f  \n", fl->get_surface_tension());
      glprintf_trcv(fl->get_color(), "sound speed: %.2e  \n", std::sqrt(fl->get_sound_speed2()));
    }
  }

  CBVolumeR<REAL> proc;
  proc.init_kernel(m_h);
  m_bound_volume_proc = &proc;
  //compute_bound_quantity< CBVolumeR<REAL> >();
  m_bound_volume_proc = NULL;
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::update_grid()
{
  clear_fluid_data();
  FTiter<REAL,SIZE,DEFAULT>::populate_fluid_data(*this); // populate grid with particles
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::clear_fluid_data()
{
  for (SIZE i = 0; i < m_gridsize[0]; ++i)
    for (SIZE j = 0; j < m_gridsize[1]; ++j)
      for (SIZE k = 0; k < m_gridsize[2]; ++k)
        m_grid[i][j][k].fluidvec.clear();
}

template<typename REAL, typename SIZE, int FT>
void FTiter<REAL,SIZE,FT>::init_fluid_processors(UniformGridRS<REAL,SIZE> &g)
{ 
  if (g.m_fluids[FT].size())
    for (auto &fl : g.m_fluids[FT])
      fl->template cast<FT>()->init_processors();

  FTiter<REAL,SIZE,FT+1>::init_fluid_processors(g);
}

template<typename REAL, typename SIZE, int FT>
void FTiter<REAL,SIZE,FT>::populate_fluid_data(UniformGridRS<REAL,SIZE> &g)
{ 
  if (g.m_fluids[FT].size())
  {
    unsigned short id = 0; // generate ids
    for (auto &fl : g.m_fluids[FT])
    {
      SIZE num_vtx = fl->get_num_vertices();
      for ( SIZE i = 0; i < num_vtx; ++i )
      {
        Vector3R<REAL> pos(fl->get_pos().col(i));
        Vector3R<REAL> vel(fl->get_vel().col(i));
        typename UniformGridRS<REAL,SIZE>::Array3Index idx = g.get_voxel_index(pos);
        g.m_grid(idx).fluidvec.push_back( 
            FluidParticleRT<REAL,FT>( pos, vel, fl->accel_at(i), fl->extern_accel_at(i), id));
      }
      id++;
    }
  }

  FTiter<REAL,SIZE,FT+1>::populate_fluid_data(g);
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::populate_bound_data()
{
  // walk through all boundary cells
  SIZE nx = m_gridsize[0];
  SIZE ny = m_gridsize[1];
  SIZE nz = m_gridsize[2];
  float h = m_h;
  Vector3R<REAL> bmin = m_bmin.template cast<REAL>();

  unsigned char inflate = 3;
  float d = h/inflate;

  SIZE i = 0;
  for (SIZE j = 0; j < ny; ++j)
    for (SIZE k = 0; k < nz; ++k)
    {
      StaticParticles &boundvec = m_grid[i][j][k].boundvec;
      for (unsigned char n = 0; n < inflate; ++n)
        for (unsigned char m = 0; m < inflate; ++m)
          boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(0, h*j+d*n, h*k+d*m) + bmin) );

      if (k == nz-1)
        for (unsigned char n = 0; n < inflate; ++n)
          boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(0, h*j+d*n, h*nz) + bmin) );
      if (j == ny-1)
        for (unsigned char n = 0; n < inflate; ++n)
          boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(0, h*ny, h*k + d*n) + bmin) );
      if (j == ny-1 && k == nz-1)
        boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(0, h*ny, h*nz) + bmin) );
    }

  // cover strip near i=0
  for (SIZE k = 0; k < nz; ++k)
  {
    for (unsigned char n = 1; n < inflate; ++n)
      for (unsigned char m = 0; m < inflate; ++m)
      {
        m_grid[0][0][k].boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(0 + d*n, 0, h*k + d*m) + bmin) );
        m_grid[0][ny-1][k].boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(0 + d*n, h*ny, h*k + d*m) + bmin) );
      }
  }

  for (unsigned char n = 1; n < inflate; ++n)
  {
    m_grid[0][0][nz-1].boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(0 + d*n, 0, h*nz) + bmin) );
    m_grid[0][ny-1][nz-1].boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(0 + d*n, h*ny, h*nz) + bmin) );
  }

  // two blocks near i=0, j=0
  for (unsigned char n = 1; n < inflate; ++n)
    for (unsigned char m = 1; m < inflate; ++m)
    {
      m_grid[0][0][0].boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(0 + d*m, 0+d*n, 0) + bmin) );
      m_grid[0][0][nz-1].boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(0 + d*m, 0+d*n, h*nz) + bmin) );
    }

  for (SIZE j = 1; j < ny; ++j)
  {
    for (unsigned char n = 0; n < inflate; ++n)
      for (unsigned char m = 1; m < inflate; ++m)
      {
        m_grid[i][j][0].boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(0 + d*m, h*j+d*n, 0) + bmin) );
        m_grid[i][j][nz-1].boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(0 + d*m, h*j+d*n, h*nz) + bmin) );
      }
  }
  
  for (i = 1; i < nx; ++i)
  {
    SIZE j = 0;
    for (SIZE k = 0; k < nz; ++k)
    {
      StaticParticles &boundvec = m_grid[i][j][k].boundvec;
      for (unsigned char n = 0; n < inflate; ++n)
        for (unsigned char m = 0; m < inflate; ++m)
          boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(h*i + d*n, 0, h*k + d*m) + bmin) );

      if (k == nz-1)
        for (unsigned char n = 0; n < inflate; ++n)
          boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(h*i + d*n, 0, h*nz) + bmin) );
    }

    // cover strip where j=0
    for (unsigned char n = 0; n < inflate; ++n)
      for (unsigned char m = 1; m < inflate; ++m)
      {
        m_grid[i][j][0].boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(h*i + d*n, 0 + d*m, 0) + bmin) );
        m_grid[i][j][nz-1].boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(h*i + d*n, 0 + d*m, h*nz) + bmin) );
      }

    for (j = 1; j < ny; ++j)
    {
      for (unsigned char n = 0; n < inflate; ++n)
      {
        for (unsigned char m = 0; m < inflate; ++m)
        {
          m_grid[i][j][0].boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(h*i+d*n, h*j+d*m, 0) + bmin) );
          m_grid[i][j][nz-1].boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(h*i+d*n, h*j+d*m, h*nz) + bmin) );
        }
      }
    }

    j = ny-1;
    for (SIZE k = 0; k < nz; ++k)
    {
      StaticParticles &boundvec = m_grid[i][j][k].boundvec;
      for (unsigned char n = 0; n < inflate; ++n)
        for (unsigned char m = 0; m < inflate; ++m)
          boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(h*i+d*n, h*ny, h*k+d*m) + bmin) );
      if (k == nz-1)
        for (unsigned char n = 0; n < inflate; ++n)
          boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(h*i+d*n, h*ny, h*nz) + bmin) );
    }

  } // for i

  i = nx-1; // cap off the right side
  for (SIZE j = 0; j < ny; ++j)
    for (SIZE k = 0; k < nz; ++k)
    {
      StaticParticles &boundvec = m_grid[i][j][k].boundvec;
      for (unsigned char n = 0; n < inflate; ++n)
        for (unsigned char m = 0; m < inflate; ++m)
          boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(h*nx, h*j+d*n, h*k+d*m) + bmin) );
      if (k == nz-1)
        for (unsigned char n = 0; n < inflate; ++n)
          boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(h*nx, h*j+d*n, h*nz) + bmin) );
      if (j == ny-1)
        for (unsigned char n = 0; n < inflate; ++n)
          boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(h*nx, h*ny, h*k+d*n) + bmin) );
      if (j == ny-1 && k == nz-1)
        boundvec.push_back( ParticleR<REAL>(Vector3R<REAL>(h*nx, h*ny, h*nz) + bmin) );
    }

#if 0
  // check to make sure we didn't add the same particle twice
  for (SIZE i = 0; i < nx; ++i)
    for (SIZE j = 0; j < ny; ++j)
      for (SIZE k = 0; k < nz; ++k)
      {
        StaticParticles &bv1 = m_grid[i][j][k].boundvec;
        for (SIZE ni = 0; ni < nx; ++ni)
          for (SIZE nj = 0; nj < ny; ++nj)
            for (SIZE nk = 0; nk < nz; ++nk)
            {
              StaticParticles &bv2 = m_grid[ni][nj][nk].boundvec;
              for (auto &p1 : bv1)
              {
                for (auto &p2 : bv2)
                {
                  if (&p1 == &p2)
                    continue;
                  if (p1.pos[0] == p2.pos[0] && p1.pos[1] == p2.pos[1] && p1.pos[2] == p2.pos[2])
                    qDebug() << "cells " << i << j << k << " and " << ni << nj << nk;

                }
              }
            }
      }
#endif
}

template<typename REAL, typename SIZE>
template<typename ProcessPairFunc, typename ParticleType>
void UniformGridRS<REAL,SIZE>::compute_quantity()
{
  SIZE nx = m_gridsize[0];
  SIZE ny = m_gridsize[1];
  SIZE nz = m_gridsize[2];
  for (SIZE i = 0; i < nx; ++i)
  {
    for (SIZE j = 0; j < ny; ++j)
    {
      for (SIZE k = 0; k < nz; ++k)
      {
        auto &pvec = 
          m_grid[i][j][k].template get_vec<ParticleType>();
        if (pvec.empty())
          continue;

        IndexRange xrange = range3(i,nx);
        IndexRange yrange = range3(j,ny);
        IndexRange zrange = range3(k,nz);
        GridView neigh_view = m_grid[ boost::indices[xrange][yrange][zrange] ];

        for ( auto &p : pvec )  // prepare data
        {
          ProcessPairFunc &process = determine_proc< ProcessPairFunc, ParticleType >(p);
          process.init_particle(p);
        }

        SIZE xrange_size = xrange.finish() - xrange.start();
        SIZE yrange_size = yrange.finish() - yrange.start();
        SIZE zrange_size = zrange.finish() - zrange.start();
        for (SIZE near_i = 0; near_i < xrange_size; ++near_i)
        {
          for (SIZE near_j = 0; near_j < yrange_size; ++near_j)
          {
            for (SIZE near_k = 0; near_k < zrange_size; ++near_k)
            {
              Cell &cell = neigh_view[near_i][near_j][near_k];
              DynamicParticles &neigh_fluidvec = cell.fluidvec;
              StaticParticles  &neigh_boundvec = cell.boundvec;
              for ( auto &p : pvec )
              {
                ProcessPairFunc &process = determine_proc< ProcessPairFunc, ParticleType>(p);
                for ( FluidParticleR<REAL> &near_p : neigh_fluidvec )
                {
                  process.fluid(p, near_p); // process neighbouring fluid data
                }
                for ( ParticleR<REAL> &near_p : neigh_boundvec )
                {
                  process.bound(p, near_p); // process neighbouring boundary data
                }
              }
            }
          }
        }

        for ( auto &p : pvec )  // finalize data
        {
          ProcessPairFunc &process = determine_proc< ProcessPairFunc, ParticleType >(p);
          process.finish_particle(p);
        }
      } // for k
    } // for j
  } // for i
}

template<typename REAL, typename SIZE, int FT>
void FTiter<REAL,SIZE,FT>::update_density(UniformGridRS<REAL,SIZE> &g,float timestep)
{
  if (g.m_fluids[FT].size())
  {
    for (auto &fl : g.m_fluids[FT])
      fl->template cast<FT>()->m_fluid_density_update_proc.init( timestep ); 

    g.template compute_density_updateT<FT>();
  }
  FTiter<REAL,SIZE,FT+1>::update_density(g,timestep); // recurse
}

template<typename REAL, typename SIZE, int FT>
void FTiter<REAL,SIZE,FT>::compute_density(UniformGridRS<REAL,SIZE> &g)
{
  if (g.m_fluids[FT].size())
  {

    REAL max_var[g.m_num_fluids];
    REAL avg_var[g.m_num_fluids];
    
    int i = 0;
    for (auto &fl : g.m_fluids[FT] )
    {
      max_var[i] = avg_var[i] = 0.0f;
      fl->template cast<FT>()->m_fluid_density_proc.init( max_var[i], avg_var[i] );
      i++;
    }

    g.template compute_densityT<FT>();

#ifdef REPORT_DENSITY_VARIATION
    i = 0;
    for (auto &fl : g.m_fluids[FT] )
    {
      avg_var[i] = avg_var[i]/fl->cast<FT>()->get_num_vertices();
      qDebug("Fluid %d:  max: %.0f, %.1f percent;    avg: %.0f, %.1f percent", i,
          max_var[i], 100.00f*max_var[i]/fl->cast<FT>()->get_rest_density(), 
          avg_var[i], 100.00f*avg_var[i]/fl->cast<FT>()->get_rest_density());
      i++;
    }
#endif
  }
  FTiter<REAL,SIZE,FT+1>::compute_density(g); // recurse
}

template<typename REAL, typename SIZE, int FT>
void FTiter<REAL,SIZE,FT>::compute_pressure(UniformGridRS<REAL,SIZE> &g)
{
  if (g.m_fluids[FT].size())
    g.template compute_pressureT<FT>();
  FTiter<REAL,SIZE,FT+1>::compute_pressure(g); // recurse
}

template<typename REAL, typename SIZE, int FT>
void FTiter<REAL,SIZE,FT>::compute_accel(UniformGridRS<REAL,SIZE> &g)
{
  if (g.m_fluids[FT].size())
  {
    for (auto &fl : g.m_fluids[FT])
      fl->template cast<FT>()->reset_accel();     // now may assume all accelerations are zero

    g.template compute_pressure_accelT<FT>();
    g.template compute_viscosity_accelT<FT>();
  }

  FTiter<REAL,SIZE,FT+1>::compute_accel(g); // recurse
}
template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::run()
{
  if (m_num_fluids < 1)
    return;

  init();

  // timestep
  float dt = 1.0f/(global::dynset.fps * global::dynset.substeps);

  glprintf_tr("step: %.2es\n", dt);
  float rdt = dt;

  bool cached = true; // true if frame is cached

  // check if first frame is cached
  for (int j = 0; j < NUMTYPES; ++j)
    for (auto &fl : m_fluids[j])
      cached &= fl->is_cached(0); // check if we have already cached the first frame

  float file_read_t = 0.0f;
  float frame_t = 0.0f;
  float substep_t = 0.0f;
  unsigned int file_reads = 0;
  unsigned int frame_count = 0;
  unsigned int frame = 1;

  for ( ; ; ++frame, ++frame_count ) // for each frame
  {
    // check if we have the next frame cached for ALL of the fluids
    clock_t s = clock();
    cached = true;
    for (int j = 0; j < NUMTYPES; ++j)
      for (auto &fl : m_fluids[j])
        cached &= fl->is_cached(frame);

    if (cached) // if frame is cached just load it up
    {
      for (int j = 0; j < NUMTYPES; ++j)
        for (auto &fl : m_fluids[j])
          fl->read_cache(frame);

      file_read_t += float(clock() - s) / CLOCKS_PER_SEC;
      file_reads += 1;

      for (int j = 0; j < NUMTYPES; ++j)
        for (auto &fl : m_fluids[j])
          fl->update_data(); // notify gl we have new positions

    }
    else // compute next frame
    {
      substep_t = 0.0f;
      for (unsigned int iter = 0; iter < global::dynset.substeps; ++iter)
      { // for each simulation substep
#if 0
        ///////// testing adaptive time step
        float fdt = std::numeric_limits<float>::infinity();
        float cvdt = std::numeric_limits<float>::infinity();
        for (int j = 0; j < NUMTYPES; ++j)
          for (auto &fl : m_fluids[j])
          {
            for (SIZE i = 0; i < fl->get_num_vertices(); ++i)
            {
              fdt = std::min(fdt,
                  float(fl->get_radius() /
                    (fl->get_mass()*Vector3R<REAL>(fl->get_extern_accel().col(i)).norm())));
              cvdt = std::min(cvdt,
                  float(fl->get_radius()/(fl->get_sound_speed2()*(1+0.6*fl->get_viscosity()))));
              // check for nan and infs
              //if (!std::isfinite(fl->vel_at(i)[0]) ||
              //    !std::isfinite(fl->vel_at(i)[1]) || 
              //    !std::isfinite(fl->vel_at(i)[2]))
              //  qDebug() << "Found NaN at " << i;
            }
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

        clock_t prev_t = clock();

        //if (count % 100)
        //  update_density(dt);
        //else
          FTiter<REAL,SIZE,DEFAULT>::compute_density(*this);
        FTiter<REAL,SIZE,DEFAULT>::compute_pressure(*this);
        FTiter<REAL,SIZE,DEFAULT>::compute_accel(*this); // update m_accel

        float factor = 1.0f; // leap-frog method has a different first step
        if (iter == 0 && frame == 1)
          factor = 0.5f;

        for (int j = 0; j < NUMTYPES; ++j)
        {
          for (auto &fl : m_fluids[j])
          {
            fl->get_vel() = fl->get_vel() + factor*dt*fl->get_accel();
            fl->get_pos() = (fl->get_pos() + dt*fl->get_vel()).eval();
            if (j == MCG03)
              fl->resolve_collisions();

            // prepare velocities for acceleration computation
            fl->get_vel() = fl->get_vel() + 0.5*dt*fl->get_accel();
          }
        }

        update_grid(); // prepare grid for next simulation step

        substep_t += float(clock() - prev_t) / CLOCKS_PER_SEC;

        if (m_stop_requested) break;

        for (int j = 0; j < NUMTYPES; ++j)
          for (auto &fl : m_fluids[j])
            fl->update_data(); // notify gl we have new positions

      } // for each substep

      for (int j = 0; j < NUMTYPES; ++j)
        for (auto &fl : m_fluids[j])
          fl->write_cache(frame); // save frame to cache

    } // if not cached

    frame_t += float(clock() - s) / CLOCKS_PER_SEC;

    glprintf_tl("\ravg time per:  substep = %.2es,   frame = %.2es,   read = %.2es",
        substep_t / float(global::dynset.substeps),
        frame_t / float(frame_count),
        file_read_t / float(file_reads));

    if (m_stop_requested) break;

    if (frame >= global::dynset.frames)
    {
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

template class UniformGridRS<double, unsigned int>;
