#include <ctime>
#include <algorithm>
#include <limits>
#include <fstream>
#include <iomanip>
#include "dynamics.h"
#include "glpointcloud.h"
#include "gltext.h"
#include "fluid.h"
#include "settings.h"
#include "quantityprocessor.h"
#define BOOST_CHRONO_HEADER_ONLY
#include <boost/chrono.hpp>

typedef boost::chrono::process_real_cpu_clock real_clock;
typedef real_clock::time_point real_t;

// From [Solenthaler and Pajarola 2008], an alternative density
// (number_density * mass) is used
template<typename REAL, typename SIZE>
class CBVolumeRS : public CBQ<REAL, SIZE, CBVolumeRS<REAL,SIZE> >
{
public:
  inline void init_kernel(float h) { m_kern.init(h); }
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
private:
  CubicSplineKernel m_kern;
};

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
  , m_pause(false)
{ }

template<typename REAL, typename SIZE>
UniformGridRS<REAL,SIZE>::~UniformGridRS() { }
  
template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::init()
{
  if (m_num_fluids < 1)
    return;

  FTiter<REAL,SIZE,DEFAULT>::init_fluid_processors(*this);

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
  
  //populate_fluid_data(); // populate grid with particles
  populate_bound_data();

  for (int j = 0; j < NUMTYPES; ++j)
  {
    for (auto &fl : m_fluids[j])
    {
      glprintf_trcv(fl->get_color(), " Fluid ---");
      switch (j)
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
  CBVolumeRS<REAL,SIZE> proc;
  proc.init_kernel(m_h);
  m_bound_volume_proc = &proc;
  compute_bound_quantity< CBVolumeRS<REAL,SIZE> >();
  m_bound_volume_proc = NULL;
}

template<typename REAL, typename SIZE>
inline void UniformGridRS<REAL,SIZE>::update_grid()
{
  clear_fluid_data();
  populate_fluid_data(); // populate grid with particles
}

template<typename REAL, typename SIZE>
inline void UniformGridRS<REAL,SIZE>::clear_fluid_data()
{
  for (SIZE i = 0; i < m_gridsize[0]; ++i)
    for (SIZE j = 0; j < m_gridsize[1]; ++j)
      for (SIZE k = 0; k < m_gridsize[2]; ++k)
        m_grid[i][j][k].fluidvec.clear();
}

template<typename REAL, typename SIZE, int FT>
inline void FTiter<REAL,SIZE,FT>::init_fluid_processors(UniformGridRS<REAL,SIZE> &g)
{ 
  if (g.m_fluids[FT].size())
    for (auto &fl : g.m_fluids[FT])
      fl->template cast<FT>()->init_processors();

  FTiter<REAL,SIZE,FT+1>::init_fluid_processors(g);
}

template<typename REAL, typename SIZE>
inline void UniformGridRS<REAL,SIZE>::populate_fluid_data()
{ 
  REAL color = 1.0f;
  for (int ft = DEFAULT; ft < NUMTYPES; ++ft)
  {
    unsigned short id = 0; // generate ids
    for (auto &fl : m_fluids[ft])
    {
      SIZE num_vtx = fl->get_num_vertices();
      for ( SIZE i = 0; i < num_vtx; ++i )
      {
        Vector3R<REAL> pos(fl->get_pos().col(i));
        Vector3R<REAL> vel(fl->get_vel().col(i));
        typename UniformGridRS<REAL,SIZE>::Array3Index idx = get_voxel_index(pos);
        m_grid(idx).fluidvec.push_back(
            FluidParticleR<REAL>( pos, vel, 
              fl->accel_at(i), 
              fl->extern_accel_at(i), 
              fl->dinv_at(i), 
              id, color));
      }
      id++;
      color += 1.0f;
    }
  }
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::populate_bound_data()
{
  unsigned char particles_per_cell_length = 4;

  // place boundary particles slightly away from the actual boundary
  float pad = m_h;

  SIZE nx = particles_per_cell_length*(m_gridsize[0]+2);
  SIZE ny = particles_per_cell_length*(m_gridsize[1]+2);
  SIZE nz = particles_per_cell_length*(m_gridsize[2]+2);
  float incx = (m_bmax[0] - m_bmin[0] + 2*pad)/nx;
  float incy = (m_bmax[1] - m_bmin[1] + 2*pad)/ny;
  float incz = (m_bmax[2] - m_bmin[2] + 2*pad)/nz;

  SIZE i,j,k; // indices
  
  auto f = [&](SIZE i, SIZE j, SIZE k)
  {
    Vector3R<REAL> pos( m_bmin[0] - pad + i*incx, 
                        m_bmin[1] - pad + j*incy,
                        m_bmin[2] - pad + k*incz );
    m_grid( get_voxel_index(pos) ).boundvec.push_back(ParticleR<REAL>(pos));
    //fprintf(stderr, "v %f %f %f\n", pos[0], pos[1], pos[2]);
  };

  k = 0;
  for (i = 0; i < nx; ++i)
    for (j = 0; j < ny; ++j)
      f(i,j,k);
  k = nz;
  for (i = nx; i > 0; --i)
    for (j = ny; j > 0; --j)
      f(i,j,k);

  j = 0;
  for (i = nx; i > 0; --i)
    for (k = nz; k > 0; --k)
      f(i,j,k);
  j = ny;
  for (i = 0; i < nx; ++i)
    for (k = 0; k < nz; ++k)
      f(i,j,k);

  i = 0;
  for (j = 0; j < ny; ++j)
    for (k = nz; k > 0; --k)
      f(i,j,k);
  i = nx;
  for (j = ny; j > 0; --j)
    for (k = 0; k < nz; ++k)
      f(i,j,k);

  // two points not filled
  f(0, ny, nz);
  f(nx, 0, 0);
}

template<typename REAL, typename SIZE> template<int FT>
inline void UniformGridRS<REAL,SIZE>::compute_accelT()
{ compute_fluid_quantity< CFAccelRST<REAL, SIZE, FT>, FT >(); }

template<typename REAL, typename SIZE> template<int FT>
inline void UniformGridRS<REAL,SIZE>::compute_surface_normalT()
{ compute_fluid_quantity< CFSurfaceNormalRST<REAL, SIZE, FT>, FT >(); }

template<typename REAL, typename SIZE> template<int FT>
inline void UniformGridRS<REAL,SIZE>::compute_surface_tensionT()
{ compute_fluid_quantity< CFSurfaceTensionRST<REAL, SIZE, FT>, FT >(); }

template<typename REAL, typename SIZE> template<int FT>
inline void UniformGridRS<REAL,SIZE>::compute_densityT()
{ compute_fluid_quantity< CFDensityRST<REAL,SIZE,FT>, FT >(); }

template<typename REAL, typename SIZE> template<int FT>
inline void UniformGridRS<REAL,SIZE>::compute_density_updateT()
{ compute_fluid_quantity< CFDensityUpdateRST<REAL,SIZE,FT>, FT >(); }

template<typename REAL, typename SIZE>
template<typename ProcessPairFunc, typename ParticleType>
inline void UniformGridRS<REAL,SIZE>::compute_quantity()
{
//  clock_t s = clock();
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

  //qDebug() << "average time" << (float(clock() - s) / CLOCKS_PER_SEC);
}

template<typename REAL, typename SIZE, int FT>
inline void FTiter<REAL,SIZE,FT>::update_density(UniformGridRS<REAL,SIZE> &g,float timestep)
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
inline void FTiter<REAL,SIZE,FT>::compute_density(UniformGridRS<REAL,SIZE> &g)
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
inline void FTiter<REAL,SIZE,FT>::compute_accel(UniformGridRS<REAL,SIZE> &g)
{
  if (g.m_fluids[FT].size())
  {
    for (auto &fl : g.m_fluids[FT])
      fl->template cast<FT>()->reset_accel();     // now may assume all accelerations are zero

    g.template compute_accelT<FT>();
  }

  FTiter<REAL,SIZE,FT+1>::compute_accel(g); // recurse
}

template<typename REAL, typename SIZE>
inline void UniformGridRS<REAL,SIZE>::jacobi_pressure_solve(float dt)
{
  if (!m_fluids[ICS13].size())
    return;
  
  compute_fluid_quantity< CFPrepareJacobiRST<REAL, SIZE, ICS13>, ICS13 >();
  compute_fluid_quantity< CFJacobiSolveFirstRST<REAL, SIZE, ICS13>, ICS13 >();
  compute_fluid_quantity< CFJacobiSolveSecondRST<REAL, SIZE, ICS13>, ICS13 >();
  compute_fluid_quantity< CFJacobiSolveFirstRST<REAL, SIZE, ICS13>, ICS13 >();
  compute_fluid_quantity< CFJacobiSolveSecondRST<REAL, SIZE, ICS13>, ICS13 >();
  compute_fluid_quantity< CFPressureAccelRST<REAL, SIZE, ICS13>, ICS13 >();

  for (auto &fl : m_fluids[ICS13])
    fl->get_vel() = fl->get_vel() + dt*fl->get_accel();
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::run()
{
  if (m_num_fluids < 1)
    return;

  // timestep
  float dt = 1.0f/(global::dynset.fps * global::dynset.substeps);

  glprintf_tr("step: %.2es\n", dt);
  float rdt = dt;

  bool all_cached = check_and_write_hash();

  if (all_cached) // try to load cached frames into the fluid objects
    for (int j = 0; j < NUMTYPES; ++j)
      for (auto &fl : m_fluids[j])
        all_cached &= fl->load_saved_cache();
  
  if (!all_cached)
  {
    // Initialize the fluid for init_steps steps before simulating
    // temporarily disable gravity
    Vector3f grav = global::dynset.gravity;
    global::dynset.gravity = Vector3f(0.0,0.0,0.0);
    for (unsigned int iter = 0; iter < global::dynset.init_steps; ++iter)
    { // for each simulation substep
      step(dt, iter == 0);
    } // for each substep

    global::dynset.gravity = grav; // restore gravity

    if (m_stop_requested)
      return;

    for (int j = 0; j < NUMTYPES; ++j)
      for (auto &fl : m_fluids[j])
        fl->cache(0);
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
    all_cached = true;
    for (int j = 0; j < NUMTYPES; ++j)
      for (auto &fl : m_fluids[j])
        all_cached &= fl->is_cached(frame);

    if (all_cached)
    {
      for (int j = 0; j < NUMTYPES; ++j)
        for (auto &fl : m_fluids[j])
          fl->load_cached(frame);

      for (int j = 0; j < NUMTYPES; ++j)
        for (auto &fl : m_fluids[j])
          fl->update_data(); // notify gl we have new positions

      file_reads += 1;
      file_read_t += float(clock() - s) / CLOCKS_PER_SEC;

      boost::chrono::nanoseconds elapsed = (real_clock::now() - start_time);
      boost::chrono::nanoseconds perframe(boost::int_least64_t(1.0e9/global::dynset.fps));
      if (elapsed < perframe)
        std::this_thread::sleep_for(std::chrono::nanoseconds((perframe - elapsed).count()));
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

      for (int j = 0; j < NUMTYPES; ++j)
        for (auto &fl : m_fluids[j])
          fl->cache(frame); // cache frame

    } // if not cached

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

template<typename REAL, typename SIZE>
inline bool UniformGridRS<REAL,SIZE>::step(float dt, bool first_step, float *substep_t)
{
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

  update_grid(); // prepare grid for simulation step

  clock_t prev_t = clock();

  //if (iter)
  // FTiter<REAL,SIZE,DEFAULT>::update_density(*this, dt);
  //else
  FTiter<REAL,SIZE,DEFAULT>::compute_density(*this);

  FTiter<REAL,SIZE,DEFAULT>::compute_accel(*this); // update m_accel

  if (m_fluids[MCG03].size()) // extra steps to compute surface tension for MCG03
  {
    compute_surface_normalT<MCG03>();
    compute_surface_tensionT<MCG03>();
  }

  if (substep_t)
    *substep_t += float(clock() - prev_t) / CLOCKS_PER_SEC;

  float factor = first_step ? 0.5f : 1.0f; // leap-frog method has a different first step

  for (int j = 0; j < NUMTYPES; ++j)
    for (auto &fl : m_fluids[j])
      fl->get_vel() = fl->get_vel() + factor*dt*fl->get_accel();

  jacobi_pressure_solve(dt); // only relevant for IISPH fluids

  for (int j = 0; j < NUMTYPES; ++j)
  {
    for (auto &fl : m_fluids[j])
    {
      fl->get_pos() = (fl->get_pos() + dt*fl->get_vel()).eval();
      if (j == MCG03)
        fl->resolve_collisions();
      //    else
      //      fl->clamp(0.0f);

      // prepare velocities for acceleration computation in next step
      fl->get_vel() = fl->get_vel() + 0.5*dt*fl->get_accel();
    }
  }

  if (m_stop_requested) return false;

  for (int j = 0; j < NUMTYPES; ++j)
    for (auto &fl : m_fluids[j])
      fl->update_data(); // notify gl we have new positions

  return true;
}


template<typename REAL, typename SIZE>
inline bool UniformGridRS<REAL,SIZE>::check_and_write_hash() 
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

template class UniformGridRS<double, unsigned int>;
