#include <ctime>
#include <algorithm>
#include <limits>
#include <fstream>
#include <iomanip>
#include "sph.h"
#include "glpointcloud.h"
#include "gltext.h"
//#include "quantityprocessor.h"
#include "fluid.h"
#include "settings.h"
#define BOOST_CHRONO_HEADER_ONLY
#include <boost/chrono.hpp>

//#define REPORT_DENSITY_VARIATION

typedef boost::chrono::process_real_cpu_clock real_clock;
typedef boost::chrono::nanoseconds ns_t;
typedef real_clock::time_point real_t;

// From [Solenthaler and Pajarola 2008], an alternative density
// (number_density * mass) is used
/*
class CBVolume : public CBQ<CBVolume>
{
public:
  inline void init_kernel(float h) { m_kern.init(h); }
  inline void init_particle(Particle &p) 
  { p.dinv = 0.0f; }
  inline void fluid(Particle &p, FluidParticle &near_p)
  { }
  inline void bound(Particle &p, Particle &near_p)
  {
    p.dinv += this->m_kern[ p.pos - near_p.pos ];
  }
  inline void finish_particle(Particle &p)
  {
    p.dinv = 1.0f/(p.dinv * this->m_kern.coef);
  }
private:
  CubicSplineKernel m_kern;
};
*/

// SPHGrid stuff


SPHGrid::SPHGrid(
    const Vector3f &bmin,
    const Vector3f &bmax,
    DynamicsManager &dynman ) 
  : m_dynman(dynman)
  , m_h(2e-3f) // minimum possible cell size
  , m_hinv(500.0f)
  , m_bmin(bmin) // smallest boundary corner
  , m_bmax(bmax) // largest boundary corner
  , m_stop_requested(false)
  , m_pause(false)
{ }


SPHGrid::~SPHGrid() { }

void
SPHGrid::init(MaterialManager &matman)
{
  if (m_dynman.get_numfluids() < 1)
    return;

  // set cell size to be the maximum of the kernel support over all fluids;
  m_h = std::max(m_dynman.get_max_radius(), m_h);

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
  
  m_dynman.populate_sph_grid(*this); // populate grid with particles

  //CBVolume proc;
  //proc.init_kernel(m_h);
  //m_bound_volume_proc = &proc;
  //compute_bound_quantity< CBVolume >();
  //m_bound_volume_proc = NULL;
}


inline void SPHGrid::update_grid()
{
  clear_fluid_data();
  populate_fluid_data(); // populate grid with particles
}


inline void SPHGrid::clear_fluid_data()
{
  for (Size i = 0; i < m_gridsize[0]; ++i)
    for (Size j = 0; j < m_gridsize[1]; ++j)
      for (Size k = 0; k < m_gridsize[2]; ++k)
        m_grid[i][j][k].clear();
}

#if 0
void SPHGrid::populate_bound_data()
{
  unsigned char particles_per_cell_length = 4;

  // place boundary particles slightly away from the actual boundary
  float pad = 0.5*m_h;

  Size nx = particles_per_cell_length*(m_gridsize[0]+2);
  Size ny = particles_per_cell_length*(m_gridsize[1]+2);
  Size nz = particles_per_cell_length*(m_gridsize[2]+2);
  float incx = (m_bmax[0] - m_bmin[0] + 2*pad)/nx;
  float incy = (m_bmax[1] - m_bmin[1] + 2*pad)/ny;
  float incz = (m_bmax[2] - m_bmin[2] + 2*pad)/nz;

  Size i,j,k; // indices
  
  auto f = [&](Size i, Size j, Size k)
  {
    Vector3R<Real> pos( m_bmin[0] - pad + i*incx, 
                        m_bmin[1] - pad + j*incy,
                        m_bmin[2] - pad + k*incz );
    m_grid( get_voxel_index(pos)
        ).get_pvec<STATIC>().push_back(ParticleT<STATIC>(pos));
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
#endif

template<int FT>
inline void SPHGrid::compute_accelT()
{ 
//  compute_fluid_quantity< CFAccelT<FT>, FT >(); 
}

template<int FT>
inline void SPHGrid::compute_surface_normalT()
{ 
//  compute_fluid_quantity< CFSurfaceNormalT<FT>, FT >(); 
}

template<int FT>
inline void SPHGrid::compute_surface_tensionT()
{ 
//  compute_fluid_quantity< CFSurfaceTensionT<FT>, FT >(); 
}

template<int FT>
inline void SPHGrid::compute_densityT()
{ 
 // compute_fluid_quantity< CFDensityT<FT>, FT >(); 
}

template<int FT>
inline void SPHGrid::compute_density_updateT()
{
//  compute_fluid_quantity< CFDensityUpdateT<FT>, FT >(); 
}

// variadic base case
template<int F, typename ParticleType>
inline void SPHGrid::interact_with_neigh_cell(
    ParticleType &p, Cell &cell)
{ }

// variadic induction
template<int F, typename ParticleType, int NPT, int... NPTs>
inline void SPHGrid::interact_with_neigh_cell(
    ParticleType &p, Cell &cell)
{
  auto &neigh_pvec = cell.template get_pvec<NPT>();
  for ( auto &neigh_p : neigh_pvec )
    p.neigh<F>(neigh_p); // process neighbouring particles

  interact_with_neigh_cell<F,ParticleType, NPTs...>(p,cell);
}

// variadic base case
template<int F>
inline void SPHGrid::compute_quantity_in_cell(
    Size i, Size j, Size k,
    Size nx, Size ny, Size nz)
{ }

// variadic induction
template<int F, int PT, int... PTs>
inline void SPHGrid::compute_quantity_in_cell(
    Size i, Size j, Size k,
    Size nx, Size ny, Size nz)
{
  auto &pvec = m_grid[i][j][k].template get_pvec<PT>();
  if (pvec.empty())
    return;

  IndexRange xrange = range3(i,nx);
  IndexRange yrange = range3(j,ny);
  IndexRange zrange = range3(k,nz);
  // TODO: why do we need the GridView? just iterate directly in the grid
  GridView neigh_view = m_grid[ boost::indices[xrange][yrange][zrange] ];

  for ( auto &p : pvec )  // prepare data
    p.init<F>();

  Size xrange_size = xrange.finish() - xrange.start();
  Size yrange_size = yrange.finish() - yrange.start();
  Size zrange_size = zrange.finish() - zrange.start();
  for (Size near_i = 0; near_i < xrange_size; ++near_i)
  {
    for (Size near_j = 0; near_j < yrange_size; ++near_j)
    {
      for (Size near_k = 0; near_k < zrange_size; ++near_k)
      {
        Cell &cell = neigh_view[near_i][near_j][near_k];
        for ( auto &p : pvec )
          interact_with_neigh_cell<F,PT,decltype(p),ALL_PARTICLE_TYPES>(p,cell);
      }
    }
  }

  for ( auto &p : pvec )  // finalize data
    p.finish<F>();

  compute_quantity_in_cell<F,PTs...>(i,j,k,nx,ny,nz);
}

template<typename F, int... PT>
inline void SPHGrid::compute_quantity()
{
//  clock_t s = clock();
  Size nx = m_gridsize[0];
  Size ny = m_gridsize[1];
  Size nz = m_gridsize[2];
  for (Size i = 0; i < nx; ++i)
  {
    for (Size j = 0; j < ny; ++j)
    {
      for (Size k = 0; k < nz; ++k)
      {
        compute_quantity_in_cell<F,PT...>(i,j,k,nx,ny,nz);
      } // for k
    } // for j
  } // for i

  //qDebug() << "average time" << (float(clock() - s) / CLOCKS_PER_SEC);
}

#if 0
template<int FT>
inline void FTiter<FT>::update_density(SPHGrid &g,float timestep)
{
  if (g.m_fluids[FT].size())
  {
    for (auto &fl : g.m_fluids[FT])
      fl->template cast<FT>()->m_fluid_density_update_proc.init( timestep ); 

    g.template compute_density_updateT<FT>();
  }
  FTiter<FT+1>::update_density(g,timestep); // recurse
}
#endif

template<int PT>
inline void SPHGrid::compute_density()
{
  FluidDataVecT<PT> &fldatavec = m_dynman.get_fluiddatas<PT>();
  if (!fldatavec.size())
    return;

#ifdef REPORT_DENSITY_VARIATION
  Real max_var[fldatavec.size()];
  Real avg_var[fldatavec.size()];
  
  int i = 0;
  for ( auto &fldata : fldatavec )
  {
    max_var[i] = avg_var[i] = 0.0f; // initialize
    fldata.m_max_var = &max_var[i];
    fldata.m_avg_var = &avg_var[i];
    i++;
  }
#endif

  compute_quantity<Density,PT>();

#ifdef REPORT_DENSITY_VARIATION
  i = 0;
  for (auto &fldata : fldatavec )
  {
    avg_var[i] = avg_var[i]/fldata.fl.get_num_vertices();
    qDebug("Fluid %d:  max: %.0f, %.1f percent;    avg: %.0f, %.1f percent", i,
        max_var[i], 100.00f*max_var[i]/fldata.fl.get_rest_density(), 
        avg_var[i], 100.00f*avg_var[i]/fldata.fl.get_rest_density());
    i++;
  }
#endif
}

template<int PT, int PT2, int... PTs>
inline void SPHGrid::compute_density()
{
  compute_density<PT>();       // do one
  compute_density<PT2, PTs>(); // recurse on the rest
}

template<int FT>
inline void FTiter<FT>::compute_accel(SPHGrid &g)
{
  if (g.m_fluids[FT].size())
  {
    for (auto &fl : g.m_fluids[FT])
      fl->template cast<FT>()->reset_accel();     // now may assume all accelerations are zero

    g.template compute_accelT<FT>();
  }

  FTiter<FT+1>::compute_accel(g); // recurse
}

#if 0
inline void SPHGrid::jacobi_pressure_solve(float dt, float factor)
{
  if (!m_fluids[ICS13].size())
    return;
  
  for (auto &fl : m_fluids[ICS13])
  {
    //fl->template cast<ICS13>()->m_fluid_prepare_jacobi_proc.init( dt );
    //fl->template cast<ICS13>()->m_fluid_jacobi_solve1_proc.init( dt );
    //fl->template cast<ICS13>()->m_fluid_jacobi_solve2_proc.init( dt, fl->get_avg_density(), fl->get_avg_pressure() );
    //fl->template cast<ICS13>()->reset_extern_accel();     // now may assume all accelerations are zero
  }

  //compute_fluid_quantity< CFPrepareJacobiT<ICS13>, ICS13 >();

  bool proceed = false;
  int iter = 0;
  do
  {
    //compute_fluid_quantity< CFJacobiSolveFirstT<ICS13>, ICS13 >();
    //compute_fluid_quantity< CFJacobiSolveSecondT<ICS13>, ICS13 >();
    proceed = false;
    for (auto &fl : m_fluids[ICS13])
    {
      proceed |= std::abs(fl->get_avg_density()/fl->get_num_vertices() -
          fl->get_rest_density()) > 1;
      qDebug() << "avg density  = " << fl->get_avg_density()/fl->get_num_vertices();
      qDebug() << "avg pressure = " << fl->get_avg_pressure()/fl->get_num_vertices();
      fl->get_avg_density() = 0.0f;
      fl->get_avg_pressure() = 0.0f;
    }
  } while(proceed);

  //compute_fluid_quantity< CFPressureAccelT<ICS13>, ICS13 >();

  for (auto &fl : m_fluids[ICS13])
    fl->get_vel() = fl->get_vel() + factor*dt*fl->get_extern_accel();
}
#endif


void SPHGrid::run()
{
  if (m_dynman.get_numfluids() < 1)
    return;

  // timestep
  float dt = 1.0f/(global::dynset.fps * global::dynset.substeps);

  glprintf_tr("step: %.2es\n", dt);
  float rdt = dt;

  bool all_cached = check_and_write_hash();

  if (all_cached) // try to load cached frames into the fluid objects
    all_cached &= m_dynman.load_saved_cache();
  
  if (!all_cached)
  {
    // Initialize the fluid for init_steps steps before simulating
    // temporarily disable gravity
    Vector3f grav = global::dynset.gravity;
    global::dynset.gravity = Vector3f(0.0,0.0,0.0);
    for (unsigned int iter = 0; iter < global::dynset.init_steps; ++iter)
    { // for each simulation substep
      if (!m_dynman.step(dt, iter == 0))
        break;
    } // for each substep

    global::dynset.gravity = grav; // restore gravity

    if (m_stop_requested)
      return;

    m_dynman.cache(0);
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
    cached = m_dynman.is_cached(frame);

    if (cached)
    {
      m_dynman.load_cached(frame); // load cached frame
      m_dynman.prepare_vis_data(); // notify gl we have new positions

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
        if (!m_dynman.step(dt, first_step, &substep_t))
          break;
      } // for each substep

      m_dynman.cache(frame);
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


inline bool SPHGrid::check_and_write_hash() 
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
