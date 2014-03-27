#include <ctime>
#include <algorithm>
#include <limits>
#include "dynamics.h"
#include "glpointcloud.h"
#include "gltext.h"
#include "fluid.h"


// UniformGrid stuff

template<typename REAL, typename SIZE>
UniformGridRS<REAL,SIZE>::UniformGridRS(
    const Vector3f &bmin,
    const Vector3f &bmax ) 
  : m_h(2e-3f) // minimum possible cell size
  , m_hinv(500.0f)
  , m_bmin(bmin) // smallest boundary corner
  , m_bmax(bmax) // largest boundary corner
  , m_stop_requested(false)
{ }

template<typename REAL, typename SIZE>
UniformGridRS<REAL,SIZE>::~UniformGridRS() { }
  
template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::add_fluid( FluidRS<REAL,SIZE> *fl )
{
  m_fluids.push_back(fl);
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::init()
{
  // set cell size to be the maximum of the kernel support over all fluids;
  for ( auto &fl : m_fluids )
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
  
  populate_fluid_data(); // populate grid with particles
  populate_bound_data();

  int count = 0;
  for (auto &fl : m_fluids)
  {
    glprintf_trcv(fl->get_color(), "--- Fluid %d ---\n", ++count);
    glprintf_trcv(fl->get_color(), "# of particles: %.2f  \n", fl->get_num_vertices());
    glprintf_trcv(fl->get_color(), "rest density: %.2f  \n", fl->get_rest_density());
    glprintf_trcv(fl->get_color(), "particle mass: %.2f  \n", fl->get_mass());
    glprintf_trcv(fl->get_color(), "viscosity: %.2f  \n", fl->get_viscosity());
    glprintf_trcv(fl->get_color(), "surface tension: %.2f  \n", fl->get_surface_tension());
    glprintf_trcv(fl->get_color(), "sound speed: %.2e  \n", std::sqrt(fl->get_sound_speed2()));
  }

  CBVolumeR<REAL> proc;
  proc.init_kernel(m_h);
  m_bound_volume_proc = &proc;
  compute_bound_quantity< CBVolumeR<REAL> >();
  m_bound_volume_proc = NULL;
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::update_grid()
{
  clear_fluid_data();
  populate_fluid_data(); // populate grid with particles
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::clear_fluid_data()
{
  for (SIZE i = 0; i < m_gridsize[0]; ++i)
    for (SIZE j = 0; j < m_gridsize[1]; ++j)
      for (SIZE k = 0; k < m_gridsize[2]; ++k)
        m_grid[i][j][k].fluidvec.clear();
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::populate_fluid_data()
{ 
  unsigned short id = 0; // generate ids
  for ( auto &fl : m_fluids )
  {
    SIZE num_vtx = fl->get_num_vertices();
    for ( SIZE i = 0; i < num_vtx; ++i )
    {
      Vector3R<REAL> pos( fl->get_pos().col(i) );
      Vector3R<REAL> vel( fl->get_vel().col(i) );
      Array3Index idx = get_voxel_index(pos);
      m_grid(idx).fluidvec.push_back( 
          FluidParticleR<REAL>(
            pos, vel, fl->accel_at(i), fl->extern_accel_at(i), id));
    }
    id++;
  }
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

#if 0
// heuristic to approximate the initial density
// currently we use the given density and use Shepard's method 
// to correct the kernel
template<typename REAL>
class ComputeInitialDensityR
{
public:
  ComputeInitialDensityR(float h, REAL m, REAL &rd) 
    : kern(h), mass(m), rest_density(rd) { }
  ~ComputeInitialDensityR() { }

  inline void init_particle(ParticleR<REAL> &p)
  {
//    qDebug() << "v" << p.pos[0] << p.pos[1] << p.pos[2];
    p.dinv = 0.0f;
  }

  inline void fluid(ParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    p.dinv += kern[ p.pos - near_p.pos ];
  }
  inline void bound(ParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    p.dinv += kern[ p.pos - near_p.pos ];
  }

  inline void finish_particle(ParticleR<REAL> &p)
  {
    REAL density = p.dinv * mass * kern.coef;
    p.dinv = 1.0f/density;
//    qDebug() << density;
    rest_density += density;
  }

private:
  Poly6Kernel kern; // used to compute pressure force
  REAL mass;
  REAL &rest_density;
};

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::compute_initial_density()
{
  ComputeInitialDensityR<REAL> proc(m_h, m_dpc->m_mass, m_dpc->m_rest_density);
  compute_bound_quantity(proc);

  m_dpc->m_rest_density = 0.0f;
  compute_fluid_quantity(proc);
  m_dpc->m_rest_density /= m_dpc->m_num_vertices;
}
#endif


template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::update_density(float timestep)
{
  for ( auto &fl : m_fluids )
    fl->m_fluid_density_update_proc.init( timestep );

  compute_fluid_quantity< CFDensityUpdateRS<REAL, SIZE> >();
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::compute_density()
{
  REAL max_var[m_fluids.size()];
  REAL avg_var[m_fluids.size()];
  
  int i = 0;
  for ( auto &fl : m_fluids )
  {
    max_var[i] = avg_var[i] = 0.0f;
    fl->m_fluid_density_proc.init( max_var[i], avg_var[i] );
    i++;
  }

  compute_fluid_quantity< CFDensityRS<REAL,SIZE> >();

#ifdef REPORT_DENSITY_VARIATION
  i = 0;
  for ( auto &fl : m_fluids )
  {
    avg_var[i] = avg_var[i]/fl->get_num_vertices();
    qDebug("Fluid %d:  max: %.0f, %.1f percent;    avg: %.0f, %.1f percent", i,
        max_var[i], 100.00f*max_var[i]/fl->get_rest_density(), 
        avg_var[i], 100.00f*avg_var[i]/fl->get_rest_density());
    i++;
  }
#endif
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::compute_pressure()
{
  // will be needed when we use alternative pressure methods (pressure projection)
  //compute_fluid_quantity(m_proc_fluid_density);

  compute_fluid_quantity< CFPressureRS<REAL, SIZE> >();
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::compute_accel()
{
  for ( auto &fl : m_fluids )
    fl->reset_accel();     // now may assume all accelerations are zero

  compute_fluid_quantity< CFPressureAccelRS<REAL, SIZE> >();
  compute_fluid_quantity< CFViscosityAccelRS<REAL, SIZE> >();
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
        std::vector<ParticleType> &pvec = 
          m_grid[i][j][k].template get_vec<ParticleType>();
        if (pvec.empty())
          continue;

        IndexRange xrange = range3(i,nx);
        IndexRange yrange = range3(j,ny);
        IndexRange zrange = range3(k,nz);
        GridView neigh_view = m_grid[ boost::indices[xrange][yrange][zrange] ];

        for ( auto &p : pvec )  // prepare data
        {
          ProcessPairFunc &process = determine_proc< ProcessPairFunc >(p);
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
                ProcessPairFunc &process = determine_proc< ProcessPairFunc >(p);
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
          ProcessPairFunc &process = determine_proc< ProcessPairFunc >(p);
          process.finish_particle(p);
        }
      } // for k
    } // for j
  } // for i

}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::run()
{
  if (m_fluids.size() < 1)
    return;

  init();

#ifdef MCG03
  float dt = 2.5e-3;
#else
  float dt = 4.52e-4;
#endif

  glprintf_tr("step: %.2es\n", dt);
  float rdt = dt;

  compute_density();
  compute_pressure();
  compute_accel(); // update m_accel

  for ( auto &fl : m_fluids )
    fl->get_vel() = fl->get_vel() + 0.5*dt*fl->get_accel(); // initial half velocity

  clock_t prev_t = clock();
  float t = 0.0f;
  unsigned int count=0;
  for ( ; ; )
  {
    ///////// testing adaptive time step
    float fdt = std::numeric_limits<float>::infinity();
    float cvdt = std::numeric_limits<float>::infinity();
    for (auto &fl : m_fluids)
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

    for (auto &fl : m_fluids)
      fl->get_pos() = (fl->get_pos() + dt*fl->get_vel()).eval();
    
    for (auto &fl : m_fluids)
      fl->resolve_collisions();

    update_grid();

    for (auto &fl : m_fluids)
      fl->update_data(); // notify gl we have new positions

    if (m_stop_requested)
      break;

    for (auto &fl : m_fluids)
      fl->get_vel() = fl->get_vel() + 0.5*dt*fl->get_accel();

    //if (count % 100)
    //  update_density(dt);
    //else
      compute_density();
    compute_pressure();
    compute_accel(); // update m_accel

    for (auto &fl : m_fluids)
      fl->get_vel() = fl->get_vel() + dt*fl->get_accel();

//    if (count == 0)
//    {
//      for (int i = 0; i < 10; ++i)
//        qDebug() << m_accel.col(i)[0] << m_accel.col(i)[1] << m_accel.col(i)[2] ;
//      return;
//    }

    clock_t cur_t = clock();
    t += float(cur_t - prev_t) / CLOCKS_PER_SEC;
    prev_t = cur_t;
    count += 1;
  }
  glprintf_tr("avg time per step: %.2es\n", t / float(count));
}

template class UniformGridRS<double, unsigned int>;
