#include "dynamics.h"
#include "glpointcloud.h"
#include "gltext.h"
#include <ctime>
#include <algorithm>
#include <limits>

// UniformGrid stuff

template<typename REAL, typename SIZE>
UniformGridRS<REAL,SIZE>::UniformGridRS( const Vector3f &bmin, const Vector3f &bmax) 
  : m_h(2e-3f) // minimum possible cell size
  , m_hinv(500.0f)
  , m_bmin(bmin) // smallest boundary corner
  , m_bmax(bmax) // largest boundary corner
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
    if (m_h < fl->get_kernel_radius())
      m_h = fl->get_kernel_radius();

  m_hinv = 1.0f/m_h;

  // determine the number of cells needed
  m_gridsize = {{
    static_cast<Index>((m_bmax[0] - m_bmin[0])*m_hinv),
    static_cast<Index>((m_bmax[1] - m_bmin[1])*m_hinv),
    static_cast<Index>((m_bmax[2] - m_bmin[2])*m_hinv) }};
  glprintf_tr("grid size: %d %d %d\n", m_gridsize[0], m_gridsize[1], m_gridsize[2]);

  m_grid.resize( boost::extents[m_gridsize[0]][m_gridsize[1]][m_gridsize[2]] );
  
  populate_fluid_data(); // populate grid with particles
  populate_bound_data();

  //compute_initial_density();

  int count = 0;
  for (auto &fl : m_fluids)
  {
    glprintf_tr("Fluid %d:\n", ++count);
    glprintf_tr("rest density: %.2f  \n", fl->get_rest_density());
    glprintf_tr("particle mass: %.2f  \n", fl->get_mass());
    glprintf_tr("viscosity: %.2f  \n", fl->get_viscosity());
    glprintf_tr("surface tension: %.2f  \n", fl->get_surface_tension());
    glprintf_tr("sound speed: %.2e  \n", std::sqrt(fl->get_sound_speed2()));
  }

  compute_bound_volume();
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::update()
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
  for ( auto &fl : m_fluids )
  {
    SIZE num_vtx = fl->get_num_vertices();
    for ( SIZE i = 0; i < num_vtx; ++i )
    {
      Vector3R<REAL> pos( fl->get_pos().col(i) );
      Vector3R<REAL> vel( fl->get_vel().col(i) );
      Array3Index idx = get_voxel_index(pos);
      m_grid(idx).fluidvec.push_back( 
          FluidParticleR<REAL>(pos, vel, fl->accel_at(i), fl->extern_accel_at(i))
          );
    }
  }
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::populate_bound_data()
{
  // walk through all boundary cells
  SIZE nx = m_gridsize[0];
  SIZE ny = m_gridsize[1];
  SIZE nz = m_gridsize[2];
  float h = 1/m_hinv; // TODO: save and use actual gridsize
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

  // check for consistency
  #if 0
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


template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::update_density(float timestep)
{
  ComputeFluidDensityUpdateR<REAL> proc(m_h);
  proc.init(m_dpc->m_mass, timestep );
  compute_fluid_quantity(proc);
}
template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::compute_density()
{
  ComputeFluidDensityR<REAL> proc(m_h);
  REAL max_var = 0.0f;
  REAL avg_var = 0.0f;
  proc.init(m_dpc->m_mass, m_dpc->m_rest_density, m_dpc->get_radius(), max_var, avg_var );
  compute_fluid_quantity(proc);
  avg_var = avg_var/m_dpc->m_num_vertices;
  qDebug("max: %.0f, %.1f percent;    avg: %.0f, %.1f percent",
      max_var, 100.00f*max_var/m_dpc->m_rest_density, 
      avg_var, 100.00f*avg_var/m_dpc->m_rest_density);
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::compute_pressure()
{
  // will be needed when we use alternative pressure methods (pressure projection)
  //compute_fluid_quantity(m_proc_fluid_density);

  compute_fluid_quantity(m_proc_pressure);
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::compute_accel()
{
  m_dpc->reset_accel();     // now may assume all accelerations are zero
  compute_fluid_quantity( m_proc_pressure_accel );
  compute_fluid_quantity( m_proc_viscosity_accel );
}

// TODO: refactor the following two routines
template<typename REAL, typename SIZE>
template<typename ProcessPairFunc>
void UniformGridRS<REAL,SIZE>::compute_bound_quantity(ProcessPairFunc process)
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
        StaticParticles &boundvec = m_grid[i][j][k].boundvec;
        if (boundvec.empty())
          continue;

        IndexRange xrange = range3(i,nx);
        IndexRange yrange = range3(j,ny);
        IndexRange zrange = range3(k,nz);
        GridView neigh_view = m_grid[ boost::indices[xrange][yrange][zrange] ];

        for ( ParticleR<REAL> &p : boundvec )  // prepare data
          process.init_particle(p);

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
              for ( ParticleR<REAL> &p : boundvec )
              {
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

        for ( ParticleR<REAL> &p : boundvec )  // finalize data
          process.finish_particle(p);
      } // for k
    } // for j
  } // for i

}

template<typename REAL, typename SIZE>
template<typename ProcessPairFunc>
void UniformGridRS<REAL,SIZE>::compute_fluid_quantity(ProcessPairFunc process)
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
        DynamicParticles &fluidvec = m_grid[i][j][k].fluidvec;
        if (fluidvec.empty())
          continue;

        IndexRange xrange = range3(i,nx);
        IndexRange yrange = range3(j,ny);
        IndexRange zrange = range3(k,nz);
        GridView neigh_view = m_grid[ boost::indices[xrange][yrange][zrange] ];

        for ( FluidParticleR<REAL> &p : fluidvec )  // prepare data
          process.init_particle(p);

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
              for ( FluidParticleR<REAL> &p : fluidvec )
              {
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

        for ( FluidParticleR<REAL> &p : fluidvec )  // finalize data
          process.finish_particle(p);
      } // for k
    } // for j
  } // for i
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::run()
{
  init();

  float dt = 4.52e-4;
  glprintf_tr("step: %.2es\n", dt);
  float rdt = dt;

  m_grid.compute_density();
  m_grid.compute_pressure();
  m_grid.compute_accel(); // update m_accel

  m_vel = m_vel + 0.5*dt*m_accel; // initial half velocity

  clock_t prev_t = clock();
  float t = 0.0f;
  unsigned int count=0;
  for ( ; ; )
  {
    // testing adaptive time step
    float fdt = std::numeric_limits<float>::infinity();
    float cvdt = std::numeric_limits<float>::infinity();
    for (SIZE i = 0; i < m_num_vertices; ++i)
    {
      fdt = std::min(fdt, float(this->get_radius()/(m_mass*Vector3R<REAL>(m_extern_accel.col(i)).norm())));
      cvdt = std::min(cvdt, float(this->get_radius()/(m_c2*(1+0.6*0.08))));
      // check for nan and infs
      if (!std::isfinite(m_vel.col(i)[0]) || !std::isfinite(m_vel.col(i)[1]) || !std::isfinite(m_vel.col(i)[2]))
        qDebug() << "Found NaN at " << i;
    }
    float new_rdt = std::min(0.25*fdt, 0.4*cvdt);
    if (new_rdt != rdt)
    {
      glprintf_tl("\rrecommended step: %.2es", new_rdt);
      rdt = new_rdt;
    }

    // end of adaptive timestep test

    this->m_pos = (this->m_pos + dt*m_vel).eval();
    
    //resolve_collisions();
    m_grid.update();

    if (m_glpc)
      m_glpc->update_data(); // notify gl we have new positions

    if (m_stop_requested)
      break;

    m_vel = m_vel + 0.5*dt*m_accel;

    if (count % 100)
      m_grid.update_density(dt);
    else
      m_grid.compute_density();
    m_grid.compute_pressure();
    m_grid.compute_accel(); // update m_accel

    m_vel = m_vel + dt*m_accel;
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


typedef UniformGridRS<double, unsigned int> UniformGrid;
template class UniformGrid;
