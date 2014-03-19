#include "dynamics.h"
#include "glpointcloud.h"
#include "gltext.h"
#include <ctime>

// UniformGrid stuff

template<typename REAL, typename SIZE>
UniformGridRS<REAL,SIZE>::UniformGridRS(
    DynamicPointCloudRS<REAL,SIZE> *dpc,
    float h,
    float grid_h,
    const Vector3f &bmin)
  : m_dpc(dpc)
  , m_h(h)
  , m_hinv(1.0f/grid_h)
  , m_proc_pressure(m_h)
  , m_proc_accel(m_h)
  , m_bmin(bmin)
{
  glprintf_tr("radius: %.2f\n", h);
}

template<typename REAL, typename SIZE>
UniformGridRS<REAL,SIZE>::UniformGridRS(DynamicPointCloudRS<REAL,SIZE> *dpc, float h, const Vector3f &bmin)
  : m_dpc(dpc)
  , m_h(h)
  , m_hinv(1.0f/h)
  , m_proc_pressure(m_h)
  , m_proc_accel(m_h)
  , m_bmin(bmin)
{
  glprintf_tr("radius: %.2f\n", h);
}

template<typename REAL, typename SIZE>
UniformGridRS<REAL,SIZE>::~UniformGridRS() 
{ }
  
template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::init()
{
  // determine the number of voxels needed
  const Vector3f &bmax = m_dpc->get_bmax();

  m_gridsize = {{
    static_cast<Index>((bmax[0] - m_bmin[0])*m_hinv),
    static_cast<Index>((bmax[1] - m_bmin[1])*m_hinv),
    static_cast<Index>((bmax[2] - m_bmin[2])*m_hinv) }};
  glprintf_tr("grid size: %d %d %d\n", m_gridsize[0], m_gridsize[1], m_gridsize[2]);

  m_grid.resize( boost::extents[m_gridsize[0]][m_gridsize[1]][m_gridsize[2]] );
  
  populate_fluid_data(); // populate grid with particles
  populate_bound_data();

  compute_initial_density();

  glprintf_tr("mass: %.2f\n", m_dpc->m_mass);

  m_proc_pressure.init(m_dpc->m_mass, m_dpc->m_rest_density, m_dpc->m_c2);
  m_proc_accel.init(m_dpc->m_mass, m_dpc->m_viscosity, m_dpc->m_st);
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
  SIZE num_vtx = m_dpc->get_num_vertices();
  for ( SIZE i = 0; i < num_vtx; ++i )
  {
    Vector3R<REAL> pos( m_dpc->m_pos.col(i) );
    Vector3R<REAL> vel( m_dpc->m_vel.col(i) );
    Array3Index idx = get_voxel_index(pos);
    m_grid(idx).fluidvec.push_back( DynamicParticleR<REAL>(pos, vel, m_dpc->accel_at(i)) );
  }
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::populate_bound_data()
{
  // walk through all boundary cells
  SIZE nx = m_gridsize[0];
  SIZE ny = m_gridsize[1];
  SIZE nz = m_gridsize[2];
  REAL h = 1/m_hinv; // TODO: save and use actual gridsize
  Vector3R<REAL> bmin = m_bmin.template cast<REAL>();

  SIZE i = 0;
  for (SIZE j = 0; j < ny; ++j)
    for (SIZE k = 0; k < nz; ++k)
    {
      StaticParticles &boundvec = m_grid[i][j][k].boundvec;
      boundvec.push_back( ParticleR<REAL>(h*Vector3R<REAL>(0, j, k) + bmin) );
      if (k == nz-1)
        boundvec.push_back( ParticleR<REAL>(h*Vector3R<REAL>(0, j, nz) + bmin) );
      if (j == ny-1)
        boundvec.push_back( ParticleR<REAL>(h*Vector3R<REAL>(0, ny, k) + bmin) );
      if (j == ny-1 && k == nz-1)
        boundvec.push_back( ParticleR<REAL>(h*Vector3R<REAL>(0, ny, nz) + bmin) );
    }

  for (i = 1; i < nx; ++i)
  {
    SIZE j = 0;
    for (SIZE k = 0; k < nz; ++k)
    {
      StaticParticles &boundvec = m_grid[i][j][k].boundvec;
      boundvec.push_back( ParticleR<REAL>(h*Vector3R<REAL>(i, 0, k) + bmin) );
      if (k == nz-1)
        boundvec.push_back( ParticleR<REAL>(h*Vector3R<REAL>(i, 0, nz) + bmin) );
    }

    for (j = 1; j < ny; ++j)
    {
      m_grid[i][j][0].boundvec.push_back( ParticleR<REAL>(h*Vector3R<REAL>(i, j, 0) + bmin) );

      m_grid[i][j][nz-1].boundvec.push_back( ParticleR<REAL>(h*Vector3R<REAL>(i, j, nz) + bmin) );
    }

    j = ny-1;
    for (SIZE k = 0; k < nz; ++k)
    {
      StaticParticles &boundvec = m_grid[i][j][k].boundvec;
      boundvec.push_back( ParticleR<REAL>(h*Vector3R<REAL>(i, ny, k) + bmin) );
      if (k == nz-1)
        boundvec.push_back( ParticleR<REAL>(h*Vector3R<REAL>(i, ny, nz) + bmin) );
    }

  } // for i

  i = nx-1; // cap off the top
  for (SIZE j = 0; j < ny; ++j)
    for (SIZE k = 0; k < nz; ++k)
    {
      StaticParticles &boundvec = m_grid[i][j][k].boundvec;
      boundvec.push_back( ParticleR<REAL>(h*Vector3R<REAL>(nx, j, k) + bmin) );
      if (k == nz-1)
        boundvec.push_back( ParticleR<REAL>(h*Vector3R<REAL>(nx, j, nz) + bmin) );
      if (j == ny-1)
        boundvec.push_back( ParticleR<REAL>(h*Vector3R<REAL>(nx, ny, k) + bmin) );
      if (j == ny-1 && k == nz-1)
        boundvec.push_back( ParticleR<REAL>(h*Vector3R<REAL>(nx, ny, nz) + bmin) );
    }
}



template<typename REAL>
class ProcessInitialDensityR
{
public:
  ProcessInitialDensityR(float h, REAL m, REAL &rd) 
    : kern(h), mass(m), rest_density(rd) { }
  ~ProcessInitialDensityR() { }

  inline void init(ParticleR<REAL> &p)
  {
    qDebug() << p.pos[0] << p.pos[1] << p.pos[2];
    p.dinv = 0.0f;
  }

  inline void operator()(ParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    p.dinv += kern[ p.pos - near_p.pos ];
  }

  inline void finish(ParticleR<REAL> &p)
  {
    REAL density = p.dinv * mass * kern.coef;
    p.dinv = 1.0f/density;
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
  ProcessInitialDensityR<REAL> proc(m_h, m_dpc->m_mass, m_dpc->m_rest_density);
  compute_bound_quantity(proc);

  m_dpc->m_rest_density = 0.0f;
  compute_fluid_quantity(proc);
  m_dpc->m_rest_density /= m_dpc->m_num_vertices;
  glprintf_tr("rest density: %.2f\n", m_dpc->m_rest_density);
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::compute_pressure()
{
  compute_bound_quantity(m_proc_pressure);
  compute_fluid_quantity(m_proc_pressure);
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::compute_accel()
{
  m_dpc->reset_accel();     // now may assume all accelerations are zero
  compute_fluid_quantity< ProcessAccelR<REAL> >( m_proc_accel );
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
          process.init(p);

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
                for ( DynamicParticleR<REAL> &near_p : neigh_fluidvec )
                {
                  process(p, near_p); // process data
                }
                for ( ParticleR<REAL> &near_p : neigh_boundvec )
                {
                  process(p, near_p); // process data
                }
              }
            }
          }
        }

        for ( ParticleR<REAL> &p : boundvec )  // finalize data
          process.finish(p);
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

        for ( DynamicParticleR<REAL> &p : fluidvec )  // prepare data
          process.init(p);

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
              for ( DynamicParticleR<REAL> &p : fluidvec )
              {
                for ( DynamicParticleR<REAL> &near_p : neigh_fluidvec )
                {
                  process(p, near_p); // process data
                }
                for ( ParticleR<REAL> &near_p : neigh_boundvec )
                {
                  process(p, near_p); // process data
                }
              }
            }
          }
        }

        for ( DynamicParticleR<REAL> &p : fluidvec )  // finalize data
          process.finish(p);
      } // for k
    } // for j
  } // for i
}

// DynamicPointCloud stuff

// inflating the size of the grid
#define INFLATE 3.0 

template<typename REAL, typename SIZE>
DynamicPointCloudRS<REAL,SIZE>::DynamicPointCloudRS(
    GLPointCloudRS<REAL, SIZE> *glpc,
    REAL density, REAL viscosity, REAL st)
  : PointCloudRS<REAL,SIZE>(*glpc->get_pointcloud())
  , m_num_vertices(this->get_num_vertices())
  , m_c2(2.2e2)
  , m_rest_density(density)
  , m_viscosity(viscosity)
  , m_st(st) // surface tension
  , m_bmin(this->m_bbox.corner(Eigen::AlignedBox3f::BottomLeftFloor))
  , m_bmax(this->m_bbox.corner(Eigen::AlignedBox3f::TopRightCeil))
  , m_kernel_radius(INFLATE*this->get_radius())
  , m_grid(this, m_kernel_radius, m_bmin)
  , m_glpc(glpc)
  , m_stop_requested(false)
{
  m_accel.resizeLike(this->m_pos);
  m_vel.resizeLike(this->m_pos);
  m_vel.setZero();
  reset_accel();

  REAL h = this->compute_mindist();
  m_mass = m_rest_density * h * h * h * 1.1f;

  //const clock_t begin_time = clock();
  //float d1 = PointCloudRS<REAL,SIZE>::compute_mindist();
  //float t1 = clock();
  //float d2 = PointCloudRS<REAL,SIZE>::compute_mindist_brute();
  //float t2 = clock();
  //qDebug() << float( t1 - begin_time ) / CLOCKS_PER_SEC << " to get " << d1;
  //qDebug() << float( t2 - t1 ) / CLOCKS_PER_SEC << " to get " << d2;
}

template<typename REAL, typename SIZE>
DynamicPointCloudRS<REAL,SIZE>::~DynamicPointCloudRS()
{
}


// clamp value d to min and max boundaries,
// return -1 if clamed with max, 1 if clamped with min, otherwise 0
template<typename REAL, typename SIZE>
inline bool DynamicPointCloudRS<REAL,SIZE>::clamp(REAL &d, REAL min, REAL max)
{
  if ( d < min )
  {
    d = min;
    return true;
  }
  else if (d > max)
  {
    d = max;
    return true;
  }
  return false;
}


template<typename REAL, typename SIZE>
void DynamicPointCloudRS<REAL,SIZE>::resolve_collisions()
{
  for (SIZE i = 0; i < m_num_vertices; ++i) // TODO: vectorize this
  {
    for (unsigned char j = 0; j < 3; ++j)
    {
      if (clamp(pos_at(i)[j], m_bmin[j], m_bmax[j]))
      {
        vel_at(i)[j] *= -1.0;
      }
    }
  }
}

template<typename REAL, typename SIZE>
void DynamicPointCloudRS<REAL,SIZE>::run()
{
  m_grid.init();
  float dt = 0.00217;
  glprintf_tr("step: %.2es\n", dt);

  m_grid.compute_pressure();
  m_grid.compute_accel(); // update m_accel

  m_vel = m_vel + 0.5*dt*m_accel; // initial half velocity

  clock_t prev_t = clock();
  float t = 0.0f;
  unsigned int count=0;
  for ( ; ; )
  {
    this->m_pos = (this->m_pos + dt*m_vel).eval();
    //resolve_collisions();
    m_grid.update();

    if (m_glpc)
      m_glpc->update_data(); // notify gl we have new positions

    if (m_stop_requested)
      break;

    m_vel = m_vel + 0.5*dt*m_accel;

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

template class DynamicPointCloudRS<double, unsigned int>;
template class UniformGridRS<double, unsigned int>;
