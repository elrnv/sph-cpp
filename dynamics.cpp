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
  glprintf_tr("cell size: %.2f\n", grid_h);
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
  glprintf_tr("cell size: %.2f\n", h);
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

  qDebug() << m_gridsize[0] << m_gridsize[1] << m_gridsize[2];
  
  m_grid.resize( boost::extents[m_gridsize[0]][m_gridsize[1]][m_gridsize[2]] );

  SIZE num_vtx = m_dpc->get_num_vertices();
  for ( SIZE i = 0; i < num_vtx; ++i )
  {
    Vector3R<REAL> pos( m_dpc->m_pos.col(i) );
    Vector3R<REAL> vel( m_dpc->m_vel.col(i) );
    Array3Index idx = get_voxel_index(pos);
    m_grid(idx).push_back( ParticleDataR<REAL>(pos, vel, m_dpc->accel_at(i)) );
  }

  compute_initial_density();

  qDebug() << "mass per particle" << m_dpc->m_mass;

  m_proc_pressure.init(m_dpc->m_mass, m_dpc->m_rest_density, m_dpc->m_c2);
  m_proc_accel.init(m_dpc->m_mass, m_dpc->m_viscosity, m_dpc->m_st);
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::update()
{
#if 1
  int size = 0;
  int count = 0;
#endif
  SIZE nx = m_gridsize[0];
  SIZE ny = m_gridsize[1];
  SIZE nz = m_gridsize[2];
  for (SIZE i = 0; i < nx; ++i)
  {
    for (SIZE j = 0; j < ny; ++j)
    {
      for (SIZE k = 0; k < nz; ++k)
      {
        DataVec &datavec = m_grid[i][j][k];
        datavec.clear();
      } // for k
    } // for j
  } // for i


  SIZE num_vtx = m_dpc->get_num_vertices();
  for ( SIZE i = 0; i < num_vtx; ++i )
  {
    Vector3R<REAL> pos( m_dpc->m_pos.col(i) );
    Vector3R<REAL> vel( m_dpc->m_vel.col(i) );
    Array3Index idx = get_voxel_index(pos);
    m_grid(idx).push_back( ParticleDataR<REAL>(pos, vel, m_dpc->accel_at(i)) );
  }
  //qDebug() << "avg particles per cell:" << size / count;
}

template<typename REAL>
class ProcessInitialDensityR
{
public:
  ProcessInitialDensityR(float h, REAL m, REAL &rd) 
    : kern(h), mass(m), rest_density(rd) { }
  ~ProcessInitialDensityR() { }

  inline void init(ParticleDataR<REAL> &p)
  {
    p.dinv = 0.0f;
  }

  inline void operator()(ParticleDataR<REAL> &p, ParticleDataR<REAL> &near_p)
  {
    p.dinv += kern[ p.pos - near_p.pos ];
  }

  inline void finish(ParticleDataR<REAL> &p)
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
  m_dpc->m_rest_density = 0.0f;
  ProcessInitialDensityR<REAL> proc(m_h, m_dpc->m_mass, m_dpc->m_rest_density);
  compute_quantity(proc);
  m_dpc->m_rest_density /= m_dpc->m_num_vertices;
  qDebug() << "rest density" << m_dpc->m_rest_density;
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::compute_pressure()
{
  compute_quantity(m_proc_pressure);
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::compute_accel()
{
  m_dpc->reset_accel();     // now may assume all accelerations are zero
  compute_quantity< ProcessAccelR<REAL> >( m_proc_accel );
}

template<typename REAL, typename SIZE>
template<typename ProcessFunc>
void UniformGridRS<REAL,SIZE>::compute_quantity(ProcessFunc process)
{
  SIZE nx = m_grid.shape()[0];
  SIZE ny = m_grid.shape()[1];
  SIZE nz = m_grid.shape()[2];
  for (SIZE i = 0; i < nx; ++i)
  {
    for (SIZE j = 0; j < ny; ++j)
    {
      for (SIZE k = 0; k < nz; ++k)
      {
        DataVec &datavec = m_grid[i][j][k];
        if (datavec.empty())
          continue;

        IndexRange xrange = range3(i,nx);
        IndexRange yrange = range3(j,ny);
        IndexRange zrange = range3(k,nz);
        GridView neigh_view = m_grid[ boost::indices[xrange][yrange][zrange] ];

        for ( ParticleDataR<REAL> &p : datavec )  // prepare data
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
              for ( ParticleDataR<REAL> &p : datavec )
              {
                for ( ParticleDataR<REAL> &near_p : neigh_view[near_i][near_j][near_k])
                {
                  process(p, near_p); // process data
                }
              }
            }
          }
        }

        for ( ParticleDataR<REAL> &p : datavec )  // finalize data
          process.finish(p);
      } // for k
    } // for j
  } // for i
}

// DynamicPointCloud stuff

// inflating the size of the grid
#define INFLATE 6.0 

template<typename REAL, typename SIZE>
DynamicPointCloudRS<REAL,SIZE>::DynamicPointCloudRS(
    GLPointCloudRS<REAL, SIZE> *glpc,
    REAL density, REAL viscosity, REAL st)
  : PointCloudRS<REAL,SIZE>(*glpc->get_pointcloud())
  , m_num_vertices(this->get_num_vertices())
  , m_c2(2.2e4)
  , m_rest_density(density)
  , m_viscosity(viscosity)
  , m_st(st) // surface tension
  , m_bmin(this->m_bbox.corner(Eigen::AlignedBox3f::BottomLeftFloor))
  , m_bmax(this->m_bbox.corner(Eigen::AlignedBox3f::TopRightCeil))
  , m_grid(this, INFLATE*this->compute_mindist(), 2.0f, m_bmin)
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
        vel_at(i)[j] *= -0.1;
      }
    }
  }
}

template<typename REAL, typename SIZE>
void DynamicPointCloudRS<REAL,SIZE>::run()
{
  m_grid.init();
  float dt = 0.000217;

  m_grid.compute_pressure();
  m_grid.compute_accel(); // update m_accel

  m_vel = m_vel + 0.5*dt*m_accel; // initial half velocity

  clock_t prev_t = clock();
  float t = 0.0f;
  unsigned int count=0;
  for ( ; ; )
  {
    this->m_pos = (this->m_pos + dt*m_vel).eval();
    resolve_collisions();
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
  qDebug() << "average time per step:" << t / float(count);
}

template class DynamicPointCloudRS<double, unsigned int>;
template class UniformGridRS<double, unsigned int>;
