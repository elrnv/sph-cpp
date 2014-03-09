#include "dynamics.h"
#include "glpointcloud.h"

// UniformGrid stuff

template<typename REAL, typename SIZE>
UniformGridRS<REAL,SIZE>::UniformGridRS(DynamicPointCloudRS<REAL,SIZE> *dpc, float h, const Vector3f &bmin)
  : m_dpc(dpc)
  , m_h(h)
  , m_hinv(1.0f/h)
  , m_proc_pressure(m_h)
  , m_proc_accel(m_h)
  , m_bmin(bmin)
{ }

template<typename REAL, typename SIZE>
UniformGridRS<REAL,SIZE>::~UniformGridRS() 
{ }
  
template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::init()
{
  // determine the number of voxels needed
  const Vector3f &bmax = m_dpc->get_bmax();

  SIZE nx = 1 + (bmax[0] - m_bmin[0])*m_hinv;
  SIZE ny = 1 + (bmax[1] - m_bmin[1])*m_hinv;
  SIZE nz = 1 + (bmax[2] - m_bmin[2])*m_hinv;
  
  m_grid.resize( boost::extents[nx][ny][nz] );

  SIZE num_vtx = m_dpc->get_num_vertices();
  for ( SIZE i = 0; i < num_vtx; ++i )
  {
    Vector3R<REAL> pos( m_dpc->pos_at(i) );
    Vector3R<REAL> vel( m_dpc->vel_at(i) );
    Index idx = get_voxel_index(pos);
    if (is_valid(idx))
      m_grid(idx).push_back( ParticleDataR<REAL>(pos, vel, m_dpc->accel_at(i)) );
    else
      qDebug() << "particle off grid at creation";
  }

  compute_initial_density();

  m_proc_pressure.init(m_dpc->m_mass, m_dpc->m_rest_density, m_dpc->m_constant);
  m_proc_accel.init(m_dpc->get_mass());
}

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::update()
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
        
        typename DataVec::iterator it = datavec.begin();
        for ( ; it != datavec.end(); )
        {
          Index idx(get_voxel_index(it->pos));
          if (!is_valid(idx))
          {
            qDebug() << "particle off grid";
            datavec.erase(it);
            continue;
          }

          if ( idx[0] != i || idx[1] != j || idx[2] != k )
          { 
            datavec.erase(it);
            m_grid(idx).push_back(*it);
            continue;
          }
          ++it;
        }
      }
    }
  }
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
                  if (&p != &near_p)
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

template<typename REAL, typename SIZE>
DynamicPointCloudRS<REAL,SIZE>::DynamicPointCloudRS(
    GLPointCloudRS<REAL, SIZE> *glpc,
    REAL mass)
  : PointCloudRS<REAL,SIZE>(*glpc->get_pointcloud())
  , m_num_vertices(this->get_num_vertices())
  , m_mass(mass)
  , m_constant(300.0f * 0.3f)
  , m_rest_density(0.0f)
  , m_bmin(this->m_bbox.corner(Eigen::AlignedBox3f::BottomLeftFloor))
  , m_bmax(this->m_bbox.corner(Eigen::AlignedBox3f::TopRightCeil))
  , m_grid(this, 0.04f, m_bmin)
  , m_glpc(glpc)
  , m_stop_requested(false)
{
  m_accel.resizeLike(this->m_pos);
  m_vel.resizeLike(this->m_pos);
  m_vel.setZero();
  reset_accel();
}

template<typename REAL, typename SIZE>
DynamicPointCloudRS<REAL,SIZE>::~DynamicPointCloudRS()
{
}


// clamp value d to min and max boundaries,
// return -1 if clamed with max, 1 if clamped with min, otherwise 0
template<typename REAL, typename SIZE>
inline int DynamicPointCloudRS<REAL,SIZE>::clamp(REAL &d, REAL min, REAL max)
{
  if ( d < min )
  {
    d = min;
    return 1;
  }
  else if (d > max)
  {
    d = max;
    return 1;
  }
  return 0;
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
        vel_at(i)[j] *= -1;
      }
    }
  }
}

template<typename REAL, typename SIZE>
void DynamicPointCloudRS<REAL,SIZE>::run()
{
  m_grid.init();
  float dt = 1e-8;

  m_grid.compute_pressure();
  m_grid.compute_accel(); // update m_accel

  m_vel = m_vel + 0.5*dt*m_accel; // initial half velocity

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
  }
}

template class DynamicPointCloudRS<double, unsigned int>;
template class UniformGridRS<double, unsigned int>;
