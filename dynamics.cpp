#include "dynamics.h"

// UniformGrid stuff

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


template<typename REAL>
class ProcessPressureR
{
public:
  ProcessPressureR(float h, REAL m, REAL rd, REAL c) 
    : kern(h), mass(m), rest_density(rd), constant(c) { }
  ~ProcessPressureR() { }

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
    p.pressure = constant * (density - rest_density);
  }

private:
  Poly6Kernel kern; // used to compute pressure force
  REAL mass;
  REAL rest_density;
  REAL constant;
};

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::compute_pressure()
{
  ProcessPressureR<REAL> proc(m_h,
      m_dpc->m_mass, m_dpc->m_rest_density, m_dpc->m_constant);
  compute_quantity(proc);
}


template<typename REAL>
class ProcessAccelR
{
public:
  ProcessAccelR(float h, REAL m)
    : kern(h), mass(m) { }
  ~ProcessAccelR() { }

  inline void init(ParticleDataR<REAL> &p) { Q_UNUSED(p); }

  inline void operator()(ParticleDataR<REAL> &p, ParticleDataR<REAL> &near_p)
  {
    Vector3R<REAL> res(0.5*near_p.dinv*(p.pressure + near_p.pressure)*kern[p.pos - near_p.pos]);
    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result
  }

  inline void finish(ParticleDataR<REAL> &p)
  {
    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] = -(mass*p.dinv*kern.coef*p.accel[i]);
  }

private:
  SpikyGradKernel kern; // used to compute pressure force
  REAL mass;
};

template<typename REAL, typename SIZE>
void UniformGridRS<REAL,SIZE>::compute_accel()
{
  m_dpc->reset_accel();     // now may assume all accelerations are zero
  ProcessAccelR<REAL> proc(m_h, m_dpc->get_mass());
  compute_quantity< ProcessAccelR<REAL> >( proc );
}

template<typename REAL, typename SIZE>
template<typename ProcessFunc>
void UniformGridRS<REAL,SIZE>::compute_quantity(ProcessFunc process)
{
  SIZE nx = m_data.shape()[0];
  SIZE ny = m_data.shape()[1];
  SIZE nz = m_data.shape()[2];
  for (SIZE i = 0; i < nx; ++i)
  {
    for (SIZE j = 0; j < ny; ++j)
    {
      for (SIZE k = 0; k < nz; ++k)
      {
        IndexRange xrange = range3(i,nx);
        IndexRange yrange = range3(j,ny);
        IndexRange zrange = range3(k,nz);
        GridView neigh_view = m_data[ boost::indices[xrange][yrange][zrange] ];
        DataVec &datavec = m_data[i][j][k];

        for ( ParticleDataR<REAL> &p : datavec )  // prepare data
          process.init(p);

        SIZE xrange_size = xrange.finish() - xrange.start();
        SIZE yrange_size = yrange.finish() - yrange.start();
        SIZE zrange_size = zrange.finish() - zrange.start();
        for (SIZE near_i = 0; near_i < xrange_size; ++i)
          for (SIZE near_j = 0; near_j < yrange_size; ++j)
            for (SIZE near_k = 0; near_k < zrange_size; ++k)
              for ( ParticleDataR<REAL> &p : datavec )
                for ( ParticleDataR<REAL> &near_p : neigh_view[near_i][near_j][near_k])
                  process(p, near_p); // process data

        for ( ParticleDataR<REAL> &p : datavec )  // finalize data
          process.finish(p);
      } // for k
    } // for j
  } // for i
}

// DynamicPointCloud stuff

template<typename REAL, typename SIZE>
DynamicPointCloudRS<REAL,SIZE>::DynamicPointCloudRS(
    const aiMesh *mesh,
    REAL mass,
    void (* update_data_callback)(void))
  : PointCloudRS<REAL,SIZE>(mesh)
  , m_num_vertices(this->get_num_vertices())
  , m_mass(mass)
  , m_constant(300.0f * 0.8f)
  , m_rest_density(0.0f)
  , m_grid(this, 0.1f)
  , m_update_data_callback(update_data_callback)
{
  m_accel.resizeLike(this->m_pos);
  m_vel.resizeLike(this->m_pos);
  m_vel.setZero();
  reset_accel();
}

template<typename REAL, typename SIZE>
DynamicPointCloudRS<REAL,SIZE>::DynamicPointCloudRS(
    const PointCloudRS<REAL, SIZE> &pc,
    REAL mass,
    void (* update_data_callback)(void))
  : PointCloudRS<REAL,SIZE>(pc)
  , m_num_vertices(this->get_num_vertices())
  , m_mass(mass)
  , m_constant(300.0f * 0.8f)
  , m_rest_density(0.0f)
  , m_grid(this, 0.1f)
  , m_update_data_callback(update_data_callback)
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

template class DynamicPointCloudRS<double, unsigned int>;
template class UniformGridRS<double, unsigned int>;
