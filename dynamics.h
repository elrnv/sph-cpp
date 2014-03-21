#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <vector>
#include <boost/multi_array.hpp>
#include "kernel.h"
#include "pointcloud.h"

#define M_G 9.81f

// Partially matrix template for convenience
template<typename REAL>
using Vector3R = Matrix<REAL, 3, 1>;

// forward declarations
template<typename REAL, typename SIZE>
class DynamicPointCloudRS;

template<typename REAL, typename SIZE>
class GLPointCloudRS;

template<typename REAL>
struct ParticleR
{
  ParticleR(Vector3R<REAL> p) : pos(p){ }
  ~ParticleR() { }

  Vector3R<REAL> pos;
  REAL dinv;
  REAL pressure;
};

template<typename REAL>
struct DynamicParticleR : public ParticleR<REAL>
{
  DynamicParticleR(Vector3R<REAL> p, Vector3R<REAL> v, REAL *a, REAL *ea)
    : ParticleR<REAL>(p), vel(v), accel(a), extern_accel(ea) { }
  ~DynamicParticleR() { }

  Vector3R<REAL> vel;
  REAL *accel;    // 3 array to which we will write total acceleration
  REAL *extern_accel;    // 3 array to which we will write acceleration due to external forces
  REAL temp;      // store temporary values at the particle during computation
  REAL vol;       // volume estimate
};


// From [Solenthaler and Pajarola 2008], an alternative density
// (number_density * mass) is used
template<typename REAL>
class ComputeBoundaryVolumeR
{
public:
  ComputeBoundaryVolumeR(float h) : kern(h) { }
  ~ComputeBoundaryVolumeR() { }

  inline void init_particle(ParticleR<REAL> &p) { p.dinv = 0.0f; }

  inline void fluid(ParticleR<REAL> &p, DynamicParticleR<REAL> &near_p)
  { }

  inline void bound(ParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    p.dinv += kern[ p.pos - near_p.pos ];
  }

  inline void finish_particle(ParticleR<REAL> &p)
  {
    p.dinv = 1.0f/(p.dinv * kern.coef);
  }

private:
  Poly6Kernel kern; // used to compute pressure force
};

template<typename REAL>
class ComputeFluidDensityR
{
public:
  ComputeFluidDensityR(float h) : kern(h) { }
  ~ComputeFluidDensityR() { }

  void init(REAL m, REAL rd, REAL r, REAL &mv, REAL &av)
  {
    mass = m; rest_density = rd; radius = r; 
    max_var = &mv; avg_var = &av;
  }

  inline void init_particle(DynamicParticleR<REAL> &p)
  { 
    p.dinv = 0.0f; p.vol = 0.0f;
  }
  inline void fluid(DynamicParticleR<REAL> &p, DynamicParticleR<REAL> &near_p)
  {
    p.dinv += kern[ p.pos - near_p.pos ];
    p.vol += kern[ p.pos - near_p.pos ];
  }
  inline void bound(DynamicParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    p.dinv += rest_density * near_p.dinv * kern[ p.pos - near_p.pos ];
    p.vol += kern[ p.pos - near_p.pos ];
  }
  inline void finish_particle(DynamicParticleR<REAL> &p)
  {
    //qDebug() << (mass * p.dinv) << "/" <<(8*radius*radius*radius*p.vol) << " = " << (mass * p.dinv)/(8*radius*radius*radius*p.vol);
    p.dinv = (8*radius*radius*radius*p.vol)/(mass * p.dinv);
    REAL var = std::abs(1.0f/p.dinv - rest_density);
    if ( var > *max_var )
      *max_var = var;
    *avg_var += var;
  }

private:
  Poly6Kernel kern; // used to compute pressure force
  REAL mass;
  REAL rest_density;
  REAL radius;
  REAL *max_var; // max variation
  REAL *avg_var; // average variation
};

template<typename REAL>
class ComputeFluidDensityUpdateR
{
public:
  ComputeFluidDensityUpdateR(float h) : kern(h) { }
  ~ComputeFluidDensityUpdateR() { }

  void init(REAL m, float ts) 
  { mass = m; timestep = ts; }

  inline void init_particle(DynamicParticleR<REAL> &p) { p.temp = 0.0f; }

  inline void fluid(DynamicParticleR<REAL> &p, DynamicParticleR<REAL> &near_p)
  {
    p.temp += (near_p.vel - p.vel).dot(kern[ p.pos - near_p.pos ]);
  }
  inline void bound(DynamicParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    p.temp += (-p.vel).dot(kern[ p.pos - near_p.pos ]);
    //p.dinv += rest_density * near_p.dinv * kern[ p.pos - near_p.pos ];
  }

  inline void finish_particle(DynamicParticleR<REAL> &p)
  {
    p.dinv += timestep * p.dinv * p.dinv * p.temp * kern.coef;
  }

private:
  Poly6GradKernel kern; // used to compute pressure force
  REAL mass;
  float timestep;
};

// process routines used to compute SPH quantities
template<typename REAL>
class ComputePressureR
{
public:
  ComputePressureR(float h) : kern(h) { }
  ~ComputePressureR() { }
  void init(REAL m, REAL rd, REAL c2) 
  { mass = m; rest_density = rd; cs2 = c2; }

  inline void init_particle(ParticleR<REAL> &p) { p.dinv = 0.0f; }

  inline void fluid(ParticleR<REAL> &p, DynamicParticleR<REAL> &near_p)
  {
    p.dinv += kern[ p.pos - near_p.pos ];
  }

  inline void bound(ParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    p.dinv += rest_density * near_p.dinv * kern[ p.pos - near_p.pos ];
  }

  inline REAL pow7(REAL x) { return x*x*x*x*x*x*x; }

  inline void finish_particle(ParticleR<REAL> &p)
  {
    REAL density = mass * p.dinv * kern.coef;
    p.dinv = 1.0f/density;
    p.pressure = rest_density * cs2 * 0.14285714285714 * (pow7(density / rest_density) - 1);
    //qDebug() << p.pressure;
  }

private:
  Poly6Kernel kern; // used to compute pressure force
  REAL mass;
  REAL rest_density;
  REAL cs2;
};


template<typename REAL>
class ComputePressureAccelR
{
public:
  ComputePressureAccelR(float h) : b_kern(h), p_kern(h), v_kern(h) { }
  ~ComputePressureAccelR() { }

  void init(REAL m, REAL rd, REAL v, REAL s, REAL c2)
  {  mass = m; rest_density = rd; viscosity = v; st = s; const2 = c2; }

  inline void init_particle(DynamicParticleR<REAL> &p) {
    Q_UNUSED(p); 
  }

  inline void fluid(DynamicParticleR<REAL> &p, DynamicParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;
    Vector3R<REAL> res(
        -mass*(p.pressure*p.dinv*p.dinv + near_p.pressure*near_p.dinv*near_p.dinv)*p_kern(p.pos - near_p.pos) // pressure contribution
       // + viscosity*near_p.dinv*(near_p.vel - p.vel)*v_kern(p.pos - near_p.pos) // viscosity contribution
        );

    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result
  }

  inline void bound(DynamicParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
   // Vector3R<REAL> res(
   //     -rest_density*near_p.dinv*p.pressure*p.dinv*p.dinv*p_kern[p.pos - near_p.pos] // pressure contribution
   //    // + viscosity*near_p.dinv*(near_p.vel - p.vel)*v_kern(p.pos - near_p.pos) // viscosity contribution
   //     );
    float massb = rest_density * near_p.dinv;
    Vector3R<REAL> res(
        const2*(massb/(mass*(mass + massb))) * (p.pos - near_p.pos) * b_kern(p.pos - near_p.pos)  // pressure contribution
       // + viscosity*near_p.dinv*(near_p.vel - p.vel)*v_kern(p.pos - near_p.pos) // viscosity contribution
        );

    for (unsigned char i = 0; i < 3; ++i)
      p.extern_accel[i] += res[i]; // copy intermediate result
  }

  inline void finish_particle(DynamicParticleR<REAL> &p)
  {
    p.extern_accel[1] -= M_G;
    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] = p.accel[i] + p.extern_accel[i];

//    qDebug() << "p.dinv:" << p.dinv;
    //qDebug(" ac: % 10.5e % 10.5e % 10.5e", p.accel[0], p.accel[1], p.accel[2]);

  }

private:
  LeonardJonesKernel b_kern; // used to compute boundary force
  SpikyGradKernel p_kern; // used to compute pressure force
  ViscLapKernel   v_kern; // used to compute viscosity force
  REAL mass;
  REAL rest_density;
  REAL viscosity;
  REAL st;
  REAL const2;
};


// Grid structure used to optimize computing particle properties using kernels
template<typename REAL, typename SIZE>
class UniformGridRS
{
public:
  typedef std::vector< DynamicParticleR<REAL> > DynamicParticles;
  typedef std::vector< ParticleR<REAL> > StaticParticles;

  struct Cell
  {
    DynamicParticles fluidvec; // fluid particles
    StaticParticles boundvec; // static boundary particles
    //DynamicParticles rigidvec; // dynamic rigid body object data
  };

  typedef boost::multi_array< Cell, 3 > Array3;
  typedef typename Array3::index Index;
  typedef boost::array<Index, 3> Array3Index;
  typedef typename Array3::template array_view<3>::type GridView;
  typedef typename Array3::index_range IndexRange;

  UniformGridRS(DynamicPointCloudRS<REAL,SIZE> *dpc, float h, float grid_h, const Vector3f &bmin);
  UniformGridRS(DynamicPointCloudRS<REAL,SIZE> *dpc, float h, const Vector3f &bmin);
  ~UniformGridRS();
  
  void init();
  inline void update();
  inline void populate_fluid_data();
  inline void populate_bound_data();
  inline void clear_fluid_data();

  inline Index clamp(Index d, Index min, Index max)
  {
    return std::max(std::min(d, max), min);
  }

  inline Array3Index get_voxel_index(const Vector3R<REAL> &pos)
  {
    return {{ 
      clamp(static_cast<Index>(m_hinv*(pos[0]-m_bmin[0])), 0, m_gridsize[0]-1),
      clamp(static_cast<Index>(m_hinv*(pos[1]-m_bmin[1])), 0, m_gridsize[1]-1),
      clamp(static_cast<Index>(m_hinv*(pos[2]-m_bmin[2])), 0, m_gridsize[2]-1) }};
  }

  void compute_accel();
  void compute_initial_density();
  void compute_pressure();
  void compute_density();
  void update_density(float timestep);

  template<typename ProcessPairFunc>
  inline void compute_bound_quantity(ProcessPairFunc process);

  template<typename ProcessPairFunc>
  inline void compute_fluid_quantity(ProcessPairFunc process);

  inline void resolve_boundary_collisions();
  inline void process_collisions_in_cell(SIZE i, SIZE j, SIZE k);

private:
  // utility function used in the constructor to get a range of two elements
  // centered at x (or 2 if x is on the boundary)
  inline IndexRange range3(SIZE x, SIZE hi)
  {
    return IndexRange(x == 0 ? 0 : x-1, x == hi-1 ? hi : x+2 );
  }

private:
  DynamicPointCloudRS<REAL, SIZE> *m_dpc; // main point cloud

  // array of cells containing xyzp (position and density) for each vertex
  Array3 m_grid;
  Array3Index  m_gridsize;

  float m_h;    // kernel size
  float m_hinv; // 1 / grid_size

  ComputeFluidDensityR<REAL>   m_proc_fluid_density;
  ComputeBoundaryVolumeR<REAL> m_proc_bound_volume;
  ComputePressureR<REAL>       m_proc_pressure;
  ComputePressureAccelR<REAL>  m_proc_pressure_accel;

  const Vector3f &m_bmin;
}; // class UniformGridRS


// A dynamic cloud of points
template<typename REAL, typename SIZE>
class DynamicPointCloudRS : public PointCloudRS<REAL,SIZE>
{
public:
  // dynamic point cloud from a regular updatable gl point cloud
  explicit DynamicPointCloudRS(GLPointCloudRS<REAL, SIZE> *glpc,
      REAL density, REAL viscosity, REAL surfacetension);
  ~DynamicPointCloudRS();

  inline REAL *pos_at(SIZE i)    { return this->m_pos.data() + i*3; }
  inline REAL *vel_at(SIZE i)    { return m_vel.data() + i*3; }
  inline REAL *accel_at(SIZE i)  { return m_accel.data() + i*3; }
  inline REAL *extern_accel_at(SIZE i)  { return m_extern_accel.data() + i*3; }

  inline const Vector3f &get_bmin() const { return m_bmin; }
  inline const Vector3f &get_bmax() const { return m_bmax; }

  inline REAL get_kernel_radius() const { return m_kernel_radius; }

  inline bool is_dynamic() const { return true; }

  inline void reset_accel() { m_accel.setZero(); m_extern_accel.setZero(); }

  inline bool clamp(REAL &d, REAL min, REAL max);
  inline void resolve_collisions();

  void run();
  void request_stop() { m_stop_requested = true; }

  friend UniformGridRS<REAL,SIZE>;

protected:
  SIZE m_num_vertices;
  REAL m_c2;
  REAL m_mass;
  REAL m_rest_density;
  REAL m_viscosity;
  REAL m_st;
  Matrix3XR<REAL> m_vel; // velocities
  Matrix3XR<REAL> m_accel; // accelerations
  Matrix3XR<REAL> m_extern_accel; // accelerations

  Vector3f m_bmin;
  Vector3f m_bmax;

  REAL m_kernel_radius;

  UniformGridRS<REAL, SIZE> m_grid;

  GLPointCloudRS<REAL,SIZE> *m_glpc; // should only used for callback

  std::atomic<bool> m_stop_requested;
}; // class DynamicPointCloudRS

// defaults
typedef DynamicPointCloudRS<double, unsigned int> DynamicPointCloud;
typedef UniformGridRS<double, unsigned int> UniformGrid;

#endif // DYNAMICS_H
