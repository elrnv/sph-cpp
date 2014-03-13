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
struct ParticleDataR
{
  ParticleDataR(Vector3R<REAL> p, Vector3R<REAL> v, REAL *a)
    : pos(p), vel(v), accel(a){ }
  ~ParticleDataR() { }

  Vector3R<REAL> pos;
  Vector3R<REAL> vel;
  REAL *accel; // 3 array to which we will write
  REAL dinv;
  REAL pressure;
};


// process routines used to compute SPH quantities
template<typename REAL>
class ProcessPressureR
{
public:
  ProcessPressureR(float h) : kern(h) { }
  ~ProcessPressureR() { }
  void init(REAL m, REAL rd, REAL c) 
  { mass = m; rest_density = rd; constant = c; }

  inline void init(ParticleDataR<REAL> &p)
  {
    p.dinv = 0.0f;
  }

  inline void operator()(ParticleDataR<REAL> &p, ParticleDataR<REAL> &near_p)
  {
    p.dinv += kern[ p.pos - near_p.pos ];
    //qDebug() << "temp p.dinv" << p.dinv;
  }

  inline REAL pow7(REAL x) { return x*x*x*x*x*x*x; }

  inline void finish(ParticleDataR<REAL> &p)
  {
    REAL density = p.dinv * mass * kern.coef;
    p.dinv = 1.0f/density;
    p.pressure = constant * (pow7(density / rest_density) - 1);
  }

private:
  Poly6Kernel kern; // used to compute pressure force
  REAL mass;
  REAL rest_density;
  REAL constant;
};


template<typename REAL>
class ProcessAccelR
{
public:
  ProcessAccelR(float h) : p_kern(h), v_kern(h) { }
  ~ProcessAccelR() { }

  void init(REAL m, REAL v, REAL s)
  {  mass = m; viscosity = v; st = s; }

  inline void init(ParticleDataR<REAL> &p) {
    Q_UNUSED(p); 
  }

  inline void operator()(ParticleDataR<REAL> &p, ParticleDataR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;
    Vector3R<REAL> res(
        -0.5*near_p.dinv*(p.pressure + near_p.pressure)*p_kern(p.pos - near_p.pos) // pressure contribution
//        + viscosity*near_p.dinv*(near_p.vel - p.vel)*v_kern(p.pos - near_p.pos) // viscosity contribution
        );

    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result
  }

  inline void finish(ParticleDataR<REAL> &p)
  {
    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] = mass*p.dinv*p.accel[i];
    p.accel[1] -= M_G;

//    qDebug() << "p.dinv:" << p.dinv;
    //qDebug(" ac: % 10.5e % 10.5e % 10.5e", p.accel[0], p.accel[1], p.accel[2]);

  }

private:
  SpikyGradKernel p_kern; // used to compute pressure force
  ViscLapKernel   v_kern; // used to compute viscosity force
  REAL mass;
  REAL viscosity;
  REAL st;
};


// Grid structure used to optimize computing particle properties using kernels
template<typename REAL, typename SIZE>
class UniformGridRS
{
public:
  typedef std::vector< ParticleDataR<REAL> > DataVec;

  typedef boost::multi_array< DataVec, 3 > Array3;
  typedef typename Array3::index Index;
  typedef boost::array<Index, 3> Array3Index;
  typedef typename Array3::template array_view<3>::type GridView;
  typedef typename Array3::index_range IndexRange;

  UniformGridRS(DynamicPointCloudRS<REAL,SIZE> *dpc, float h, float grid_h, const Vector3f &bmin);
  UniformGridRS(DynamicPointCloudRS<REAL,SIZE> *dpc, float h, const Vector3f &bmin);
  ~UniformGridRS();
  
  void init();

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

  template<typename ProcessFunc>
  inline void compute_quantity(ProcessFunc process);

  inline void update();

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

  ProcessPressureR<REAL> m_proc_pressure;
  ProcessAccelR<REAL> m_proc_accel;

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

  inline const Vector3f &get_bmin() const { return m_bmin; }
  inline const Vector3f &get_bmax() const { return m_bmax; }

  inline REAL get_kernel_radius() const { return m_kernel_radius; }

  inline bool is_dynamic() const { return true; }

  inline void reset_accel() { m_accel.setZero(); }

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
