#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <vector>
#include <boost/multi_array.hpp>
#include "kernel.h"
#include "pointcloud.h"

// Partially matrix template for convenience
template<typename REAL>
using Vector3R = Matrix<REAL, 3, 1>;

// forward declaration
template<typename REAL, typename SIZE>
class DynamicPointCloudRS;


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

// Grid structure used to optimize computing particle properties using kernels
template<typename REAL, typename SIZE>
class UniformGridRS
{
public:
  typedef std::vector< ParticleDataR<REAL> > DataVec;

  typedef boost::multi_array< DataVec, 3 > Grid;
  typedef boost::array<typename Grid::index, 3> Index;
  typedef typename Grid::template array_view<3>::type GridView;
  typedef typename Grid::index_range IndexRange;

  UniformGridRS(DynamicPointCloudRS<REAL,SIZE> *dpc, float h)
    : m_dpc(dpc)
    , m_h(h)
    , m_hinv(1.0f/h)
  { }
  
  void init()
  {
    // determine the number of voxels needed
    AlignedBox3f &bbox = m_dpc->get_bbox();

    bmin = bbox.corner(Eigen::AlignedBox3f::BottomLeftFloor);
    bmax = bbox.corner(Eigen::AlignedBox3f::TopRightCeil);

    SIZE nx = 1 + (bmax[0] - bmin[0])*m_hinv;
    SIZE ny = 1 + (bmax[1] - bmin[1])*m_hinv;
    SIZE nz = 1 + (bmax[2] - bmin[2])*m_hinv;
    
    m_data.resize( boost::extents[nx][ny][nz] );

    SIZE num_vtx = m_dpc->get_num_vertices();
    for ( SIZE i = 0; i < num_vtx; ++i )
    {
      Vector3R<REAL> pos( m_dpc->pos_at(i) );
      Vector3R<REAL> vel( m_dpc->vel_at(i) );
      m_data(get_voxel_index(pos)).push_back( ParticleDataR<REAL>(pos, vel, m_dpc->accel_at(i)) );
    }

    compute_initial_density();
  }

  ~UniformGridRS() { }

  inline Index get_voxel_index(const Vector3R<REAL> &pos)
  {
    return {{ static_cast<typename Grid::index>(m_hinv*(pos[0]-bmin[0])),
              static_cast<typename Grid::index>(m_hinv*(pos[1]-bmin[1])),
              static_cast<typename Grid::index>(m_hinv*(pos[2]-bmin[2])) }};
  }

  void compute_accel();
  void compute_initial_density();
  void compute_pressure();

  template<typename ProcessFunc>
  void compute_quantity(ProcessFunc process);

private:
  // utility function used in the constructor to get a range of two elements
  // centered at x (or 2 if x is on the boundary)
  inline IndexRange range3(SIZE x, SIZE hi)
  {
    return IndexRange(x == 0 ? 0 : x-1, x == hi-1 ? x+1 : hi-1);
  }

private:
  DynamicPointCloudRS<REAL, SIZE> *m_dpc; // main point cloud

  // array of cells containing xyzp (position and density) for each vertex
  Grid  m_data;

  float m_h;    // grid size
  float m_hinv; // 1 / h
  Vector3f bmin;
  Vector3f bmax;
}; // class UniformGridRS


// A dynamic cloud of points
template<typename REAL, typename SIZE>
class DynamicPointCloudRS : public PointCloudRS<REAL,SIZE>
{
public:
  // dynamic point cloud from a mesh
  explicit DynamicPointCloudRS(const aiMesh *pc, REAL mass, void (* update_data_callback)(void));

  // dynamic point cloud from a regular point cloud
  explicit DynamicPointCloudRS(const PointCloudRS<REAL, SIZE> &pc, REAL mass, void (* update_data_callback)(void));
  ~DynamicPointCloudRS();

  inline REAL get_mass(SIZE idx = -1) { Q_UNUSED(idx); return m_mass; }
  inline REAL get_constant()     { return m_constant; }
  inline REAL get_rest_density() { return m_rest_density; }
  inline REAL *pos_at(SIZE i)    { return this->m_pos.data() + i*3; }
  inline REAL *vel_at(SIZE i)    { return m_vel.data() + i*3; }
  inline REAL *accel_at(SIZE i)  { return m_accel.data() + i*3; }

  inline bool is_dynamic() const { return true; }

  inline void reset_accel() { m_accel.setZero(); }

  void run() 
  {
    std::cerr << "HELLOW" << std::endl; 
    if (m_update_data_callback)
      m_update_data_callback();
  }

  friend UniformGridRS<REAL,SIZE>;

protected:
  SIZE m_num_vertices;
  REAL m_mass; // uniform constant mass for each particle
  REAL m_constant;
  REAL m_rest_density;
  Matrix3XR<REAL> m_vel; // velocities
  Matrix3XR<REAL> m_accel; // accelerations
  UniformGridRS<REAL, SIZE> m_grid;

  // callback to GL class to update new positions
  void (* m_update_data_callback)(void);
}; // class DynamicPointCloudRS

// defaults
typedef DynamicPointCloudRS<double, unsigned int> DynamicPointCloud;
typedef UniformGridRS<double, unsigned int> UniformGrid;

#endif // DYNAMICS_H
