#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include "pointcloud.h"
#include "quantityprocessor.h"
#include "dynparams.h"

// Fluid Stuff

// Forward declaration
template<typename REAL, typename SIZE>
class GLPointCloudRS;
template<typename REAL, typename SIZE, int FT>
class FluidRST;

// A dynamic cloud of points, this is a purely abstract class
template<typename REAL, typename SIZE>
class FluidRS : public PointCloudRS<REAL,SIZE>
{
public:
  // explicitly state that we use some base class members (convenience)
  using PointCloudRS<REAL,SIZE>::get_radius;
  using PointCloudRS<REAL,SIZE>::pos_at;
  using PointCloudRS<REAL,SIZE>::m_pos;
  using PointCloudRS<REAL,SIZE>::m_bbox;

  struct __attribute__ ((__packed__)) CachedFrame
  {
    Matrix3XR<REAL> pos;
    Matrix3XR<REAL> vel;
    bool valid;

    CachedFrame() : valid(false) { }
    CachedFrame(const Matrix3XR<REAL> &p, const Matrix3XR<REAL> &v, bool good)
      : pos(p), vel(v), valid(good) { }
  };

  // Define a set of cached frames
  typedef std::vector< CachedFrame > Cache;

  // dynamic point cloud from a regular updatable gl point cloud
  explicit FluidRS(const PointCloudRS<REAL,SIZE> *pc, FluidParamsPtr params);
  explicit FluidRS(const aiMesh *pc, FluidParamsPtr params);
  ~FluidRS();

  void init(GLPointCloudRS<REAL, SIZE> *glpc);

  // get pointers to be able to externally evolve data
  inline REAL *vel_at(SIZE i)    { return m_vel.data() + i*3; }
  inline REAL *accel_at(SIZE i)  { return m_accel.data() + i*3; }
  inline REAL *extern_accel_at(SIZE i)  { return m_extern_accel.data() + i*3; }
  inline REAL *dinv_at(SIZE i)  { return m_dinv.data() + i; }

  inline REAL &get_avg_density()  { return m_avg_density; }
  inline REAL &get_avg_pressure()  { return m_avg_pressure; }

  inline Matrix3XR<REAL> &get_vel()          { return m_vel; }
  inline Matrix3XR<REAL> &get_accel()        { return m_accel; }
  inline Matrix3XR<REAL> &get_extern_accel() { return m_extern_accel; }

  // boundary getters
  inline const Vector3f &get_bmin() const { return m_bmin; }
  inline const Vector3f &get_bmax() const { return m_bmax; }

  // kernel support radius
  inline REAL get_kernel_radius()
  { 
    return m_params->kernel_inflation * this->get_radius(); 
  }

  inline REAL get_halo_radius()   { return get_kernel_radius(); }

  inline REAL get_mass() const            { return m_mass; }
  inline REAL get_rest_density() const    { return m_rest_density; }
  inline REAL get_viscosity() const       { return m_viscosity; }
  inline REAL get_surface_tension() const { return m_st; }
  inline REAL get_recoil_velocity_damping() const { return m_recoil_velocity_damping; }
  inline REAL get_sound_speed2() const    { return m_c2; }
  inline REAL get_compressibility() const { return m_params->compressibility; }
  inline REAL get_friction() const        { return m_params->friction; }
  inline FluidType get_type() const       { return m_params->fluid_type; }

  friend std::size_t hash_value( const FluidRS<REAL,SIZE> &fl ) 
  { return hash_value(*(fl.m_params)); }

  inline Vector3f get_color() const 
  { 
    const QVector3D &diffuse = m_glpc->get_diffuse();
    return Vector3f(diffuse[0], diffuse[1], diffuse[2]);
  }

  inline bool is_dynamic() const { return true; }

  inline void reset_accel() { m_accel.setZero(); m_extern_accel.setZero(); }
  inline void reset_extern_accel() { m_extern_accel.setZero(); }

  bool clamp(REAL &d, REAL min, REAL max, REAL tol);
  void clamp(float adjust, float push);
  void resolve_collisions();

  void update_data(); // propagate changes to some viewer

  template<int FT>
  inline FluidRST<REAL,SIZE,FT> *cast() 
  {
    return static_cast<FluidRST<REAL,SIZE,FT> *>(this); 
  }

  void clear_saved();
  void save(unsigned int frame);
  bool is_saved(unsigned int frame);
  bool load_saved_cache();
  bool read_saved(unsigned int frame);

  void clear_cache();
  void cache(unsigned int frame);
  bool is_cached(unsigned int frame);
  bool load_cached(unsigned int frame);

protected:
  FluidParamsPtr m_params;

  Vector3f m_bmin;
  Vector3f m_bmax;

  REAL m_kernel_radius;
  REAL m_rest_density;
  REAL m_viscosity;
  REAL m_st;
  REAL m_mass;
  REAL m_recoil_velocity_damping; // only for impulse collision model
  REAL m_c2;
  REAL m_avg_density; // used in ICS13 computations
  REAL m_avg_pressure; // used in ICS13 computations

  Matrix3XR<REAL> m_vel;   // velocities
  Matrix3XR<REAL> m_accel; // accelerations
  Matrix3XR<REAL> m_extern_accel; // accelerations
  VectorXR<REAL>  m_dinv;  // density inverses

  GLPointCloudRS<REAL,SIZE> *m_glpc; // should only used for callback

  std::string m_savefmt;

  Cache m_cache;
}; // class FluidRS


// Typed Fluid
template<typename REAL, typename SIZE, int FT>
class FluidRST : public FluidRS<REAL,SIZE>
{
public:
  explicit FluidRST(const PointCloudRS<REAL,SIZE> *pc, FluidParamsPtr params);
  explicit FluidRST(const aiMesh *pc, FluidParamsPtr params);
  ~FluidRST();

  // explicitly state that we use some base class members (convenience)

  using FluidRS<REAL,SIZE>::m_kernel_radius;
  void init_processors();

  // compile time type checked proc getters
#define GET_PROC(proc_type) \
  template <typename ProcessPairFunc> \
    inline typename std::enable_if< \
    std::is_same< ProcessPairFunc, proc_type<REAL,SIZE,FT>  >::value, \
    proc_type<REAL,SIZE,FT> >::type &get_proc()

  GET_PROC( CFDensityRST )        { return m_fluid_density_proc; }
  GET_PROC( CFDensityUpdateRST )  { return m_fluid_density_update_proc; }
  GET_PROC( CFAccelRST )          { return m_fluid_accel_proc; }
  GET_PROC( CFSurfaceNormalRST )  { return m_fluid_surface_normal_proc; }
  GET_PROC( CFSurfaceTensionRST ) { return m_fluid_surface_tension_proc; }
  GET_PROC( CFPrepareJacobiRST )  { return m_fluid_prepare_jacobi_proc; }
  GET_PROC( CFJacobiSolveFirstRST )  { return m_fluid_jacobi_solve1_proc; }
  GET_PROC( CFJacobiSolveSecondRST ) { return m_fluid_jacobi_solve2_proc; }
  GET_PROC( CFPressureAccelRST )     { return m_fluid_pressure_accel_proc; }

  // Quantity Processors
  CFDensityRST<REAL,SIZE,FT>        m_fluid_density_proc;
  CFDensityUpdateRST<REAL,SIZE,FT>  m_fluid_density_update_proc;
  CFAccelRST<REAL,SIZE,FT>          m_fluid_accel_proc;
  CFSurfaceNormalRST<REAL,SIZE,FT>  m_fluid_surface_normal_proc;
  CFSurfaceTensionRST<REAL,SIZE,FT> m_fluid_surface_tension_proc;
  CFPrepareJacobiRST<REAL,SIZE,FT>  m_fluid_prepare_jacobi_proc;
  CFJacobiSolveFirstRST<REAL,SIZE,FT>   m_fluid_jacobi_solve1_proc;
  CFJacobiSolveSecondRST<REAL,SIZE,FT>  m_fluid_jacobi_solve2_proc;
  CFPressureAccelRST<REAL,SIZE,FT>      m_fluid_pressure_accel_proc;
}; // FluidRST

template<int FT>
using FluidT = FluidRST<double, unsigned int, FT>;

#endif // FLUID_H
