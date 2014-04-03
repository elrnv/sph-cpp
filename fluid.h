#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include "pointcloud.h"
#include "quantityprocessor.h"
#include "dynparams.h"

//#define MCG03
//
//#if defined MCG03
//  #define BOUNDARY_IMPULSE
//  #define IDEAL_GAS_PRESSURE
//  #define MCG03_VISCOSITY_FORCE
//  #define MCG03_PRESSURE_FORCE
//  #define VELOCITY_DAMPING 0.4f
//#elif defined BT07
//  #define BOUNDARY_PARTICLE
//  #define BOUNDARY_PENALTY_FORCE
//
//#endif

//#define REPORT_DENSITY_VARIATION

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
  // dynamic point cloud from a regular updatable gl point cloud
  explicit FluidRS(const PointCloudRS<REAL,SIZE> *pc, FluidParamsPtr params);
  explicit FluidRS(const aiMesh *pc, FluidParamsPtr params);
  ~FluidRS();

  // explicitly state that we use some base class members (convenience)
  using PointCloudRS<REAL,SIZE>::get_radius;
  using PointCloudRS<REAL,SIZE>::pos_at;
  using PointCloudRS<REAL,SIZE>::m_pos;
  using PointCloudRS<REAL,SIZE>::m_bbox;

  void init(GLPointCloudRS<REAL, SIZE> *glpc);

  // get pointers to be able to externally evolve data
  inline REAL *vel_at(SIZE i)    { return m_vel.data() + i*3; }
  inline REAL *accel_at(SIZE i)  { return m_accel.data() + i*3; }
  inline REAL *extern_accel_at(SIZE i)  { return m_extern_accel.data() + i*3; }
  inline REAL *dinv_at(SIZE i)  { return m_dinv.data() + i; }

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
  inline FluidType get_type() const       { return m_params->fluid_type; }

  inline Vector3f get_color() const 
  { 
    const QVector3D &diffuse = m_glpc->get_diffuse();
    return Vector3f(diffuse[0], diffuse[1], diffuse[2]);
  }

  inline bool is_dynamic() const { return true; }

  inline void reset_accel() { m_accel.setZero(); m_extern_accel.setZero(); }

  inline bool clamp(REAL &d, REAL min, REAL max, REAL tol);
  inline void clamp();
  inline void resolve_collisions();

  inline void update_data(); // propagate changes to some viewer

  template<int FT>
  inline FluidRST<REAL,SIZE,FT> *cast() 
  {
    return static_cast<FluidRST<REAL,SIZE,FT> *>(this); 
  }

  void clear_cache();
  inline void write_cache(unsigned int frame);
  inline bool is_cached(unsigned int frame);
  inline bool read_cache(unsigned int frame);

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

  Matrix3XR<REAL> m_vel;   // velocities
  Matrix3XR<REAL> m_accel; // accelerations
  Matrix3XR<REAL> m_extern_accel; // accelerations
  VectorXR<REAL>  m_dinv;  // density inverses

  GLPointCloudRS<REAL,SIZE> *m_glpc; // should only used for callback

  std::string m_cachefmt;
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
  using FluidRS<REAL,SIZE>::m_rest_density;
  using FluidRS<REAL,SIZE>::m_viscosity;
  using FluidRS<REAL,SIZE>::m_st;
  using FluidRS<REAL,SIZE>::m_mass;
  using FluidRS<REAL,SIZE>::m_recoil_velocity_damping;
  using FluidRS<REAL,SIZE>::m_c2;

  inline void init_processors();

  // compile time type checked proc getters
#define GET_PROC(proc_type) \
  template <typename ProcessPairFunc> \
    inline typename std::enable_if< \
    std::is_same< ProcessPairFunc, proc_type<REAL,SIZE,FT>  >::value, \
    proc_type<REAL,SIZE,FT> >::type &get_proc()

  GET_PROC( CFDensityRST ) { return m_fluid_density_proc; }
  GET_PROC( CFDensityUpdateRST ) { return m_fluid_density_update_proc; }
  GET_PROC( CFAccelRST ) { return m_fluid_accel_proc; }

  // Quantity Processors
  CFDensityRST<REAL,SIZE,FT>              m_fluid_density_proc;
  CFDensityUpdateRST<REAL,SIZE,FT>        m_fluid_density_update_proc;
  CFAccelRST<REAL,SIZE,FT>             m_fluid_accel_proc;
}; // FluidRST

template<int FT>
using FluidT = FluidRST<double, unsigned int, FT>;

#endif // FLUID_H
