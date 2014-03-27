#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include "pointcloud.h"
#include "dynamics.h"
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

// A dynamic cloud of points
template<typename REAL, typename SIZE>
class FluidRS : public PointCloudRS<REAL,SIZE>
{
public:
  // dynamic point cloud from a regular updatable gl point cloud
  explicit FluidRS(const PointCloudRS<REAL,SIZE> *pc, FluidParamsPtr params);
  explicit FluidRS(const aiMesh *pc, FluidParamsPtr params);
  ~FluidRS();

  inline void init(GLPointCloudRS<REAL, SIZE> *glpc);

  // get pointers to be able to externally evolve data
  inline REAL *pos_at(SIZE i)    { return this->m_pos.data() + i*3; }
  inline REAL *vel_at(SIZE i)    { return m_vel.data() + i*3; }
  inline REAL *accel_at(SIZE i)  { return m_accel.data() + i*3; }
  inline REAL *extern_accel_at(SIZE i)  { return m_extern_accel.data() + i*3; }

  inline Matrix3XR<REAL> &get_vel()          { return m_vel; }
  inline Matrix3XR<REAL> &get_accel()        { return m_accel; }
  inline Matrix3XR<REAL> &get_extern_accel() { return m_extern_accel; }

  // boundary getters
  inline const Vector3f &get_bmin() const { return m_bmin; }
  inline const Vector3f &get_bmax() const { return m_bmax; }

  // kernel support radius
  inline REAL get_kernel_radius() const { return m_kernel_radius; }

  inline REAL get_mass() const            { return m_mass; }
  inline REAL get_rest_density() const    { return m_rest_density; }
  inline REAL get_viscosity() const       { return m_viscosity; }
  inline REAL get_surface_tension() const { return m_st; }
  inline REAL get_sound_speed2() const    { return m_c2; }

  inline Vector3f get_color() const 
  { 
    Vector4f diffuse = m_glpc->get_diffuse();
    return Vector3f(diffuse[0], diffuse[1], diffuse[2]);
  }

  inline bool is_dynamic() const { return true; }

  inline void reset_accel() { m_accel.setZero(); m_extern_accel.setZero(); }

  inline bool clamp(REAL &d, REAL min, REAL max);
  inline void resolve_collisions();

  inline void update_data(); // propagate changes to some viewer

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

  Matrix3XR<REAL> m_vel; // velocities
  Matrix3XR<REAL> m_accel; // accelerations
  Matrix3XR<REAL> m_extern_accel; // accelerations

  GLPointCloudRS<REAL,SIZE> *m_glpc; // should only used for callback

public:
  // Quantity Processors
  CFDensityRS<REAL,SIZE>        m_fluid_density_proc;
  CFDensityUpdateRS<REAL,SIZE>  m_fluid_density_update_proc;
  CFPressureRS<REAL,SIZE>       m_fluid_pressure_proc;
  CFViscosityAccelRS<REAL,SIZE> m_fluid_viscosity_accel_proc;
  CFPressureAccelRS<REAL,SIZE>  m_fluid_pressure_accel_proc;

  // routine to copy fluid properties to processors above
  template<class OutputType, class KernelType, class ComputeType>
  inline void copy_properties_to_proc(
    CFQ<REAL,SIZE,OutputType,KernelType,ComputeType> &cfq_proc);

  // compile time type checked proc getters
#define GET_PROC(proc_type) \
  template <typename ProcessPairFunc> \
    inline typename std::enable_if< \
    std::is_same< ProcessPairFunc, proc_type<REAL,SIZE>  >::value, \
    proc_type<REAL,SIZE> >::type &get_proc()

  GET_PROC( CFDensityRS )
  { return m_fluid_density_proc; }

  GET_PROC( CFDensityUpdateRS )
  { return m_fluid_density_update_proc; }

  GET_PROC( CFPressureRS )
  { return m_fluid_pressure_proc; }

  GET_PROC( CFViscosityAccelRS )
  { return m_fluid_viscosity_accel_proc; }

  GET_PROC( CFPressureAccelRS )
  { return m_fluid_pressure_accel_proc; }


}; // class FluidRS

// defaults
typedef FluidRS<double, unsigned int> Fluid;

#endif // FLUID_H
