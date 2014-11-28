#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include "types.h"
//#include "quantityprocessor.h"
#include "dynparams.h"

// Fluid Stuff

// Forward declaration
class PointCloud;
class SPHGrid;

// A dynamic cloud of points

class Fluid
{
public:
  struct __attribute__ ((__packed__)) CachedFrame
  {
    Matrix3XR<Real> pos;
    Matrix3XR<Real> vel;
    bool valid;

    CachedFrame() : valid(false) { }
    CachedFrame(const Matrix3XR<Real> &p, const Matrix3XR<Real> &v, bool good)
      : pos(p), vel(v), valid(good) { }
  };

  // Define a set of cached frames
  typedef std::vector< CachedFrame > Cache;

  // dynamic point cloud from a regular updatable gl point cloud
  explicit Fluid(const aiMesh *pc, Index matidx, FluidParamsPtr params);
  ~Fluid();

  void init(const SPHGrid &grid);

  // get pointers to be able to externally evolve data
  inline Real *vel_at(Size i)    { return m_vel.data() + i*3; }
  inline Real *accel_at(Size i)  { return m_accel.data() + i*3; }
  inline Real *extern_accel_at(Size i)  { return m_extern_accel.data() + i*3; }
  inline Real *dinv_at(Size i)  { return m_dinv.data() + i; }

  inline Real &get_avg_density()  { return m_avg_density; }
  inline Real &get_avg_pressure()  { return m_avg_pressure; }

  inline Matrix3XR<Real> &get_pos()          { return m_pc.get_pos(); }
  inline const Matrix3XR<Real> &get_pos() const  { return m_pc.get_pos(); }
  inline Matrix3XR<Real> &get_vel()          { return m_vel; }
  inline Matrix3XR<Real> &get_accel()        { return m_accel; }
  inline Matrix3XR<Real> &get_extern_accel() { return m_extern_accel; }

  // boundary getters
  inline const Vector3f &get_bmin() const { return m_bmin; }
  inline const Vector3f &get_bmax() const { return m_bmax; }

  // kernel support radius
  inline Real get_kernel_radius()
  { 
    return m_params->kernel_inflation * m_pc.get_radius();
  }

  inline Size get_num_vertices() const { return m_pc.get_num_vertices(); }
  inline Real get_halo_radius()   { return get_kernel_radius(); }

  inline Real get_mass() const            { return m_mass; }
  inline Real get_rest_density() const    { return m_rest_density; }
  inline Real get_viscosity() const       { return m_viscosity; }
  inline Real get_surface_tension() const { return m_st; }
  inline Real get_recoil_velocity_damping() const { return m_recoil_velocity_damping; }
  inline Real get_sound_speed2() const    { return m_c2; }
  inline Real get_compressibility() const { return m_params->compressibility; }
  inline Real get_friction() const        { return m_params->friction; }
  inline ParticleType get_type() const    { return m_params->fluid_type; }

  friend std::size_t hash_value( const Fluid &fl ) 
  { 
    return hash_value(*(fl.m_params)); 
  }

  inline void reset_accel() { m_accel.setZero(); m_extern_accel.setZero(); }
  inline void reset_extern_accel() { m_extern_accel.setZero(); }

  bool clamp(Real &d, Real min, Real max, Real tol);
  void clamp(float adjust, float push);
  void resolve_collisions();

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

  Real m_kernel_radius;
  Real m_rest_density;
  Real m_viscosity;
  Real m_st;
  Real m_mass;
  Real m_recoil_velocity_damping; // only for impulse collision model
  Real m_c2;
  Real m_avg_density; // used in ICS13 computations
  Real m_avg_pressure; // used in ICS13 computations

  Matrix3XR<Real> m_vel;   // velocities
  Matrix3XR<Real> m_accel; // accelerations
  Matrix3XR<Real> m_extern_accel; // accelerations
  VectorXR<Real>  m_dinv;  // density inverses

  PointCloud     &m_pc;    // The underlying dynamic cloud of points

  std::string m_savefmt;

  Cache m_cache; // collection of cached frames
}; // class Fluid


typedef boost::shared_ptr< Fluid > FluidPtr;
typedef std::vector< Fluid > FluidVec;

#endif // FLUID_H
