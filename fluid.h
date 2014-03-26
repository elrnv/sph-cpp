#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include "pointcloud.h"
#include "quantityprocessor.h"

#define MCG03

#if defined MCG03
  #define BOUNDARY_IMPULSE
  #define IDEAL_GAS_PRESSURE
  #define MCG03_VISCOSITY_FORCE
  #define MCG03_PRESSURE_FORCE
  #define VELOCITY_DAMPING 0.4f
#elif defined BT07
  #define BOUNDARY_PARTICLE
  #define BOUNDARY_PENALTY_FORCE

#endif

//#define REPORT_DENSITY_VARIATION

#include "dynamics.h"

#define M_G 9.81f

// Compute SPH Quantity Processors
template<typename REAL, typename SIZE>
class CFDensityRS :
  public CFQPoly6RS< REAL, SIZE, CFDensityRS<REAL,SIZE> >
{
public:
  CFDensityRS(float h) 
    : CFQPoly6RS<REAL, SIZE, CFDensityRS<REAL,SIZE> >(h) { }

  void init(REAL &mv, REAL &av) { max_var = &mv; avg_var = &av; }

  inline void init_particle(FluidParticleR<REAL> &p)
  { 
    p.dinv = 0.0f; p.vol = 0.0f;
  }
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    p.dinv += this->m_kern[ p.pos - near_p.pos ];
    p.vol += this->m_kern[ p.pos - near_p.pos ];
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
#ifdef BOUNDARY_PARTICLE
    p.dinv += this->m_rest_density * near_p.dinv * this->m_kern[ p.pos - near_p.pos ];
    p.vol += this->m_kern[ p.pos - near_p.pos ];
#endif
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
#ifdef DENSITY_KERNEL_CORRECTION
    p.dinv =
      (8*this->m_radius*this->m_radius*this->m_radius*p.vol)/(this->m_mass * p.dinv);
#if 0
    qDebug() << (this->m_mass * p.dinv) << "/"
    << (8*this->m_radius*this->m_radius*this->m_radius*p.vol) << " = " <<
      (this->m_mass * p.dinv)/(8*this->m_radius*this->m_radius*this->m_radius*p.vol);
#endif
#else
    p.dinv = 1.0f/(this->m_mass * p.dinv * this->m_kern.coef);
#if 0
    qDebug() << (this->m_mass * p.dinv * this->m_kern.coef);
#endif
#endif

    REAL var = std::abs(1.0f/p.dinv - this->m_rest_density);
    if ( var > *max_var )
      *max_var = var;
    *avg_var += var;
  }

private:
  REAL *max_var; // max variation
  REAL *avg_var; // average variation
};

template<typename REAL, typename SIZE>
class CFDensityUpdateRS : 
  public CFQPoly6GradRS<REAL,SIZE,CFDensityUpdateRS<REAL,SIZE> >
{
public:
  CFDensityUpdateRS(float h) 
    : CFQPoly6GradRS<REAL,SIZE,CFDensityUpdateRS<REAL, SIZE> >(h) { }

  void init(float ts) { m_timestep = ts; }

  inline void init_particle(FluidParticleR<REAL> &p)
  { p.temp = 0.0f; }
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    p.temp += (near_p.vel - p.vel).dot(this->m_kern[ p.pos - near_p.pos ]);
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    p.temp += (-p.vel).dot(this->m_kern[ p.pos - near_p.pos ]);
    //p.dinv += rest_density * near_p.dinv * kern[ p.pos - near_p.pos ];
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    p.dinv += m_timestep * p.dinv * p.dinv * p.temp * this->m_kern.coef;
    qDebug() << "p.dinv: " << p.dinv;
  }

private:
  float m_timestep;
};

template<typename REAL, typename SIZE>
class CFPressureRS :
  public CFQPoly6RS<REAL,SIZE,CFPressureRS<REAL,SIZE> >
{
public:
  CFPressureRS(float h) 
    : CFQPoly6RS<REAL,SIZE,CFPressureRS<REAL,SIZE> >(h) { }

  inline REAL pow7(REAL x) { return (x*x)*(x*x)*(x*x)*x; }

  inline void init_particle(ParticleR<REAL> &p)
  { 
//    p.dinv = 0.0f;
  }
  inline void fluid(ParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
 //   p.dinv += this->m_kern[ p.pos - near_p.pos ];
  }
  inline void bound(ParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
  //  p.dinv += this->m_rest_density * near_p.dinv * this->m_kern[ p.pos - near_p.pos ];
  }
  inline void finish_particle(ParticleR<REAL> &p)
  {
    //REAL density = this->m_mass * p.dinv * this->m_kern.coef;
    //p.dinv = 1.0f/density;
#if defined TAIT_PRESSURE
    p.pressure =
      this->m_rest_density * this->m_cs2 * 0.14285714285714 * 
      (pow7(1.0f / (p.dinv * this->m_rest_density)) - 1);
#elif defined IDEAL_GAS_PRESSURE
    p.pressure =
      this->m_cs2 *
      (1.0f / p.dinv  -  this->m_rest_density);
#endif
  //  qDebug() << "density:" << 1.0f/p.dinv << " rest_density:" <<
   //   this->m_rest_density << "  pressure:" << p.pressure;
  }
};

#ifdef MCG03_VISCOSITY_FORCE
  CFQ_TYPEDEF CFQViscKernRS = CFQViscLapRS<REAL,SIZE,ComputeType>; // used to compute viscosity force
#else
  CFQ_TYPEDEF CFQViscKernRS = CFQSpikyGradRS<REAL,SIZE,ComputeType>; // used to compute viscosity force
#endif

template<typename REAL, typename SIZE>
class CFViscosityAccelRS : 
  public CFQViscKernRS<REAL,SIZE,CFViscosityAccelRS<REAL,SIZE> >
{
public:
  CFViscosityAccelRS(float h) 
    : CFQViscKernRS<REAL,SIZE,CFViscosityAccelRS<REAL,SIZE> >(h) { }

  inline void init_particle(FluidParticleR<REAL> &p)
  { Q_UNUSED(p); }
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;
#ifdef MCG03_VISCOSITY_FORCE
    Vector3R<REAL> res(
          this->m_mass*p.dinv*this->m_viscosity*near_p.dinv*(near_p.vel - p.vel)*this->m_kern(p.pos - near_p.pos) // viscosity contribution
        );
#else
    Vector3R<REAL> x_ab = p.pos - near_p.pos;
    REAL vx = x_ab.dot(p.vel - near_p.vel);
    if (vx >= 0)
      return;

    REAL nu = 2*this->m_viscosity*this->m_radius;
    //REAL pab = -
    Vector3R<REAL> res(
          this->m_mass*this->m_kern(x_ab) // viscosity contribution
        );
#endif

    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    //Vector3R<REAL> res(
    //    + viscosity*near_p.dinv*(near_p.vel - p.vel)*kern(p.pos - near_p.pos) // viscosity contribution
    //    );

    //for (unsigned char i = 0; i < 3; ++i)
    //  p.extern_accel[i] += res[i]; // copy intermediate result
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    //for (unsigned char i = 0; i < 3; ++i)
    //  p.accel[i] = p.accel[i];

//    qDebug() << "p.dinv:" << p.dinv;
    //qDebug(" ac: % 10.5e % 10.5e % 10.5e", p.accel[0], p.accel[1], p.accel[2]);

  }
};

template<typename REAL, typename SIZE>
class CFPressureAccelRS :
  public CFQSpikyGradRS<REAL,SIZE,CFPressureAccelRS<REAL,SIZE> >
{
public:
  CFPressureAccelRS(float h) 
    : CFQSpikyGradRS<REAL,SIZE,CFPressureAccelRS<REAL,SIZE> >(h)
    , m_bound_kern(h) { }

  inline void init_particle(FluidParticleR<REAL> &p)
  { Q_UNUSED(p); }
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;
#ifdef MCG03_PRESSURE_FORCE
    Vector3R<REAL> res(
        -p.dinv*this->m_mass*(p.pressure +
          near_p.pressure)*0.5*near_p.dinv*this->m_kern(p.pos - near_p.pos) // pressure contribution
        );
#else
    Vector3R<REAL> res(
        -this->m_mass*(p.pressure*p.dinv*p.dinv +
          near_p.pressure*near_p.dinv*near_p.dinv)*this->m_kern(p.pos - near_p.pos) // pressure contribution
        );
#endif

    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
   // Vector3R<REAL> res(
   //     -rest_density*near_p.dinv*p.pressure*p.dinv*p.dinv*p_kern[p.pos - near_p.pos] // pressure contribution
   //     );
#ifdef BOUNDARY_PENALTY_FORCE
    float massb = this->m_rest_density * near_p.dinv;
    Vector3R<REAL> res(
        this->m_cs2*(massb/(this->m_mass*(this->m_mass + massb))) * (p.pos - near_p.pos) * m_bound_kern(p.pos - near_p.pos)  // pressure contribution
        );

    for (unsigned char i = 0; i < 3; ++i)
      p.extern_accel[i] += res[i]; // copy intermediate result
#endif
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    p.extern_accel[1] -= M_G;
    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += p.extern_accel[i];

//    qDebug() << "p.dinv:" << p.dinv;
    //qDebug(" ac: % 10.5e % 10.5e % 10.5e", p.accel[0], p.accel[1], p.accel[2]);
  }
private:
  MKI04Kernel m_bound_kern;
};


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
  explicit FluidRS(GLPointCloudRS<REAL, SIZE> *glpc,
      REAL density, REAL viscosity, REAL surfacetension);
  ~FluidRS();

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
  Vector3f m_bmin;
  Vector3f m_bmax;

  REAL m_kernel_radius;
  REAL m_rest_density;
  REAL m_viscosity;
  REAL m_st;
  REAL m_mass;
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
    typename std::enable_if< \
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
