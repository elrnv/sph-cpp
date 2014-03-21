#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <vector>
#include <boost/multi_array.hpp>
#include "kernel.h"
#include "pointcloud.h"

#define M_G 9.81f

// forward declarations
template<typename REAL, typename SIZE>
class FluidRS;

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
struct FluidParticleR : public ParticleR<REAL>
{
  FluidParticleR(Vector3R<REAL> p, Vector3R<REAL> v, REAL *a, REAL *ea)
    : ParticleR<REAL>(p), vel(v), accel(a), extern_accel(ea) { }
  ~FluidParticleR() { }

  Vector3R<REAL> vel;
  REAL *accel;        // 3-array to which we write total acceleration
  REAL *extern_accel; // 3-array to which we write acceleration from external forces
  REAL temp;          // store temporary values at the particle during computation
  REAL vol;           // volume estimate
};

template<typename REAL, typename OutputType, class KernelType>
class ComputeQuantityR
{
public:
  ComputeQuantityR(float h) : kern(h) { }
  ~ComputeQuantityR() { }

protected:
  // Smoothing kernel used to interpolate data
  Kernel<OutputType, KernelType> m_kern;
};

template<typename REAL, typename OutputType, class KernelType>
class ComputeFluidQuantityR : public ComputeQuantityR<REAL, OutputType, KernelType>
{
public:
  ComputeFluidQuantityR(float h, const FluidRS *dpc) 
    : ComputeQuantityR<REAL, OutputType, KernelType>(h)
    , m_mass(dpc->get_mass())
    , m_radius(dpc->get_radius())
    , m_rest_density(dpc->get_rest_density())
    , m_viscosity(dpc->get_viscosity())
    , m_st(dpc->get_surface_tension())
    , m_cs2(dpc->get_sound_speed2())
    , m_cs(std::sqrt(m_cs2))
  { }
  ~ComputeFluidQuantityR() { }

protected:
  // global quantities acquired from the current observed object
  REAL m_mass;
  REAL m_radius;
  REAL m_rest_density;
  REAL m_viscosity;
  REAL m_st;
  REAL m_cs2;
  REAL m_cs;
};

template<typename REAL, typename OutputType, class KernelType>
class ComputeBoundaryQuantityR : public ComputeQuantityR<REAL, OutputType, KernelType>
{
public:
  ComputeBoundaryQuantityR(float h) : ComputeQuantityR<REAL, OutputType, KernelType>(h) { }
  ~ComputeBoundaryQuantityR() { }
};

typedef ComputeBoundaryQuantityR<REAL, double, Poly6Kernel> CBQPoly6;
template class CBQPoly6;

typedef ComputeFluidQuantityR<REAL, double, Poly6Kernel> CFQPoly6;
template class CFQPoly6;

// From [Solenthaler and Pajarola 2008], an alternative density
// (number_density * mass) is used
template<typename REAL>
class ComputeBoundaryVolumeR : public CBQPoly6
{
public:
  ComputeBoundaryVolumeR(float h) : CBQPoly6(h) { }
  ~ComputeBoundaryVolumeR() { }

  inline void init_particle(ParticleR<REAL> &p) { p.dinv = 0.0f; }

  inline void fluid(ParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  { }

  inline void bound(ParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    p.dinv += kern[ p.pos - near_p.pos ];
  }

  inline void finish_particle(ParticleR<REAL> &p)
  {
    p.dinv = 1.0f/(p.dinv * kern.coef);
  }
};

template<typename REAL>
class ComputeFluidDensityR : public CFQPoly6
{
public:
  ComputeFluidDensityR(float h) : CFQPoly6 { }
  ~ComputeFluidDensityR() { }

  void init(REAL m, REAL rd, REAL r, REAL &mv, REAL &av)
  {
    mass = m; rest_density = rd; radius = r; 
    max_var = &mv; avg_var = &av;
  }

  inline void init_particle(FluidParticleR<REAL> &p)
  { 
    p.dinv = 0.0f; p.vol = 0.0f;
  }
  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    p.dinv += kern[ p.pos - near_p.pos ];
    p.vol += kern[ p.pos - near_p.pos ];
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    p.dinv += rest_density * near_p.dinv * kern[ p.pos - near_p.pos ];
    p.vol += kern[ p.pos - near_p.pos ];
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
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

  inline void init_particle(FluidParticleR<REAL> &p) { p.temp = 0.0f; }

  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    p.temp += (near_p.vel - p.vel).dot(kern[ p.pos - near_p.pos ]);
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    p.temp += (-p.vel).dot(kern[ p.pos - near_p.pos ]);
    //p.dinv += rest_density * near_p.dinv * kern[ p.pos - near_p.pos ];
  }

  inline void finish_particle(FluidParticleR<REAL> &p)
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

  inline void fluid(ParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
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
class ComputeViscosityAccelR
{
public:
  ComputeViscosityAccelR(float h) : kern(h) { }
  ~ComputeViscosityAccelR() { }

  void init(REAL r, REAL m, REAL rd, REAL v, REAL s, REAL c)
  { radius = r;  mass = m; rest_density = rd; viscosity = v; st = s; cs = c; }

  inline void init_particle(FluidParticleR<REAL> &p) { Q_UNUSED(p); }

  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;
#ifdef MCG03
    Vector3R<REAL> res(
          viscosity*near_p.dinv*(near_p.vel - p.vel)*kern(p.pos - near_p.pos) // viscosity contribution
        );
#else
    Vector3R<REAL> x_ab = p.pos - near_p.pos;
    REAL vx = x_ab.dot(p.vel - near_p.vel);
    if (vx >= 0)
      return;

    REAL nu = 2*viscosity*radius;
    //REAL pab = -
    Vector3R<REAL> res(
          mass*kern(x_ab) // viscosity contribution
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

private:
#ifdef MCG03
  ViscLapKernel   kern; // used to compute viscosity force
#else
  SpikyGradKernel kern; // used to compute viscosity force
#endif
  REAL radius;
  REAL mass;
  REAL rest_density;
  REAL viscosity;
  REAL st;
  REAL cs;
};

template<typename REAL>
class ComputePressureAccelR
{
public:
  ComputePressureAccelR(float h) : b_kern(h), p_kern(h) { }
  ~ComputePressureAccelR() { }

  void init(REAL m, REAL rd, REAL v, REAL s, REAL c2)
  {  mass = m; rest_density = rd; viscosity = v; st = s; const2 = c2; }

  inline void init_particle(FluidParticleR<REAL> &p) {
    Q_UNUSED(p); 
  }

  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;
#ifdef MCG03
    Vector3R<REAL> res(
        -p.dinv*mass*(p.pressure + near_p.pressure)*0.5*near_p.dinv*p_kern(p.pos - near_p.pos) // pressure contribution
        );
#else
    Vector3R<REAL> res(
        -mass*(p.pressure*p.dinv*p.dinv + near_p.pressure*near_p.dinv*near_p.dinv)*p_kern(p.pos - near_p.pos) // pressure contribution
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
    float massb = rest_density * near_p.dinv;
    Vector3R<REAL> res(
        cs2*(massb/(mass*(mass + massb))) * (p.pos - near_p.pos) * b_kern(p.pos - near_p.pos)  // pressure contribution
        );

    for (unsigned char i = 0; i < 3; ++i)
      p.extern_accel[i] += res[i]; // copy intermediate result
  }

  inline void finish_particle(FluidParticleR<REAL> &p)
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
};


// Grid structure used to optimize computing particle properties using kernels
template<typename REAL, typename SIZE>
class UniformGridRS
{
public:
  // Definitions
  typedef std::vector< FluidParticleR<REAL> > DynamicParticles;
  typedef std::vector< ParticleR<REAL> > StaticParticles;

  struct Cell
  {
    DynamicParticles fluidvec; // fluid particles
    StaticParticles boundvec; // static boundary particles
    //DynamicParticles rigidvec; // dynamic rigid body object data
  };

  typedef boost::multi_array< Cell, 3 > Array3;
  typedef typename Array3::index        Index;
  typedef boost::array<Index, 3>        Array3Index;
  typedef typename Array3::index_range  IndexRange;

  typedef typename Array3::template array_view<3>::type GridView;

  typedef std::vector< FluidRS<REAL, SIZE> * > FluidVec;

  // Constructors/Destructor
  UniformGridRS(float cell_size, const Vector3f &bmin);
  ~UniformGridRS();
  
  void add_fluid();
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

private: // member functions
  // utility function used in the constructor to get a range of two elements
  // centered at x (or 2 if x is on the boundary)
  inline IndexRange range3(SIZE x, SIZE hi)
  {
    return IndexRange(x == 0 ? 0 : x-1, x == hi-1 ? hi : x+2 );
  }

private: // member variables
  FluidVec m_fluids; // array of simulated interacting fluids

  // array of cells containing xyzp (position and density) for each vertex
  Array3      m_grid;
  Array3Index m_gridsize;

  float m_h;    // cell size
  float m_hinv; // 1 / cell_size

  // the following two values define the axis aligned bounding box
  Vector3f m_bmin; // min boundary corner
  Vector3f m_bmax; // max boundary corner

  ComputeFluidDensityR<REAL>   m_proc_fluid_density;
  ComputeBoundaryVolumeR<REAL> m_proc_bound_volume;
  ComputePressureR<REAL>       m_proc_pressure;
  ComputePressureAccelR<REAL>  m_proc_pressure_accel;
  ComputeViscosityAccelR<REAL> m_proc_viscosity_accel;

}; // class UniformGridRS

#endif // DYNAMICS_H
