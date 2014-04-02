#include "quantityprocessor.h"

// Density

template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::init(REAL &mv, REAL &av)
{
  max_var = &mv; avg_var = &av;
}
template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{ 
  p.dinv = 0.0f; p.vol = 0.0f;
}
template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{
  p.dinv += this->m_kern[ p.pos - near_p.pos ];
  p.vol += this->m_kern[ p.pos - near_p.pos ];
}
template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
{
  if (FT == MCG03)
    return;

  p.dinv += this->m_rest_density * near_p.dinv * this->m_kern[ p.pos - near_p.pos ];
  p.vol += this->m_kern[ p.pos - near_p.pos ];
}
template<typename REAL, typename SIZE, int FT>
void CFDensityRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p)
{
  if (FT == MCG03)
  {
    p.dinv = 1.0f/(this->m_mass * p.dinv * this->m_kern.coef);
  }
  else
  {
    // By default, use kernel correction
    p.dinv =
      (8*this->m_radius*this->m_radius*this->m_radius*p.vol)/(this->m_mass * p.dinv);
#if 0
  qDebug() << (this->m_mass * p.dinv) << "/"
    << (8*this->m_radius*this->m_radius*this->m_radius*p.vol) << " = " <<
    (this->m_mass * p.dinv)/(8*this->m_radius*this->m_radius*this->m_radius*p.vol);
#endif
  }

  REAL var = std::abs(1.0f/p.dinv - this->m_rest_density);
  if ( var > *max_var )
    *max_var = var;
  *avg_var += var;
}


// Density Update

template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::init(float ts) 
{ 
  m_timestep = ts;
}
template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{
  p.temp = 0.0f;
}
template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{
  p.temp += (near_p.vel - p.vel).dot(this->m_kern[ p.pos - near_p.pos ]);
}
template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
{
  p.temp += (-p.vel).dot(this->m_kern[ p.pos - near_p.pos ]);
  //p.dinv += rest_density * near_p.dinv * kern[ p.pos - near_p.pos ];
}
template<typename REAL, typename SIZE, int FT>
void CFDensityUpdateRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p)
{
  p.dinv += m_timestep * p.dinv * p.dinv * p.temp * this->m_kern.coef;
}


// Pressure

template<typename REAL, typename SIZE, int FT>
REAL CFPressureRST<REAL,SIZE,FT>::pow7(REAL x)
{
  return (x*x)*(x*x)*(x*x)*x; 
}

template<typename REAL, typename SIZE, int FT>
void CFPressureRST<REAL,SIZE,FT>::init_particle(ParticleR<REAL> &p)
{ 
//    p.dinv = 0.0f;
}
template<typename REAL, typename SIZE, int FT>
void CFPressureRST<REAL,SIZE,FT>::fluid(ParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{
//   p.dinv += this->m_kern[ p.pos - near_p.pos ];
}
template<typename REAL, typename SIZE, int FT>
void CFPressureRST<REAL,SIZE,FT>::bound(ParticleR<REAL> &p, ParticleR<REAL> &near_p)
{
//  p.dinv += this->m_rest_density * near_p.dinv * this->m_kern[ p.pos - near_p.pos ];
}

template<typename REAL, typename SIZE, int FT>
void CFPressureRST<REAL,SIZE,FT>::finish_particle(ParticleR<REAL> &p)
{
  if (FT == MCG03)
  {
    p.pressure = this->m_cs2 * (1.0f/p.dinv - this->m_rest_density);
                 
  }
  else // Enable tait pressure equation by default
    p.pressure =
      this->m_rest_density * this->m_cs2 * 0.14285714285714 * 
      (pow7(1.0f / (p.dinv * this->m_rest_density)) - 1);
}


// Pressure Acceleration

template<typename REAL, typename SIZE, int FT>
void CFPressureAccelRST<REAL,SIZE,FT>::init()
{ 
  m_bound_kern.init(this->m_kern.h);
}

template<typename REAL, typename SIZE, int FT>
void CFPressureAccelRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{ 
  Q_UNUSED(p);
}

template<typename REAL, typename SIZE, int FT>
void CFPressureAccelRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{
  if (&p == &near_p)
    return;

  Vector3R<REAL> res =
    (FT == MCG03 
     ?  -p.dinv*this->m_mass*(p.pressure +
        near_p.pressure)*0.5*near_p.dinv*this->m_kern(p.pos - near_p.pos)
     : -this->m_mass*(p.pressure*p.dinv*p.dinv +
        near_p.pressure*near_p.dinv*near_p.dinv)*this->m_kern(p.pos - near_p.pos));

  //qDebug(" ac: % 10.5e % 10.5e % 10.5e", res[0], res[1], res[2]);
  for (unsigned char i = 0; i < 3; ++i)
    p.accel[i] += res[i]; // copy intermediate result

}

template<typename REAL, typename SIZE, int FT>
void CFPressureAccelRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
{
  if (FT != BT07)
    return;

  float massb = this->m_rest_density * near_p.dinv;
  Vector3R<REAL> res(
      this->m_cs2*(massb/(this->m_mass*(this->m_mass + massb))) * (p.pos - near_p.pos) * m_bound_kern(p.pos - near_p.pos)
      );

  for (unsigned char i = 0; i < 3; ++i)
    p.extern_accel[i] += res[i]; // copy intermediate result
}

template<typename REAL, typename SIZE, int FT>
void CFPressureAccelRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p)
{
  p.extern_accel[1] -= M_G;
  for (unsigned char i = 0; i < 3; ++i)
    p.accel[i] += p.extern_accel[i];

  //qDebug() << "p.dinv:" << p.dinv;
}



// Viscosity Acceleration

template<typename REAL, typename SIZE, int FT>
void CFViscosityAccelRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{
  Q_UNUSED(p);
}

template<typename REAL, typename SIZE, int FT>
void CFViscosityAccelRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{
  if (&p == &near_p)
    return;

  Vector3R<REAL> x_ab = p.pos - near_p.pos;
  REAL vx = x_ab.dot(p.vel - near_p.vel);
  if (vx >= 0)
    return;

  REAL nu = 2*this->m_viscosity*this->m_radius;
  //REAL pab = -
  Vector3R<REAL> res(this->m_mass*this->m_kern(x_ab)); // viscosity contribution

  for (unsigned char i = 0; i < 3; ++i)
    p.accel[i] += res[i]; // copy intermediate result
}

template<typename REAL, typename SIZE, int FT>
void CFViscosityAccelRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
{
  //Vector3R<REAL> res(
  //    + viscosity*near_p.dinv*(near_p.vel - p.vel)*kern(p.pos - near_p.pos) // viscosity contribution
  //    );

  //for (unsigned char i = 0; i < 3; ++i)
  //  p.extern_accel[i] += res[i]; // copy intermediate result
}
template<typename REAL, typename SIZE, int FT>
void CFViscosityAccelRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p)
{
  //for (unsigned char i = 0; i < 3; ++i)
  //  p.accel[i] = p.accel[i];

//    qDebug() << "p.dinv:" << p.dinv;
  //qDebug(" ac: % 10.5e % 10.5e % 10.5e", p.accel[0], p.accel[1], p.accel[2]);

}

// Viscosity specialization 
template<typename REAL, typename SIZE>
class CFViscosityAccelRST<REAL,SIZE,MCG03> :
  public CFQViscLapRS<REAL,SIZE,CFViscosityAccelRST<REAL,SIZE,MCG03> >
{
public:
  inline void init_particle(FluidParticleR<REAL> &p)
  { Q_UNUSED(p); }

  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    Vector3R<REAL> res(this->m_mass * p.dinv * this->m_viscosity *
        near_p.dinv * (near_p.vel - p.vel) * this->m_kern(p.pos - near_p.pos));

    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p) { }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    //qDebug() << "p.accel = " << p.accel[0] << p.accel[1] << p.accel[2];
  }
};


// Surface Tension Acceleration

template<typename REAL, typename SIZE, int FT>
void CFSurfaceTensionAccelRST<REAL,SIZE,FT>::init()
{
  m_lap_kern.init(this->m_kern.h);
}

template<typename REAL, typename SIZE, int FT>
void CFSurfaceTensionAccelRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p)
{ 
  Q_UNUSED(p);
}

template<typename REAL, typename SIZE, int FT>
void CFSurfaceTensionAccelRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
{
  if (&p == &near_p)
    return;

  Vector3R<REAL> res;

///  for (unsigned char i = 0; i < 3; ++i)
///    p.accel[i] += res[i]; // copy intermediate result
}

template<typename REAL, typename SIZE, int FT>
void CFSurfaceTensionAccelRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p) { }

template<typename REAL, typename SIZE, int FT>
void CFSurfaceTensionAccelRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p) { }

// MCG03 surface tension specialization 
template<typename REAL, typename SIZE>
class CFSurfaceTensionAccelRST<REAL,SIZE,MCG03> : 
  public CFQPoly6GradRS<REAL,SIZE,CFSurfaceTensionAccelRST<REAL,SIZE,MCG03> >
{
public:
  inline void init() { m_lap_kern.init(this->m_kern.h); }

  inline void init_particle(FluidParticleR<REAL> &p)
  { Q_UNUSED(p); }

  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;

    Vector3R<REAL> res;

 //   for (unsigned char i = 0; i < 3; ++i)
//      p.accel[i] += res[i]; // copy intermediate result
  }

  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
  }
private:
  Poly6LapKernel m_lap_kern;
};

template class CFDensityRST<double, unsigned int, 0>;
template class CFDensityRST<double, unsigned int, 1>;
template class CFDensityRST<double, unsigned int, 2>;
template class CFDensityRST<double, unsigned int, 3>;
template class CFDensityUpdateRST<double, unsigned int, 0>;
template class CFDensityUpdateRST<double, unsigned int, 1>;
template class CFDensityUpdateRST<double, unsigned int, 2>;
template class CFDensityUpdateRST<double, unsigned int, 3>;
template class CFPressureRST<double, unsigned int, 0>;
template class CFPressureRST<double, unsigned int, 1>;
template class CFPressureRST<double, unsigned int, 2>;
template class CFPressureRST<double, unsigned int, 3>;
template class CFPressureAccelRST<double, unsigned int, 0>;
template class CFPressureAccelRST<double, unsigned int, 1>;
template class CFPressureAccelRST<double, unsigned int, 2>;
template class CFPressureAccelRST<double, unsigned int, 3>;
template class CFViscosityAccelRST<double, unsigned int, 0>;
template class CFViscosityAccelRST<double, unsigned int, 1>;
template class CFViscosityAccelRST<double, unsigned int, 2>;
template class CFViscosityAccelRST<double, unsigned int, 3>;
template class CFSurfaceTensionAccelRST<double, unsigned int, 0>;
template class CFSurfaceTensionAccelRST<double, unsigned int, 1>;
template class CFSurfaceTensionAccelRST<double, unsigned int, 2>;
template class CFSurfaceTensionAccelRST<double, unsigned int, 3>;


