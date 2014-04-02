#include "quantityprocessor.h"

// utility


// Acceleration specialization 

template<typename REAL, typename SIZE>
class CFAccelRST<REAL,SIZE,MCG03> : public CFQ<REAL,SIZE,CFAccelRST<REAL,SIZE,MCG03> >
{
public:
  inline void init_kernel(float h)
  {
    m_spikygrad_kern.init(h);
    m_visclap_kern.init(h);
    m_colorgrad_kern.init(h);
  }
  inline void init_particle(FluidParticleR<REAL> &p)
  { Q_UNUSED(p); }

  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;
    Vector3R<REAL> res(
    // pressure
         -p.dinv*this->m_mass*(p.pressure +
          near_p.pressure)*0.5*near_p.dinv*m_spikygrad_kern(p.pos -
            near_p.pos)
         +
    // viscosity
          this->m_mass * p.dinv * this->m_viscosity *
        near_p.dinv * (near_p.vel - p.vel) * this->m_visclap_kern(p.pos - near_p.pos));

    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result
  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p) { }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    p.extern_accel[1] -= M_G;
    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += p.extern_accel[i];
    //qDebug() << "p.accel = " << p.accel[0] << p.accel[1] << p.accel[2];
  }

private:
  SpikyGradKernel m_spikygrad_kern; // for pressure
  Poly6GradKernel m_colorgrad_kern;     // surface tension
  ViscLapKernel   m_visclap_kern;   // for viscosity
};

template<typename REAL, typename SIZE>
class CFAccelRST<REAL,SIZE,BT07> : public CFQ<REAL,SIZE,CFAccelRST<REAL,SIZE,BT07> >
{
public:
  inline void init_kernel(float h)
  {
    m_poly6grad_kern.init(h);
    m_spikygrad_kern.init(h);
    m_bound_kern.init(h);
  }
  inline void init_particle(FluidParticleR<REAL> &p)
  { Q_UNUSED(p); }

  inline void fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p)
  {
    if (&p == &near_p)
      return;

    // pressure

    Vector3R<REAL> res(
        -this->m_mass*(p.pressure*p.dinv*p.dinv +
          near_p.pressure*near_p.dinv*near_p.dinv)*this->m_spikygrad_kern(p.pos - near_p.pos));

    //qDebug(" ac: % 10.5e % 10.5e % 10.5e", res[0], res[1], res[2]);
    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result

    // viscosity

    Vector3R<REAL> x_ab = p.pos - near_p.pos;
    REAL vx = x_ab.dot(p.vel - near_p.vel);
    if (vx >= 0)
      return;

    REAL nu = 2*this->m_viscosity*this->m_radius;
    //REAL pab = -
    res = this->m_mass*this->m_poly6grad_kern(x_ab); // viscosity contribution

    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += res[i]; // copy intermediate result

  }
  inline void bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p)
  {
    // pressure contribution from boundary particles
    float massb = this->m_rest_density * near_p.dinv;
    Vector3R<REAL> res(
        this->m_cs2*(massb/(this->m_mass*(this->m_mass + massb))) * (p.pos - near_p.pos) * m_bound_kern(p.pos - near_p.pos)
        );

    for (unsigned char i = 0; i < 3; ++i)
      p.extern_accel[i] += res[i]; // copy intermediate result
  }
  inline void finish_particle(FluidParticleR<REAL> &p)
  {
    p.extern_accel[1] -= M_G;
    for (unsigned char i = 0; i < 3; ++i)
      p.accel[i] += p.extern_accel[i];
  }
private:
  SpikyGradKernel m_spikygrad_kern;
  Poly6GradKernel m_poly6grad_kern;
  MKI04Kernel     m_bound_kern; // for boundary particles
}; // CFAccel

template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::init_kernel(float h) { }
template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::init_particle(FluidParticleR<REAL> &p) { }
template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::fluid(FluidParticleR<REAL> &p, FluidParticleR<REAL> &near_p) { }
template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::bound(FluidParticleR<REAL> &p, ParticleR<REAL> &near_p) { }
template<typename REAL, typename SIZE, int FT>
void CFAccelRST<REAL,SIZE,FT>::finish_particle(FluidParticleR<REAL> &p) { }

template class CFAccelRST<double, unsigned int, 0>;
template class CFAccelRST<double, unsigned int, 1>;
template class CFAccelRST<double, unsigned int, 2>;
template class CFAccelRST<double, unsigned int, 3>;
