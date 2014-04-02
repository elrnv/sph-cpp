#include "quantityprocessor.h"

// utility


// Acceleration specialization 


#if 0
template<typename REAL, typename SIZE>
class CFAccelRST<REAL,SIZE,BT07> : public CFQ<REAL,SIZE,CFAccelRST<REAL,SIZE,BT07> >
{
public:
  CFAccelRST(float h)
  {
  //  m_poly6_kern.init(h);
  //  m_poly6grad_kern.init(h);
  //  m_spiky_kern.init(h);
  //  m_bound_kern.init(h);
  }
  inline void init_kernel(float h)
  {
    m_poly6_kern.init(h);
    m_poly6grad_kern.init(h);
    m_spiky_kern.init(h);
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
          near_p.pressure*near_p.dinv*near_p.dinv)*this->m_spiky_kern(p.pos - near_p.pos));

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
  Poly6Kernel     m_poly6_kern;
  Poly6GradKernel m_poly6grad_kern;
  SpikyGradKernel m_spiky_kern;
  MKI04Kernel     m_bound_kern; // for boundary particles
}; // CFAccel
#endif

