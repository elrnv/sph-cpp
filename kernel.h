#ifndef KERNEL_H
#define KERNEL_H

#include <cmath>
#include "eigen.h"

template<typename OutputType, class KernelType>
class Kernel
{
public:
  Kernel() { }

  inline void init(float _h) 
  { 
    h = _h; h2 = h*h; h3 = h2*h; h4 = h3*h; hinv = 1.0f/h; hinv2 = hinv*hinv;
    static_cast<KernelType*>(this)->init_coef();
  } 

  // by convention operator[] -> don't premultiply by coef
  inline OutputType operator[](const Vector3d &r)
  {
    return static_cast<KernelType*>(this)->kern(r);
  }
  inline OutputType operator()(const Vector3d &r)
  {
    return coef * static_cast<KernelType*>(this)->kern(r);
  }

  // main coefficient
  float coef;

  inline double pow3(double x) { return x*x*x; }
  inline double pow2(double x) { return x*x; }

  // precompute values
  float h, h2, h3, h4, hinv, hinv2;
};

template<typename KernelType>
using Kerneld = Kernel<double, KernelType>;

template<typename KernelType>
using Kernel3d = Kernel<Vector3d, KernelType>;

struct MKI04Kernel : public Kerneld<MKI04Kernel>
{
  inline void init_coef() { coef = (27.0f/11.0f)*hinv; }
  
  inline double kern(const Vector3d &r)
  {
    // only use the first coordinate (reusing Kernel class for something else)
    double q = r[0]*hinv; 

    if (q < 1)
    {
      if (3.0f*q < 1)
        return 1.0f/3.0f;
      if (q < 0.5)
        return 2.0f*q - 3.0f*q*q;

      return pow2(1 - q);
    }

    return 0.0f;
  }

}; // MKI04Kernel


struct CubicSplineKernel : public Kerneld<CubicSplineKernel>
{
  inline void init_coef() { coef = 16.0f/(M_PI*h3); }
  
  inline double kern(const Vector3d &r)
  {
    double q = r.norm()*hinv;

    if (q < 1)
    {
      if (q <= 0.5)
        return 3.0f*pow3(q) - 3.0f*pow2(q) + 0.5f;

      return pow3(1 - q);
    }

    return 0;
  }

}; // CubicSplineKernel

struct CubicSplineGradKernel : public Kernel3d<CubicSplineGradKernel>
{
  inline void init_coef() { coef = 48.0f/(M_PI*h4*h); }
  
  inline Vector3d kern(const Vector3d &r)
  {
    double q = r.norm()*hinv;

    if (q > 0 && q < 1)
    {
      if (q <= 0.5)
        return r*(3.0f*q - 2.0f);

      return -r*(pow2(1 - q)/q);
    }

    return Vector3d(0.0f, 0.0f, 0.0f);
  }

}; // CubicSplineGradKernel

struct CubicSplineLapKernel : public Kerneld<CubicSplineLapKernel>
{
  inline void init_coef() { coef = 96.0f/(M_PI*h4*h); }
  
  inline double kern(const Vector3d &r)
  {
    double q = r.norm()*hinv;

    if (q > 0 && q < 1)
    {
      if (q <= 0.5)
        return 3.0f*(2.0f*q - 1.0f);

      return (1.0f - 2*q)*(q - 1.0f)/q;
    }

    return 0;
  }

}; // CubicSplineLapKernel

struct Poly6Kernel : public Kerneld<Poly6Kernel>
{
  inline void init_coef() { coef = 315.0f/(64*M_PI*h3); }
  
  // main kernel without coefficient
  inline double kern(const Vector3d &r)
  {
    double q2 = r.squaredNorm()*hinv2;
    return q2 <= 1 ? pow3(1 - q2) : 0;
  }
}; // Poly6Kernel

struct Poly6GradKernel : public Kernel3d<Poly6GradKernel>
{
  inline void init_coef() { coef = 1890.0f/(64*M_PI*h3*h2); }

  inline Vector3d kern(const Vector3d &r)
  {
    double q2 = r.squaredNorm()*hinv2;
    return q2 <= 1 ? -r*pow2(1-q2) : Vector3d(0,0,0);
  }
}; // Poly6GradKernel

struct Poly6LapKernel : public Kerneld<Poly6LapKernel>
{
  inline void init_coef() { coef = 1890.0f/(64*M_PI*h3*h2); }

  inline double kern(const Vector3d &r)
  {
    double q2 = r.squaredNorm()*hinv2;
    return q2 <= 1 ? (1.0f-q2)*(7*q2 - 3.0f) : 0;
  }
}; // Poly6LapKernel



struct SpikyKernel : public Kerneld<SpikyKernel>
{
  inline void init_coef() { coef = 15.0f/(M_PI*h3*h3); }

  inline double kern(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? pow3(h - std::sqrt(r2)) : 0;
  }
}; // SpikyKernel

struct SpikyGradKernel : public Kernel3d<SpikyGradKernel>
{
  inline void init_coef() { coef = 45.0f/(M_PI*h3*h3); }

  // gradient of kernel
  inline Vector3d kern(const Vector3d &r)
  {
    double rn = r.norm();
    if (rn == 0.0f) // handle degeneracy
      return Vector3d(h2,h2,h2); // chose one sided limit

    return rn <= h ? -r*(pow2(h - rn) / rn) : Vector3d(0,0,0);
  }
}; // SpikyGradKernel

// laplacian of spiky kernel
struct SpikyLapKernel : public Kerneld<SpikyLapKernel>
{
  inline void init_coef() { coef = 90.0f/(M_PI*h3*h3); }

  inline double kern(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    if (r2 > h2)
      return 0;
    else
    { 
      double rn = std::sqrt(r2);
      return (h-rn)*(2-h/rn);
    }
  }
}; // SpikyLapKern



struct ViscKernel : public Kerneld<ViscKernel>
{
  inline void init_coef() { coef = 15.0f/(2*M_PI*h3); }

  inline double kern(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    if (r2 >= h2)
      return 0;
    else
    { 
      double rn = std::sqrt(r2);
      return (rn*2*h*(r2-h2) - r2*r2 + h4)/(2*h3*rn);
    }
  }
}; // ViscKernel

struct ViscGradKernel : public Kernel3d<ViscGradKernel>
{
  inline void init_coef() { coef = 15.0f/(2*M_PI*h3); }

  inline Vector3d kern(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    double rn = std::sqrt(r2);
    double r3 = r2*rn;
    return r2 < h2 ? Vector3d(r*((4.0f*h*r3 - 3.f*r2*r2 - h4) / (2.f*h3*r3))) : Vector3d(0,0,0);
  }
}; // ViscGradKernel

struct ViscLapKernel : public Kerneld<ViscLapKernel>
{
  inline void init_coef() { coef = 45.0f/(M_PI*h3*h3); }

  inline double kern(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? (h-std::sqrt(r2)) : 0;
  }
}; // ViscLapKernel


#endif // KERNEL_H
