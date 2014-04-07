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
    h = _h; h2 = h*h; h3 = h2*h; h4 = h3*h; 
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
  float h, h2, h3, h4;
};

template<typename KernelType>
using Kerneld = Kernel<double, KernelType>;

template<typename KernelType>
using Kernel3d = Kernel<Vector3d, KernelType>;

struct MKI04Kernel : public Kerneld<MKI04Kernel>
{
  inline void init_coef() { coef = 0.02; }
  
  inline double kern(const Vector3d &r)
  {
    double q = r.norm()/h;

    if (q >= 0 && q < 2)
    {
      if (3.0f*q < 2.0f)
        return 2.0f/3.0f;
      if (q < 1)
        return (2*q - (3.0f/2.0f)*q*q);

      return 0.5*pow2(2 - q);
    }

    return 0;
  }

}; // MKI04Kernel


struct CubicSplineKernel : public Kerneld<CubicSplineKernel>
{
  inline void init_coef() { coef = 1.0f/(4.0f*M_PI*h3); }
  
  inline double kern(const Vector3d &r)
  {
    double q = r.norm();

    if (q >= 0 && q < 2)
    {
      if (q <= 1)
        return pow3(2 - q) - 4.0f*pow3(1 - q);

      return pow3(2 - q);
    }

    return 0;
  }

}; // CubicSplineKernel

struct CubicSplineGradKernel : public Kernel3d<CubicSplineGradKernel>
{
  inline void init_coef() { coef = 3.0f/(4.0f*M_PI*h3); }
  
  inline Vector3d kern(const Vector3d &r)
  {
    double q = r.norm();

    if (q > 0 && q < 2)
    {
      if (q <= 1)
        return -r*((pow2(2 - q) - 4.0f*pow2(1 - q))/q);

      return -r*(pow2(2 - q)/q);
    }

    return Vector3d(0.0f, 0.0f, 0.0f);
  }

}; // CubicSplineGradKernel

struct Poly6Kernel : public Kerneld<Poly6Kernel>
{
  inline void init_coef() { coef = 315.0f/(64*M_PI*h3*h3*h3); }
  
  // main kernel without coefficient
  inline double kern(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? pow3(h2 - r2) : 0;
  }
}; // Poly6Kernel

struct Poly6GradKernel : public Kernel3d<Poly6GradKernel>
{
  inline void init_coef() { coef = 1890.0f/(64*M_PI*h3*h3*h3); }

  inline Vector3d kern(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? -r*pow2(h2-r2) : Vector3d(0,0,0);
  }
}; // Poly6GradKernel

struct Poly6LapKernel : public Kerneld<Poly6LapKernel>
{
  inline void init_coef() { coef = 1890.0f/(64*M_PI*h3*h3*h3); }

  inline double kern(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? (h2-r2)*(7*r2 - 3*h2) : 0;
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
