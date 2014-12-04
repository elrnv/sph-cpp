#ifndef KERNEL_H
#define KERNEL_H

#include <cmath>
#include "eigen.h"
#include "extramath.h"

template<typename InputType, typename OutputType, class KernelType>
class Kernel
{
public:
  Kernel() { }

  inline void init(InputType _h) 
  { 
    h = _h; h2 = h*h; h3 = h2*h; h4 = h3*h;
    hinv = InputType(1.0)/h; hinv2 = hinv*hinv;
    coef = KernelType::compute_coef(h);
  } 

  // by convention operator[] -> don't premultiply by coef
  inline OutputType operator[](const Vector3T<InputType> &r)
  {
    return static_cast<KernelType*>(this)->kern(r);
  }
  inline OutputType operator()(const Vector3T<InputType> &r)
  {
    return coef * static_cast<KernelType*>(this)->kern(r);
  }

  /// core kernel computation without requiring an instance
  /// this function allows the programmer to sacrifice a bit of performance for
  /// ease of use
  inline static OutputType compute_full(const Vector3T<InputType> &r, InputType h)
  {
    return KernelType::compute_coef(h) * KernelType::compute(r,h);
  }

  // main coefficient
  InputType coef;

  // precompute values
  InputType h, h2, h3, h4, hinv, hinv2;
};

template<typename KernelType>
using Kerneld = Kernel<double, double, KernelType>;

template<typename KernelType>
using Kernel3d = Kernel<double, Vector3d, KernelType>;

// This is a pseudo kernel (doesn't actually have kernel properties but it
// is used like a kernel for penalty based boundary forces)
struct MKI04Kernel : public Kerneld<MKI04Kernel>
{
  inline double kern(const Vector3d &r)
  {
    // only use the first coordinate (reusing Kernel class for something else)
    double q = r[0]*hinv; 
    return kern_base(q);
  }

  inline static double compute(const Vector3d &r, double h)
  {
    return kern_base(r[0]/h);
  }

  inline static double compute_coef(double h) { return (27.0f/(h*11.0f)); }

  friend class Kernel;
private:
  inline static double kern_base(double q)
  {
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
  inline double kern(const Vector3d &r)
  {
    return kern_base(r.norm()*hinv);
  }
  
  inline static double compute(const Vector3d &r, double h)
  {
    return kern_base(r.norm()/h);
  }

  inline static double compute_coef(double h) { return 16.0/(M_PI*pow3(h)); }
  
  friend class Kernel;
private:
  inline static double kern_base(double q)
  {
    if (q < 1)
    {
      if (q <= 0.5)
        return 3.0*pow3(q) - 3.0*pow2(q) + 0.5;

      return pow3(1.0 - q);
    }

    return 0.0;
  }
}; // CubicSplineKernel

struct CubicSplineGradKernel : public Kernel3d<CubicSplineGradKernel>
{
  inline Vector3d kern(const Vector3d &r)
  {
    double q = r.norm()*hinv;
    return kern_base(r,q);
  }
  inline static Vector3d compute(const Vector3d &r, double h)
  {
    double q = r.norm()/h;
    return kern_base(r,q);
  }

  inline static double compute_coef(double h) { return 48.0f/(M_PI*pow5(h)); }
  
  friend class Kernel;
private:
  inline static Vector3d kern_base(const Vector3d &r, double q)
  {
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
  inline double kern(const Vector3d &r)
  {
    return kern_base(r.norm()*hinv);
  }
  inline static double compute(const Vector3d &r, double h)
  {
    return kern_base(r.norm()/h);
  }

  inline static double compute_coef(double h) { return 96.0f/(M_PI*pow5(h)); }
  
  friend class Kernel;
private:
  inline static double kern_base(double q)
  {
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
  inline double kern(const Vector3d &r)
  {
    return kern_base(r.squaredNorm()*hinv2);
  }
  inline static double compute(const Vector3d &r, double h)
  {
    return kern_base(r.squaredNorm()/pow2(h));
  }

  inline static double compute_coef(double h) { return 315.0f/(64.0f*M_PI*pow3(h)); }
  
  friend class Kernel;
private:
  inline static double kern_base(double q2)
  {
    return q2 <= 1 ? pow3(1 - q2) : 0;
  }
}; // Poly6Kernel

struct Poly6GradKernel : public Kernel3d<Poly6GradKernel>
{
  inline Vector3d kern(const Vector3d &r)
  {
    return kern_base(r, r.squaredNorm()*hinv2);
  }
  inline static Vector3d compute(const Vector3d &r, double h)
  {
    return kern_base(r, r.squaredNorm()/pow2(h));
  }

  inline static double compute_coef(double h) { return 945.0f/(32.0f*M_PI*pow5(h)); }

  friend class Kernel;
private:
  inline static Vector3d kern_base(const Vector3d &r, double q2)
  {
    return q2 <= 1 ? -r*pow2(1-q2) : Vector3d(0,0,0);
  }
}; // Poly6GradKernel

struct Poly6LapKernel : public Kerneld<Poly6LapKernel>
{
  inline double kern(const Vector3d &r)
  {
    return kern_base(r.squaredNorm()*hinv2);
  }
  inline static double compute(const Vector3d &r, double h)
  {
    return kern_base(r.squaredNorm()/pow2(h));
  }

  inline static double compute_coef(double h) { return 945.0f/(32.0f*M_PI*pow5(h)); }

  friend class Kernel;
private:
  inline static double kern_base(double q2)
  {
    return q2 <= 1 ? (1.0f-q2)*(7*q2 - 3.0f) : 0;
  }
}; // Poly6LapKernel



struct SpikyKernel : public Kerneld<SpikyKernel>
{
  inline double kern(const Vector3d &r)
  {
    return compute(r, h);
  }
  inline static double compute(const Vector3d &r, double h)
  {
    double rn = r.norm();
    return rn <= h ? pow3(h - rn) : 0;
  }

  inline static double compute_coef(double h) { return 15.0f/(M_PI*pow6(h)); }

  friend class Kernel;
}; // SpikyKernel

struct SpikyGradKernel : public Kernel3d<SpikyGradKernel>
{
  // gradient of kernel
  inline Vector3d kern(const Vector3d &r)
  {
    double rn = r.norm();
    if (rn == 0.0f) // handle degeneracy
      return Vector3d(h2,h2,h2); // chose one sided limit
    return rn <= h ? -r*(pow2(h - rn) / rn) : Vector3d(0,0,0);
  }

  inline static Vector3d compute(const Vector3d &r, double h)
  {
    double rn = r.norm();
    if (rn == 0.0f) // handle degeneracy
    {
      double h2 = pow2(h);
      return Vector3d(h2,h2,h2); // chose one sided limit
    }
    return rn <= h ? -r*(pow2(h - rn) / rn) : Vector3d(0,0,0);
  }

  inline static double compute_coef(double h) { return 45.0f/(M_PI*pow6(h)); }

  friend class Kernel;
}; // SpikyGradKernel

// laplacian of spiky kernel
struct SpikyLapKernel : public Kerneld<SpikyLapKernel>
{
  inline double kern(const Vector3d &r)
  {
    return compute(r, h);
  }
  inline static double compute(const Vector3d &r, double h)
  {
    double rn = r.norm();
    if (rn > h)
      return 0;
    else
    { 
      return (h-rn)*(2-h/rn);
    }
  }

  inline static double compute_coef(double h) { return 90.0f/(M_PI*pow6(h)); }
  
  friend class Kernel;
}; // SpikyLapKern



struct ViscKernel : public Kerneld<ViscKernel>
{
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
  inline static double compute(const Vector3d &r, double h)
  {
    double r2 = r.squaredNorm();
    double h2 = pow2(h);
    if (r2 >= h2)
      return 0;
    else
    { 
      double h4 = pow2(h2);
      double rn = std::sqrt(r2);
      return (rn*2*h*(r2-h2) - r2*r2 + h4)/(2*h2*h*rn);
    }
  }

  inline static double compute_coef(double h) { return 15.0f/(2*M_PI*pow3(h)); }

  friend class Kernel;
}; // ViscKernel

struct ViscGradKernel : public Kernel3d<ViscGradKernel>
{
  inline Vector3d kern(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    double rn = std::sqrt(r2);
    double r3 = r2*rn;
    return r2 < h2 ? Vector3d(r*((4.0f*h*r3 - 3.f*r2*r2 - h4) / (2.f*h3*r3))) : Vector3d(0,0,0);
  }
  inline static Vector3d compute(const Vector3d &r, double h)
  {
    double r2 = r.squaredNorm();
    double rn = std::sqrt(r2);
    double r3 = r2*rn;
    double h3 = pow3(h);
    return rn < h ? Vector3d(r*((4.0f*h*r3 - 3.f*r2*r2 - h3*h) / (2.f*h3*r3))) : Vector3d(0,0,0);
  }

  inline static double compute_coef(double h) { return 15.0f/(2*M_PI*pow3(h)); }

  friend class Kernel;
}; // ViscGradKernel

struct ViscLapKernel : public Kerneld<ViscLapKernel>
{
  inline double kern(const Vector3d &r)
  {
    return compute(r,h);
  }
  inline static double compute(const Vector3d &r, double h)
  {
    double rn = r.norm();
    return rn <= h ? (h-rn) : 0;
  }

  inline static double compute_coef(double h) { return 45.0f/(M_PI*pow6(h)); }

  friend class Kernel;
}; // ViscLapKernel


#endif // KERNEL_H
