#ifndef KERNEL_H
#define KERNEL_H

#include <cmath>
#include "eigen.h"

class Kernel
{
public:
  Kernel(float h, float coef) 
    : coef(coef), h(h), h2(h*h), h3(h2*h), h4(h3*h) { } 
  ~Kernel() { }

  // main coefficient
  float coef;

protected:
  inline float pow3(double x) { return x*x*x; }
  inline float pow2(double x) { return x*x; }

  // precompute values
  float h, h2, h3, h4;
};

// two types of kernel functions ( scalar and vector ) for convenience
struct KernelScalar : public Kernel
{
  KernelScalar(float h, float coef) : Kernel(h, coef) { }
  ~KernelScalar() {}
  // by convention operator[] -> don't premultiply by coef
  virtual inline float operator[](const Vector3d &r) = 0;
  virtual inline float operator()(const Vector3d &r) = 0;
};

struct KernelVector : public Kernel
{
  KernelVector(float h, float coef) : Kernel(h, coef) { }
  ~KernelVector() {}
  // by convention operator[] -> don't premultiply by coef
  virtual inline Vector3d operator[](const Vector3d &r) = 0;
  virtual inline Vector3d operator()(const Vector3d &r) = 0;
};



struct Poly6Kernel : public KernelScalar
{
  Poly6Kernel(float h) : KernelScalar(h, 315.0f/(64*M_PI*h3)) { }
  ~Poly6Kernel() {}
  
  // main kernel without coefficient
  inline float operator[](const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? pow3(h2 - r2) : 0;
  }

  // main kernel
  inline float operator()(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? coef*pow3(h2 - r2) : 0;
  }
}; // Poly6Kernel

struct Poly6GradKernel : public KernelVector
{
  Poly6GradKernel(float h) : KernelVector(h, 1890.0f/(64*M_PI*h3)) { }
  ~Poly6GradKernel() { }

  inline Vector3d operator[](const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? -r*pow2(h2-r2) : Vector3d(0.0f, 0.0f, 0.0f);
  }

  inline Vector3d operator()(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? -r*coef*pow2(h2-r2) : Vector3d(0.0f, 0.0f, 0.0f);
  }
}; // Poly6GradKernel

struct Poly6LapKernel : public KernelScalar
{
  Poly6LapKernel(float h) : KernelScalar(h, 1890.0f/(64*M_PI*h3)) { }
  ~Poly6LapKernel() {}

  inline float operator[](const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? (h2-r2)*(7*r2 - 3*h2) : 0;
  }
  inline float operator()(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? coef*(h2-r2)*(7*r2 - 3*h2) : 0;
  }
}; // Poly6LapKernel



struct SpikyKernel : public KernelScalar
{
  SpikyKernel(float h) : KernelScalar(h, 15.0f/(M_PI*h3*h3)) { }
  ~SpikyKernel() {}

  inline float operator[](const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? pow3(h - std::sqrt(r2)) : 0;
  }

  inline float operator()(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? coef*pow3(h - std::sqrt(r2)) : 0;
  }
}; // SpikyKernel

struct SpikyGradKernel : public KernelVector
{
  SpikyGradKernel(float h) : KernelVector(h, 45.0f/(M_PI*h3*h3)) { }
  ~SpikyGradKernel() { }

  // gradient of kernel
  inline Vector3d operator[](const Vector3d &r)
  {
    double rn = r.norm();
    return rn <= h ? Vector3d((-pow2(h - rn) / rn)*r) : Vector3d(0.0f, 0.0f, 0.0f);
  }
  inline Vector3d operator()(const Vector3d &r)
  {
    double rn = r.norm();
    return rn <= h ? Vector3d((-coef*pow2(h - rn) / rn)*r) : Vector3d(0.0f, 0.0f, 0.0f);
  }
}; // SpikyGradKernel

// laplacian of spiky kernel
struct SpikyLapKernel : public KernelScalar
{
  SpikyLapKernel(float h) : KernelScalar(h, 90.0f/(M_PI*h3*h3)) { }
  ~SpikyLapKernel() {}

  inline float operator[](const Vector3d &r)
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
  inline float operator()(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    if (r2 > h2)
      return 0;
    else
    { 
      double rn = std::sqrt(r2);
      return coef*(h-rn)*(2-h/rn);
    }
  }
}; // SpikyLapKern



struct ViscKernel : public KernelScalar
{
  ViscKernel(float h) : KernelScalar(h, 15.0f/(2*M_PI*h3)) { }
  ~ViscKernel() {}

  inline float operator[](const Vector3d &r)
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

  inline float operator()(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    if (r2 >= h2)
      return 0;
    else
    { 
      double rn = std::sqrt(r2);
      return coef*(rn*2*h*(r2-h2) - r2*r2 + h4)/(2*h3*rn);
    }
  }
}; // ViscKernel

struct ViscGradKernel : public KernelVector
{
  ViscGradKernel(float h) : KernelVector(h, 15.0f/(2*M_PI*h3)) { }
  ~ViscGradKernel() { }

  inline Vector3d operator[](const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    double rn = std::sqrt(r2);
    double r3 = r2*rn;
    return r2 < h2 ? r*(4*h*r3 - 3*r2*r2 - h4) / (2*h3*r3) : Vector3d(0.0f, 0.0f, 0.0f);
  }
  inline Vector3d operator()(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    double rn = std::sqrt(r2);
    double r3 = r2*rn;
    return r2 < h2 ? r*coef*(4*h*r3 - 3*r2*r2 - h4) / (2*h3*r3) : Vector3d(0.0f, 0.0f, 0.0f);
  }
}; // ViscGradKernel

struct ViscLapKernel : public KernelScalar
{
  ViscLapKernel(float h) : KernelScalar(h, 45.0f/(M_PI*h3*h3)) { }
  ~ViscLapKernel() { }

  inline float operator[](const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? (h-std::sqrt(r2)) : 0;
  }
  inline float operator()(const Vector3d &r)
  {
    double r2 = r.squaredNorm();
    return r2 <= h2 ? coef*(h-std::sqrt(r2)) : 0;
  }
}; // ViscLapKernel


#endif // KERNEL_H
