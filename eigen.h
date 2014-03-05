#ifndef EIGEN_H
#define EIGEN_H

// This file contains opens the part of the namespace of the Eigen linear
// algebra library used by this project, so we don't have to do this every time

#include <cmath> // needed for the EigenTransformPlugin.h below
#include <QtGui/QMatrix4x4>

// Plugin implementing the perspective tranformation
#define EIGEN_TRANSFORM_PLUGIN "EigenTransformPlugin.h"

// One radian measured in degrees (conversion constant)
#define RADIAN 0.017453292519943

#include <Eigen/Geometry>
using Eigen::Vector4d;
using Eigen::Vector4f;
using Eigen::Vector3d;
using Eigen::Vector3f;
using Eigen::Vector2d;
using Eigen::Vector2f;

template<typename REAL>
using Vector3R = Eigen::Matrix<REAL, 3, 1>;

using Eigen::Matrix3d;
using Eigen::Matrix3f;
using Eigen::Matrix4d;
using Eigen::Matrix4f;

using Eigen::Matrix3Xd;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::NoChange;

using Eigen::Affine3d;
using Eigen::Affine3f;
using Eigen::AffineCompact3d;
using Eigen::AffineCompact3f;
using Eigen::Projective3d;
using Eigen::Projective3f;
using Eigen::Isometry3d;
using Eigen::Isometry3f;
using Eigen::AlignedScaling3d;
using Eigen::AlignedScaling3f;
using Eigen::AngleAxisf;

using Eigen::AlignedBox3f;

// Other Eigen Extensions

// Class to represent colors
class Color : public Vector4f
{
public:
  Color(float r, float g, float b, float a)
    : Vector4f(r, g, b, a) { }
  float r() { return (*this)[0]; }
  float g() { return (*this)[1]; }
  float b() { return (*this)[2]; }
  float a() { return (*this)[3]; }
};

#endif // EIGEN_H
