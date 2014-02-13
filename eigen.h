#ifndef EIGEN_H
#define EIGEN_H

// This file contains opens the part of the namespace of the Eigen linear
// algebra library used by this project, so we don't have to do this every time

#include <cmath> // needed for the EigenTransformPlugin.h below
#include <QtGui/QMatrix4x4>

// Plugin implementing the perspective tranformation
#define EIGEN_TRANSFORM_PLUGIN "EigenTransformPlugin.h"

#include <Eigen/Geometry>
using Eigen::Vector3d;
using Eigen::Vector3f;

using Eigen::Matrix3d;
using Eigen::Matrix3f;
using Eigen::Matrix4d;
using Eigen::Matrix4f;

using Eigen::Affine3d;
using Eigen::Affine3f;
using Eigen::AffineCompact3d;
using Eigen::AffineCompact3f;
using Eigen::Projective3d;
using Eigen::Projective3f;
using Eigen::Isometry3d;
using Eigen::AlignedScaling3d;
using Eigen::AngleAxisf;

using Eigen::AlignedBox3f;

// Other Eigen Extensions

// Class to represent colors
class Color : public Vector3f
{
public:
  Color(float r, float g, float b)
    : Vector3f(r, g, b) { }
  float r() { return (*this)[0]; }
  float g() { return (*this)[1]; }
  float b() { return (*this)[2]; }
};

#endif // EIGEN_H
