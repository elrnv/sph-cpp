#ifndef EIGEN_H
#define EIGEN_H

// This file contains opens the part of the namespace of the Eigen linear
// algebra library used by this project, so we don't have to do this every time

#include <Eigen/Geometry>
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using Eigen::Affine3d;
using Eigen::AffineCompact3d;
using Eigen::Projective3d;
using Eigen::Isometry3d;
using Eigen::AlignedScaling3d;

using Eigen::AlignedBox3f;

class Color : public Eigen::Vector3f
{
public:
  Color(float r, float g, float b)
    : Eigen::Vector3f(r, g, b) { }
  float r() { return (*this)[0]; }
  float g() { return (*this)[1]; }
  float b() { return (*this)[2]; }
};

#endif // EIGEN_H
