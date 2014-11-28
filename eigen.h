#ifndef EIGEN_H
#define EIGEN_H

// This file contains opens the part of the namespace of the Eigen linear
// algebra library used by this project, so we don't have to do this every time

#include <cmath> // needed for the EigenTransformPlugin.h below
#include <QtGui/QMatrix4x4>

// Plugin implementing custom tranformations
#define EIGEN_TRANSFORM_PLUGIN "EigenTransformPlugin.h"

// #define EIGEN_USE_MKL_ALL

#include <Eigen/Geometry>
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::NoChange;

using Eigen::Vector4d;
using Eigen::Vector4f;
using Eigen::Vector3d;
using Eigen::Vector3f;
using Eigen::Vector3i;
using Eigen::Vector2d;
using Eigen::Vector2f;

template<typename T>
using VectorXT = Matrix<T, Dynamic, 1>;

template<typename T>
using Vector3T = Matrix<T, 3, 1>;
template<typename REAL>
using Vector4T = Matrix<T, 3, 1>;

using Eigen::Matrix3d;
using Eigen::Matrix3f;
using Eigen::Matrix3Xd;
using Eigen::Matrix3Xf;
using Eigen::Matrix3Xi;
using Eigen::Matrix4d;
using Eigen::Matrix4f;

template<typename T>
using Matrix3XT = Matrix<T, 3, Dynamic>;
template<typename T>
using Matrix4XT = Matrix<T, 3, Dynamic>;

using Eigen::PermutationMatrix;

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
using Eigen::Scaling;
using Eigen::AngleAxisf;
using Eigen::Translation3f;

using Eigen::AlignedBox3f;

#endif // EIGEN_H
