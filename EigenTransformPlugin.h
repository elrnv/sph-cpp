// The following are members of the Transform class in Eigen

// perspective projection transform
Transform<Scalar,Dim,Mode,Options>&
perspective(float verticalAngle, float aspectRatio, float nearPlane, float farPlane)
{
  EIGEN_STATIC_ASSERT(Mode!=int(Isometry), THIS_METHOD_IS_ONLY_FOR_SPECIFIC_TRANSFORMATIONS)

  // Check that projection volume is non-zero.
  if (nearPlane == farPlane || aspectRatio == 0.0f)
    return *this;

  // Construct the projection.
  MatrixType m;
  m.setIdentity();

  float radians = (verticalAngle / 2.0f) * M_PI / 180.0f;
  float sine = std::sin(radians);
  if (sine == 0.0f)
    return *this;

  float cotan = std::cos(radians) / sine;
  float clip = farPlane - nearPlane;
  m << cotan / aspectRatio, 0.0f, 0.0f, 0.0f,
    0.0f, cotan, 0.0f, 0.0f,
    0.0f, 0.0f, -(nearPlane + farPlane)/clip, -(2.0f*nearPlane*farPlane)/clip,
    0.0f, 0.0f, -1.0f, 0.0f;

  // Apply the projection.
  m_matrix *= m;
  return *this;
}

QMatrix4x4 toQMatrix4x4(void) const
{
  check_template_params();
  EIGEN_STATIC_ASSERT(Dim==3, YOU_MADE_A_PROGRAMMING_MISTAKE)

  if (Mode == AffineCompact)
    return QMatrix4x4(
        m_matrix.coeff(0,0), m_matrix.coeff(0,1), m_matrix.coeff(0,2), m_matrix.coeff(0,3), 
        m_matrix.coeff(1,0), m_matrix.coeff(1,1), m_matrix.coeff(1,2), m_matrix.coeff(1,3), 
        m_matrix.coeff(2,0), m_matrix.coeff(2,1), m_matrix.coeff(2,2), m_matrix.coeff(2,3), 
        0.0f, 0.0f, 0.0f, 1.0f);
  else
    return QMatrix4x4(
        m_matrix.coeff(0,0), m_matrix.coeff(0,1), m_matrix.coeff(0,2), m_matrix.coeff(0,3), 
        m_matrix.coeff(1,0), m_matrix.coeff(1,1), m_matrix.coeff(1,2), m_matrix.coeff(1,3), 
        m_matrix.coeff(2,0), m_matrix.coeff(2,1), m_matrix.coeff(2,2), m_matrix.coeff(2,3), 
        m_matrix.coeff(3,0), m_matrix.coeff(3,1), m_matrix.coeff(3,2), m_matrix.coeff(3,3));
}

QMatrix3x3 toQMatrix3x3(void) const
{
  check_template_params();
  EIGEN_STATIC_ASSERT(Dim==3, YOU_MADE_A_PROGRAMMING_MISTAKE)
  EIGEN_STATIC_ASSERT(!(Options & RowMajor), YOU_MADE_A_PROGRAMMING_MISTAKE)
    return QMatrix3x3(linear().data()).transposed();
}
 
