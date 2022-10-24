#pragma once

#include <math.h>
#include <ostream>

class Math
{
public:
  static constexpr double Clamp(double, double = 0.0, double = 1.0);

  // Minimum and maximum
  static constexpr double Min(double, double);
  static constexpr double Max(double, double);
  static constexpr double Min(double, double, double);
  static constexpr double Max(double, double, double);

  static constexpr double DegreeToRadian(double);
  static constexpr double RadianToDegree(double);
};

/*!
\brief Clamp a double value between two bounds.
\param x Input value.
\param a, b Lower and upper bounds.
*/
inline constexpr double Math::Clamp(double x, double a, double b)
{
  return (x < a ? a : (x > b ? b : x));
}

/*!
\brief Minimum of two reals.
\param a, b Real values.
*/
inline constexpr double Math::Min(double a, double b)
{
  return (a < b ? a : b);
}

/*!
\brief Maximum of two reals.
\param a, b Real values.
*/
inline constexpr double Math::Max(double a, double b)
{
  return (a > b ? a : b);
}

/*!
\brief Maximum of three reals.
\param a, b, c Real values.
*/
inline constexpr double Math::Max(double a, double b, double c)
{
  return Math::Max(Math::Max(a, b), c);
}

/*!
\brief Minimum of three reals.
\param a, b, c Real values.
*/
inline constexpr double Math::Min(double a, double b, double c)
{
  return Math::Min(Math::Min(a, b), c);
}

/*!
\brief Convert degrees to randians.
\param a Angle in degrees.
*/
inline constexpr double Math::DegreeToRadian(double a)
{
  return a * 3.14159265358979323846 / 180.0;
}

/*!
\brief Convert radian to degrees.
\param a Angle in radian.
*/
inline constexpr double Math::RadianToDegree(double a)
{
  return a * 180.0 / 3.14159265358979323846;
}

// Class
class Vector
{
protected:
  double c[3]; //!< Components.
public:
  //! Empty
  Vector() {}

  explicit Vector(double);
  explicit Vector(double, double, double);

  // Access members
  double& operator[] (int);
  double operator[] (int) const;

  // Unary operators
  Vector operator+ () const;
  Vector operator- () const;

  // Assignment operators
  Vector& operator+= (const Vector&);
  Vector& operator-= (const Vector&);
  Vector& operator*= (const Vector&);
  Vector& operator/= (const Vector&);
  Vector& operator*= (double);
  Vector& operator/= (double);

  // Binary operators
  friend int operator> (const Vector&, const Vector&);
  friend int operator< (const Vector&, const Vector&);

  friend int operator>= (const Vector&, const Vector&);
  friend int operator<= (const Vector&, const Vector&);

  // Binary operators
  friend Vector operator+ (const Vector&, const Vector&);
  friend Vector operator- (const Vector&, const Vector&);

  friend constexpr double operator* (const Vector&, const Vector&);

  friend Vector operator* (const Vector&, double);
  friend Vector operator* (double, const Vector&);
  friend Vector operator/ (const Vector&, double);

  friend Vector operator/ (const Vector&, const Vector&);

  // Boolean functions
  friend int operator==(const Vector&, const Vector&);
  friend int operator!=(const Vector&, const Vector&);

  // Norm
  friend double Norm(const Vector&);
  friend double SquaredNorm(const Vector&);

  friend void Normalize(Vector&);
  friend Vector Normalized(const Vector&);

  // Compare functions
  static Vector Min(const Vector&, const Vector&);
  static Vector Max(const Vector&, const Vector&);

  // Abs
  friend Vector Abs(const Vector&);

  // Orthogonal and orthonormal vectors
  Vector Orthogonal() const;
  void Orthonormal(Vector&, Vector&) const;

  friend Vector Lerp(const Vector&, const Vector&, double);
  static Vector Bilinear(const Vector&, const Vector&, const Vector&, const Vector&, double, double);

  // Scale
  Vector Scaled(const Vector&) const;
  Vector Inverse() const;

  friend std::ostream& operator<<(std::ostream&, const Vector&);

public:
  static const Vector Null; //!< Null vector.
  static const Vector X; //!< Vector(1,0,0).
  static const Vector Y; //!< Vector(0,1,0).
  static const Vector Z; //!< Vector(0,0,1).
};

/*!
\brief Create a vector with the same coordinates.
\param a Real.
*/
inline Vector::Vector(double a)
{
  c[0] = c[1] = c[2] = a;
}

/*!
\brief Create a vector with argument coordinates.
\param a,b,c Coordinates.
*/
inline Vector::Vector(double a, double b, double c)
{
  Vector::c[0] = a;
  Vector::c[1] = b;
  Vector::c[2] = c;
}

//! Gets the i-th coordinate of vector.
inline double& Vector::operator[] (int i)
{
  return c[i];
}

//! Returns the i-th coordinate of vector.
inline double Vector::operator[] (int i) const
{
  return c[i];
}

// Unary operators

//! Overloaded.
inline Vector Vector::operator+ () const
{
  return *this;
}

//! Overloaded.
inline Vector Vector::operator- () const
{
  return Vector(-c[0], -c[1], -c[2]);
}

// Assignment unary operators

//! Destructive addition.
inline Vector& Vector::operator+= (const Vector& u)
{
  c[0] += u.c[0]; c[1] += u.c[1]; c[2] += u.c[2];
  return *this;
}

//! Destructive subtraction.
inline Vector& Vector::operator-= (const Vector& u)
{
  c[0] -= u.c[0]; c[1] -= u.c[1]; c[2] -= u.c[2];
  return *this;
}

//! Destructive scalar multiply.
inline Vector& Vector::operator*= (double a)
{
  c[0] *= a; c[1] *= a; c[2] *= a;
  return *this;
}

/*!
\brief Scale a vector.
\param a Scaling vector.
*/
inline Vector Vector::Scaled(const Vector& a) const
{
  return Vector(c[0] * a[0], c[1] * a[1], c[2] * a[2]);
}

/*!
\brief Inverse of a vector.

This function inverses the components of the vector. This is the same as:
\code
Vector v=Vector(1.0/u[0],1.0/u[1],1.0/u[2]);
\endcode
*/
inline Vector Vector::Inverse() const
{
  return Vector(1.0 / c[0], 1.0 / c[1], 1.0 / c[2]);
}

//! Destructive division by a scalar.
inline Vector& Vector::operator/= (double a)
{
  c[0] /= a; c[1] /= a; c[2] /= a;
  return *this;
}

/*!
\brief Destructively scale a vector by another vector.

This is the same as Scale:
\code
Vector u(2.0,-1.0,1.0);
u=u.Scaled(Vector(3.0,1.0,2.0)); // u*=Vector(3.0,1.0,2.0);
\endcode
*/
inline Vector& Vector::operator*= (const Vector& u)
{
  c[0] *= u.c[0]; c[1] *= u.c[1]; c[2] *= u.c[2];
  return *this;
}

//! Destructively divide the components of a vector by another vector.
inline Vector& Vector::operator/= (const Vector& u)
{
  c[0] /= u.c[0]; c[1] /= u.c[1]; c[2] /= u.c[2];
  return *this;
}

//! Compare two vectors.
inline int operator> (const Vector& u, const Vector& v)
{
  return ((u.c[0] > v.c[0]) && (u.c[1] > v.c[1]) && (u.c[2] > v.c[2]));
}

//! Compare two vectors.
inline int operator< (const Vector& u, const Vector& v)
{
  return ((u.c[0] < v.c[0]) && (u.c[1] < v.c[1]) && (u.c[2] < v.c[2]));
}

//! Overloaded
inline int operator>= (const Vector& u, const Vector& v)
{
  return ((u.c[0] >= v.c[0]) && (u.c[1] >= v.c[1]) && (u.c[2] >= v.c[2]));
}

//! Overloaded
inline int operator<= (const Vector& u, const Vector& v)
{
  return ((u.c[0] <= v.c[0]) && (u.c[1] <= v.c[1]) && (u.c[2] <= v.c[2]));
}

//! Adds up two vectors.
inline Vector operator+ (const Vector& u, const Vector& v)
{
  return Vector(u.c[0] + v.c[0], u.c[1] + v.c[1], u.c[2] + v.c[2]);
}

//! Difference between two vectors.
inline Vector operator- (const Vector& u, const Vector& v)
{
  return Vector(u.c[0] - v.c[0], u.c[1] - v.c[1], u.c[2] - v.c[2]);
}

//! Scalar product.
inline constexpr double operator* (const Vector& u, const Vector& v)
{
  return (u.c[0] * v.c[0] + u.c[1] * v.c[1] + u.c[2] * v.c[2]);
}

//! Right multiply by a scalar.
inline Vector operator* (const Vector& u, double a)
{
  return Vector(u.c[0] * a, u.c[1] * a, u.c[2] * a);
}

//! Left multiply by a scalar.
inline Vector operator* (double a, const Vector& v)
{
  return v * a;
}

//! Cross product.
inline Vector operator/ (const Vector& u, const Vector& v)
{
  return Vector(u.c[1] * v.c[2] - u.c[2] * v.c[1], u.c[2] * v.c[0] - u.c[0] * v.c[2], u.c[0] * v.c[1] - u.c[1] * v.c[0]);
}

//! Left multiply by a scalar
inline Vector operator/ (const Vector& u, double a)
{
  return Vector(u.c[0] / a, u.c[1] / a, u.c[2] / a);
}

// Boolean functions

//! Strong equality test.
inline int operator== (const Vector& u, const Vector& v)
{
  return ((u.c[0] == v.c[0]) && (u.c[1] == v.c[1]) && (u.c[2] == v.c[2]));
}

//! Strong difference test.
inline int operator!= (const Vector& u, const Vector& v)
{
  return (!(u == v));
}

/*!
\brief Compute the Euclidean norm of a vector.

This function involves a square root computation, it is in general more efficient to rely on
the squared norm of a vector instead.
\param u %Vector.
\sa SquaredNorm
*/
inline double Norm(const Vector& u)
{
  return sqrt(u.c[0] * u.c[0] + u.c[1] * u.c[1] + u.c[2] * u.c[2]);
}

/*!
\brief Compute the squared Euclidean norm of a vector.
\param u %Vector.
\sa Norm
*/
inline double SquaredNorm(const Vector& u)
{
  return (u.c[0] * u.c[0] + u.c[1] * u.c[1] + u.c[2] * u.c[2]);
}

/*!
\brief Return a normalized vector.

Compute the inverse of its norm and scale the components.

This function does not check if the vector is null.
\param u %Vector.
*/
inline Vector Normalized(const Vector& u)
{
  return u * (1.0 / Norm(u));
}

/*!
\brief Computes the absolute value of a vector.
\param u %Vector.
*/
inline Vector Abs(const Vector& u)
{
  return Vector(u[0] > 0.0 ? u[0] : -u[0], u[1] > 0.0 ? u[1] : -u[1], u[2] > 0.0 ? u[2] : -u[2]);
}

/*!
\brief Return a vector with coordinates set to the minimum coordinates
of the two argument vectors.
*/
inline Vector Vector::Min(const Vector& a, const Vector& b)
{
  return Vector(a[0] < b[0] ? a[0] : b[0], a[1] < b[1] ? a[1] : b[1], a[2] < b[2] ? a[2] : b[2]);
}

/*!
\brief Return a vector with coordinates set to the maximum coordinates
of the two argument vectors.
*/
inline Vector Vector::Max(const Vector& a, const Vector& b)
{
  return Vector(a[0] > b[0] ? a[0] : b[0], a[1] > b[1] ? a[1] : b[1], a[2] > b[2] ? a[2] : b[2]);
}

/*!
\brief Linear interpolation between two vectors.
\param a,b Interpolated points.
\param t Interpolant.
*/
inline Vector Lerp(const Vector& a, const Vector& b, double t)
{
  return a + t * (b - a);
}

/*!
\brief Bi-linear interpolation between four vectors.

The values are given in trigonometric order.

\param a00,a10,a11,a01 Interpolated vectors.
\param u,v Interpolation coefficients.

\sa Math::Bilinear
*/
inline Vector Vector::Bilinear(const Vector& a00, const Vector& a10, const Vector& a11, const Vector& a01, double u, double v)
{
  return (1 - u) * (1 - v) * a00 + (1 - u) * (v)*a01 + (u) * (1 - v) * a10 + (u) * (v)*a11;
}



// Class Matrix (for homotheties and rotations)

class Matrix
{
protected:
  double comp[3][3]; //!< Components.
  bool isRotation; //!< True if the matrix is a rotation matrix.
  bool isHomothety; //!< True if the matrix is a homothety matrix.
public:
  //! Empty
  Matrix() {}

  // Common matrix initializations
  explicit Matrix(double (*)[3], bool, bool);
  explicit Matrix(double, double, double, double, double, double, double, double, double, bool, bool);

  // Rotation matrices
  explicit Matrix(const Vector&, double);

  // Homothetie matrices
  explicit Matrix(double);


  // Access members
  double* operator[] (int);
  const double* operator[] (int) const;


  // Binary operators
  friend Matrix operator* (const Matrix&, const Matrix&);

  friend Matrix operator* (const Matrix&, double);
  friend Vector operator* (const Matrix&, const Vector&);
  friend Matrix operator* (double, const Matrix&);
  friend Matrix operator/ (const Matrix&, double);


  // Norm (Fronebenius)
  friend double Norm(const Matrix&);
  friend double SquaredNorm(const Matrix&);

  // Transpose and inverse
  friend Matrix Transpose(const Matrix&);
  friend Matrix Inverse(const Matrix&);


  friend std::ostream& operator<<(std::ostream&, const Matrix&);

public:
  static const Matrix Null; //!< Null matrix.
  static const Matrix Id; //!< Identity matrix.
};

/*!
\brief Create a matrix.
*/
inline Matrix::Matrix(double a[3][3], bool r, bool h)
{
  comp[0][0] = a[0][0];
  comp[0][1] = a[0][1];
  comp[0][2] = a[0][2];
  comp[1][0] = a[1][0];
  comp[1][1] = a[1][1];
  comp[1][2] = a[1][2];
  comp[2][0] = a[2][0];
  comp[2][1] = a[2][1];
  comp[2][2] = a[2][2];
  isRotation = r;
  isHomothety = h;
}

/*!
\brief Create a matrix.
*/
inline Matrix::Matrix(double a, double b, double c, double d, double e, double f, double g, double h, double i, bool rot, bool hom)
{
  comp[0][0] = a;
  comp[0][1] = b;
  comp[0][2] = c;
  comp[1][0] = d;
  comp[1][1] = e;
  comp[1][2] = f;
  comp[2][0] = g;
  comp[2][1] = h;
  comp[2][2] = i;
  isRotation = rot;
  isHomothety = hom;
}

/*!
\brief Create a matrix with the same coordinates.
\param a Real.
*/
inline Matrix::Matrix(double a)
{
  isRotation = false;
  isHomothety = true;
  comp[0][0] = comp[1][1] = comp[2][2] = a;
  comp[0][1] = comp[0][2] = comp[1][0] = comp[1][2] = comp[2][0] = comp[2][1] = 0.0;
}

/*!
\brief Create a rotation matrix from an angle and an axe.
\param a Real.
\param axe Axe of rotation.
*/
inline Matrix::Matrix(const Vector& axe, double a)
{
  isRotation = true;
  isHomothety = false;
  Vector u = Normalized(axe);
  double c = cos(a);
  double s = sin(a);
  double t = 1.0 - c;
  double x = u[0];
  double y = u[1];
  double z = u[2];

  comp[0][0] = t * x * x + c;
  comp[0][1] = t * x * y - s * z;
  comp[0][2] = t * x * z + s * y;
  comp[1][0] = t * x * y + s * z;
  comp[1][1] = t * y * y + c;
  comp[1][2] = t * y * z - s * x;
  comp[2][0] = t * x * z - s * y;
  comp[2][1] = t * y * z + s * x;
  comp[2][2] = t * z * z + c;
}

//! Gets the i-th coordinate of matrix.
inline double* Matrix::operator[] (int i)
{
  return comp[i];
}

//! Returns the i-th coordinate of matrix.
inline const double* Matrix::operator[] (int i) const
{
  return comp[i];
}

// Unary operators


//! Matrix product.
inline Matrix operator* (const Matrix& u, const Matrix& v)
{
  bool r = u.isRotation && v.isRotation;
  bool h = u.isHomothety && v.isHomothety;
  double tmp_comp[3][3];
  tmp_comp[0][0] = u.comp[0][0] * v.comp[0][0] + u.comp[0][1] * v.comp[1][0] + u.comp[0][2] * v.comp[2][0];
  tmp_comp[0][1] = u.comp[0][0] * v.comp[0][1] + u.comp[0][1] * v.comp[1][1] + u.comp[0][2] * v.comp[2][1];
  tmp_comp[0][2] = u.comp[0][0] * v.comp[0][2] + u.comp[0][1] * v.comp[1][2] + u.comp[0][2] * v.comp[2][2];
  tmp_comp[1][0] = u.comp[1][0] * v.comp[0][0] + u.comp[1][1] * v.comp[1][0] + u.comp[1][2] * v.comp[2][0];
  tmp_comp[1][1] = u.comp[1][0] * v.comp[0][1] + u.comp[1][1] * v.comp[1][1] + u.comp[1][2] * v.comp[2][1];
  tmp_comp[1][2] = u.comp[1][0] * v.comp[0][2] + u.comp[1][1] * v.comp[1][2] + u.comp[1][2] * v.comp[2][2];
  tmp_comp[2][0] = u.comp[2][0] * v.comp[0][0] + u.comp[2][1] * v.comp[1][0] + u.comp[2][2] * v.comp[2][0];
  tmp_comp[2][1] = u.comp[2][0] * v.comp[0][1] + u.comp[2][1] * v.comp[1][1] + u.comp[2][2] * v.comp[2][1];
  tmp_comp[2][2] = u.comp[2][0] * v.comp[0][2] + u.comp[2][1] * v.comp[1][2] + u.comp[2][2] * v.comp[2][2];
  return Matrix(tmp_comp, r, h);
}

inline Vector operator* (const Matrix& u, const Vector& v)
{
  return Vector(u.comp[0][0] * v[0] + u.comp[0][1] * v[1] + u.comp[0][2] * v[2],
                u.comp[1][0] * v[0] + u.comp[1][1] * v[1] + u.comp[1][2] * v[2],
                u.comp[2][0] * v[0] + u.comp[2][1] * v[1] + u.comp[2][2] * v[2]);
}

//! Right multiply by a scalar.
inline Matrix operator* (const Matrix& u, double a)
{
  return Matrix(u.comp[0][0] * a, u.comp[0][1] * a, u.comp[0][2] * a, u.comp[1][0] * a, u.comp[1][1] * a, u.comp[1][2] * a, u.comp[2][0] * a, u.comp[2][1] * a, u.comp[2][2] * a, u.isRotation, u.isHomothety);
}

//! Left multiply by a scalar.
inline Matrix operator* (double a, const Matrix& u)
{
  return u * a;
}


//! Left multiply by a scalar
inline Matrix operator/ (const Matrix& u, double a)
{
  return u * (1.0 / a);
}


/*!
\brief Compute the Euclidean norm of a matrix.

This function involves a square root computation, it is in general more efficient to rely on
the squared norm of a matrix instead.
\param u %Matrix.
\sa SquaredNorm
*/
inline double Norm(const Matrix& u)
{
  return sqrt(SquaredNorm(u));
}

/*!
\brief Compute the squared Euclidean norm of a matrix.
\param u %Matrix.
\sa Norm
*/
inline double SquaredNorm(const Matrix& u)
{
  return u.comp[0][0] * u.comp[0][0] + u.comp[0][1] * u.comp[0][1] + u.comp[0][2] * u.comp[0][2] + u.comp[1][0] * u.comp[1][0] + u.comp[1][1] * u.comp[1][1] + u.comp[1][2] * u.comp[1][2] + u.comp[2][0] * u.comp[2][0] + u.comp[2][1] * u.comp[2][1] + u.comp[2][2] * u.comp[2][2];
}


inline Matrix Transpose(const Matrix& u)
{
  return Matrix(u.comp[0][0], u.comp[1][0], u.comp[2][0], u.comp[0][1], u.comp[1][1], u.comp[2][1], u.comp[0][2], u.comp[1][2], u.comp[2][2], u.isRotation, u.isHomothety);
}

inline Matrix Inverse(const Matrix& u)
{
  if (u.isHomothety){
    return Matrix(1.0 / u.comp[0][0], 0.0, 0.0, 0.0, 1.0 / u.comp[1][1], 0.0, 0.0, 0.0, 1.0 / u.comp[2][2], u.isRotation, u.isHomothety);
  }
  else if (u.isRotation){
    return Transpose(u);
  }
  else {
    throw std::invalid_argument("Matrix::Inverse: not implemented");
  }
}
