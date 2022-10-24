// Sphere

#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"

class Sphere
{
protected:
  Vector c; //!< Center vertex.
  double r; //!< Radius.
public:
  //! Empty.
  Sphere() {}
  explicit Sphere(const Vector&, double);
  explicit Sphere(double);

  explicit Sphere(const std::vector<Vector>&);

  //! Empty.
  ~Sphere() {}

  // Comparison
  friend int operator==(const Sphere&, const Sphere&);
  friend int operator!=(const Sphere&, const Sphere&);

  bool Inside(const Sphere&) const;
  bool Inside(const Vector&) const;

  // Acces to properties
  Vector Center() const;
  double Radius() const;

  double Volume() const;
  double Area() const;

  // Translation, scale
  void Translate(const Vector&);
  void Scale(double);

  // Wrap
  Vector Wrap(const Vector&) const;

  friend std::ostream& operator<<(std::ostream&, const Sphere&);

public:
  static const double epsilon; //!< Internal \htmlonly\epsilon;\endhtmlonly for ray intersection tests.
  static const Sphere Null; //!< Empty Sphere.
};

//! Returns the center of the sphere.
inline Vector Sphere::Center() const
{
  return c;
}

/*!
\brief Returns the radius of the sphere.
*/
inline double Sphere::Radius() const
{
  return r;
}


//! Compute the volume of a Sphere.
inline double Sphere::Volume() const
{
  return 4.0 / 3.0 * M_1_PI * r * r * r;
}

/*!
\brief Compute the surface area of a Sphere.
*/
inline double Sphere::Area() const
{
  return 4.0 * M_1_PI * r * r;
}

/*!
\brief Check if an argument Sphere is inside the Sphere.
\param Sphere The Sphere.
*/
inline bool Sphere::Inside(const Sphere& sphere) const
{
  return (Norm(sphere.c-c) + sphere.r <= r);
}

/*!
\brief Check if a point is inside the Sphere.
\param p Point.
*/
inline bool Sphere::Inside(const Vector& p) const
{
  return (Norm(p-c) < r);
}

/*!
\brief Check if two Spheres are (strictly) equal.
\param a, b Spheres.
*/
inline int operator==(const Sphere& a, const Sphere& b)
{
  return (a.c == b.c) && (a.r == b.r);
}

/*!
\brief Check if two Spheres are (strictly) different.
\param a, b Spheres.
*/
inline int operator!=(const Sphere& a, const Sphere& b)
{
  return !(a == b);
}
