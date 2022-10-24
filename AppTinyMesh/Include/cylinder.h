// Cylinder

#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"

class Cylinder
{
protected:
  Vector c; //!< Center vertex.
  double h; //!< Height.
  double r; //!< Radius.
public:
  //! Empty.
  Cylinder() {}
  explicit Cylinder(const Vector&, double, double);
  explicit Cylinder(double, double);
  explicit Cylinder(const Vector&);

  //! Empty.
  ~Cylinder() {}

  // Comparison
  friend int operator==(const Cylinder&, const Cylinder&);
  friend int operator!=(const Cylinder&, const Cylinder&);

  // Acces to properties
  Vector Center() const;
  double Height() const;
  double Radius() const;

  double Volume() const;
  double Area() const;

  // Translation, scale
  void Translate(const Vector&);
  void Scale(double);

  friend std::ostream& operator<<(std::ostream&, const Cylinder&);

public:
  static const double epsilon; //!< Internal \htmlonly\epsilon;\endhtmlonly for ray intersection tests.
  static const Cylinder Null; //!< Empty Cylinder.
};

//! Returns the center of the cylinder.
inline Vector Cylinder::Center() const
{
  return c;
}

/*!
\brief Returns the height of the cylinder.
*/
inline double Cylinder::Height() const
{
  return h;
}

/*!
\brief Returns the radius of the cylinder.
*/
inline double Cylinder::Radius() const
{
  return r;
}


//! Compute the volume of a Cylinder.
inline double Cylinder::Volume() const
{
  return M_1_PI * r * r * h;
}

/*!
\brief Compute the surface area of a Cylinder.
*/
inline double Cylinder::Area() const
{
  return 2.0 * M_1_PI * r * (r + h);
}



/*!
\brief Check if two Cylinders are (strictly) equal.
\param a, b Cylinders.
*/
inline int operator==(const Cylinder& a, const Cylinder& b)
{
  return (a.c == b.c) && (a.r == b.r) && (a.h == b.h);
}

/*!
\brief Check if two Cylinders are (strictly) different.
\param a, b Cylinders.
*/
inline int operator!=(const Cylinder& a, const Cylinder& b)
{
  return !(a == b);
}
