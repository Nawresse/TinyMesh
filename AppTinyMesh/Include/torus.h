// Torus

#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"

class Torus
{
protected:
  Vector c; //!< Center vertex.
  double ri; //!< Inner radius.
  double ro; //!< Outer radius.
public:
  //! Empty.
  Torus() {}
  explicit Torus(const Vector&, double, double);
  explicit Torus(double, double);
  explicit Torus(const Vector&);

  //! Empty.
  ~Torus() {}

  // Comparison
  friend int operator==(const Torus&, const Torus&);
  friend int operator!=(const Torus&, const Torus&);

  // Acces to properties
  Vector Center() const;
  double InnerRadius() const;
  double OuterRadius() const;

  double Volume() const;
  double Area() const;

  // Translation, scale
  void Translate(const Vector&);
  void Scale(double);

  friend std::ostream& operator<<(std::ostream&, const Torus&);

public:
  static const double epsilon; //!< Internal \htmlonly\epsilon;\endhtmlonly for ray intersection tests.
  static const Torus Null; //!< Empty Torus.
};

//! Returns the center of the torus.
inline Vector Torus::Center() const
{
  return c;
}

/*!
\brief Returns the inner radius of the torus.
*/
inline double Torus::InnerRadius() const
{
  return ri;
}

/*!
\brief Returns the radius of the torus.
*/
inline double Torus::OuterRadius() const
{
  return ro;
}


//! Compute the volume of a Torus.
inline double Torus::Volume() const
{
  return 2 * M_1_PI * M_1_PI * ri * ro * ro;
}

/*!
\brief Compute the surface area of a Torus.
*/
inline double Torus::Area() const
{
  return 4 * M_1_PI * M_1_PI * ri * ro;
}



/*!
\brief Check if two Toruss are (strictly) equal.
\param a, b Toruss.
*/
inline int operator==(const Torus& a, const Torus& b)
{
  return (a.c == b.c) && (a.ri == b.ri) && (a.ro == b.ro);
}

/*!
\brief Check if two Toruss are (strictly) different.
\param a, b Toruss.
*/
inline int operator!=(const Torus& a, const Torus& b)
{
  return !(a == b);
}
