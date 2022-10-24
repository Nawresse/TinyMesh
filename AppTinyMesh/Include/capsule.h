// Capsule

#pragma once

#include <vector>
#include <iostream>

#include "mathematics.h"

class Capsule
{
protected:
  Vector c; //!< Center vertex.
  double h; //!< Height.
  double r; //!< Radius.
public:
  //! Empty.
  Capsule() {}
  explicit Capsule(const Vector&, double, double);
  explicit Capsule(double, double);
  explicit Capsule(const Vector&);

  //! Empty.
  ~Capsule() {}

  // Comparison
  friend int operator==(const Capsule&, const Capsule&);
  friend int operator!=(const Capsule&, const Capsule&);

  // Acces to properties
  Vector Center() const;
  double Height() const;
  double Radius() const;

  double Volume() const;
  double Area() const;

  // Translation, scale
  void Translate(const Vector&);
  void Scale(double);

  friend std::ostream& operator<<(std::ostream&, const Capsule&);

public:
  static const double epsilon; //!< Internal \htmlonly\epsilon;\endhtmlonly for ray intersection tests.
  static const Capsule Null; //!< Empty Capsule.
};

//! Returns the center of the capsule.
inline Vector Capsule::Center() const
{
  return c;
}

/*!
\brief Returns the height of the capsule.
*/
inline double Capsule::Height() const
{
  return h;
}

/*!
\brief Returns the radius of the capsule.
*/
inline double Capsule::Radius() const
{
  return r;
}


//! Compute the volume of a Capsule.
inline double Capsule::Volume() const
{
  return M_1_PI * r * r * h + 4.0 / 3.0 * M_1_PI * r * r * r;
}

/*!
\brief Compute the surface area of a Capsule.
*/
inline double Capsule::Area() const
{
  return 2.0 * M_1_PI * r * (h + 2.0 * r);
}



/*!
\brief Check if two Capsules are (strictly) equal.
\param a, b Capsules.
*/
inline int operator==(const Capsule& a, const Capsule& b)
{
  return (a.c == b.c) && (a.r == b.r) && (a.h == b.h);
}

/*!
\brief Check if two Capsules are (strictly) different.
\param a, b Capsules.
*/
inline int operator!=(const Capsule& a, const Capsule& b)
{
  return !(a == b);
}
