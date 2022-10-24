// Cylinder

// Self include
#include "cylinder.h"

/*!
\class Cylinder cylinder.h
\brief An axis aligned cylinder.

The class stores the center, the height and the radius.

\code
Cylinder cylinder(Vector(0.0,0.0,0.0), 1.0, 1.0); // Unit cylinder
\endcode
*/

const double Cylinder::epsilon = 1.0e-5; //!< Epsilon value used to check intersections and some round off errors.
const Cylinder Cylinder::Null(0.0, 0.0); //!< Null cylinder, equivalent to: \code Cylinder(Vector(0.0)); \endcode


/*!
\brief Create a cylinder given a center point, the height and the radius.
\param center Center.
\param radius radius.
*/
Cylinder::Cylinder(const Vector& c, double h, double r)
{
  Cylinder::c = c;
  Cylinder::h = abs(h);
  Cylinder::r = abs(r);
}



/*!
\brief Create a centered cylinder given the radius.

This is equivalent to:
\code
Cylinder cylinder(Vector(0.0),1.0, 1.0);  // Simplified constructor Cylinder(1.0, 1.0);
\endcode
\param r Half side length.
*/
Cylinder::Cylinder(double h, double r)
{
  Cylinder::c = Vector(0.0);
  Cylinder::h = abs(h);
  Cylinder::r = abs(r);
}


/*!
\brief Create a unit cylinder given the center.

This is equivalent to:
\code
Cylinder cylinder(Vector(1.0),1.0, 1.0);  // Simplified constructor Cylinder(Vector(1.0));
\endcode
\param r Half side length.
*/
Cylinder::Cylinder(const Vector& c)
{
  Cylinder::c = c;
  Cylinder::h = 1.0;
  Cylinder::r = 1.0;
}


/*!
\brief Overloaded.
\param s Stream.
\param cylinder The cylinder.
*/
std::ostream& operator<<(std::ostream& s, const Cylinder& cylinder)
{
  s << "Cylinder(c:" << cylinder.c << ", h:" << cylinder.h << ", r:" << cylinder.r << ")";
  return s;
}


/*!
\brief Translates a cylinder.
\param t Translation vector.
*/
void Cylinder::Translate(const Vector& t)
{
  c += t;
}

/*!
\brief Scales a cylinder.
\param s Scaling.
*/
void Cylinder::Scale(double s)
{
  h *= abs(s);
  r *= abs(s);
}
