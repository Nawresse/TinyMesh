// Capsule

// Self include
#include "capsule.h"

/*!
\class Capsule capsule.h
\brief An axis aligned capsule.

The class stores the center, the height and the radius.

\code
Capsule capsule(Vector(0.0,0.0,0.0), 1.0, 1.0); // Unit capsule
\endcode
*/

const double Capsule::epsilon = 1.0e-5; //!< Epsilon value used to check intersections and some round off errors.
const Capsule Capsule::Null(0.0, 0.0); //!< Null capsule, equivalent to: \code Capsule(Vector(0.0)); \endcode


/*!
\brief Create a capsule given a center point, the height and the radius.
\param center Center.
\param radius radius.
*/
Capsule::Capsule(const Vector& c, double h, double r)
{
  Capsule::c = c;
  Capsule::h = abs(h);
  Capsule::r = abs(r);
}



/*!
\brief Create a centered capsule given the radius.

This is equivalent to:
\code
Capsule capsule(Vector(0.0),1.0, 1.0);  // Simplified constructor Capsule(1.0, 1.0);
\endcode
\param r Half side length.
*/
Capsule::Capsule(double h, double r)
{
  Capsule::c = Vector(0.0);
  Capsule::h = abs(h);
  Capsule::r = abs(r);
}


/*!
\brief Create a unit capsule given the center.

This is equivalent to:
\code
Capsule capsule(Vector(1.0),1.0, 1.0);  // Simplified constructor Capsule(Vector(1.0));
\endcode
\param r Half side length.
*/
Capsule::Capsule(const Vector& c)
{
  Capsule::c = c;
  Capsule::h = 1.0;
  Capsule::r = 1.0;
}


/*!
\brief Overloaded.
\param s Stream.
\param capsule The capsule.
*/
std::ostream& operator<<(std::ostream& s, const Capsule& capsule)
{
  s << "Capsule(c:" << capsule.c << ", h:" << capsule.h << ", r:" << capsule.r << ")";
  return s;
}


/*!
\brief Translates a capsule.
\param t Translation vector.
*/
void Capsule::Translate(const Vector& t)
{
  c += t;
}

/*!
\brief Scales a capsule.
\param s Scaling.
*/
void Capsule::Scale(double s)
{
  h *= abs(s);
  r *= abs(s);
}
