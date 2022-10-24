// Torus

// Self include
#include "torus.h"

/*!
\class Torus torus.h
\brief An axis aligned torus.

The class storuss the center, the height and the radius.

\code
Torus torus(Vector(0.0,0.0,0.0), 1.0, 1.0); // Unit torus
\endcode
*/

const double Torus::epsilon = 1.0e-5; //!< Epsilon value used to check intersections and some round off errors.
const Torus Torus::Null(0.0, 0.0); //!< Null torus, equivalent to: \code Torus(Vector(0.0)); \endcode


/*!
\brief Create a torus given a center point, the height and the radius.
\param center Center.
\param radius radius.
*/
Torus::Torus(const Vector& c, double ri, double ro)
{
  Torus::c = c;
  Torus::ri = abs(ri);
  Torus::ro = abs(ro);
}



/*!
\brief Create a centered torus given the radius.

This is equivalent to:
\code
Torus torus(Vector(0.0),1.0, 1.0);  // Simplified constructor Torus(1.0, 1.0);
\endcode
\param r Half side length.
*/
Torus::Torus(double h, double r)
{
  Torus::c = Vector(0.0);
  Torus::ri = abs(ri);
  Torus::ro = abs(ro);
}


/*!
\brief Create a unit torus given the center.

This is equivalent to:
\code
Torus torus(Vector(1.0),1.0, 1.0);  // Simplified constructor Torus(Vector(1.0));
\endcode
\param r Half side length.
*/
Torus::Torus(const Vector& c)
{
  Torus::c = c;
  Torus::ri = 1.0;
  Torus::ro = 1.0;
}


/*!
\brief Overloaded.
\param s Stream.
\param torus The torus.
*/
std::ostream& operator<<(std::ostream& s, const Torus& torus)
{
  s << "Torus(c:" << torus.c << ", ri:" << torus.ri << ", ro:" << torus.ro << ")";
  return s;
}


/*!
\brief Translates a torus.
\param t Translation vector.
*/
void Torus::Translate(const Vector& t)
{
  c += t;
}

/*!
\brief Scales a torus.
\param s Scaling.
*/
void Torus::Scale(double s)
{
  ri *= abs(s);
  ro *= abs(s);
}
