// Sphere

// Self include
#include "sphere.h"

/*!
\class Sphere sphere.h
\brief An axis aligned sphere.

The class stores the center and the radius.

\code
Sphere sphere(Vector(0.0,0.0,0.0), 1.0); // Unit sphere
\endcode
*/

const double Sphere::epsilon = 1.0e-5; //!< Epsilon value used to check intersections and some round off errors.
const Sphere Sphere::Null(0.0); //!< Null sphere, equivalent to: \code Sphere(Vector(0.0)); \endcode


/*!
\brief Create a sphere given a center point and the radius.
\param center Center.
\param radius radius.
*/
Sphere::Sphere(const Vector& c, double r)
{
  Sphere::c = c;
  Sphere::r = abs(r);
}



/*!
\brief Create a centered sphere given the radius.

This is equivalent to:
\code
Sphere sphere(Vector(0.0),1.0);  // Simplified constructor Sphere(1.0);
\endcode
\param r Half side length.
*/
Sphere::Sphere(double r)
{
  Sphere::c = Vector(0.0);
  Sphere::r = abs(r);
}


/*!
\brief Create a unit sphere given the center.

This is equivalent to:
\code
Sphere sphere(Vector(1.0),1.0);  // Simplified constructor Sphere(Vector(1.0));
\endcode
\param r Half side length.
*/
Sphere::Sphere(const Vector& c)
{
  Sphere::c = c;
  Sphere::r = 1.0;
}

/*!
\brief Creates the bounding sphere of a set of points.
\param v Array of vertices.
*/
Sphere::Sphere(const std::vector<Vector>& v)
{
  int n = v.size();
  c =  Vector(0.0);
  for (int j = 0; j < n; j++)
  {
    c += v[j];
  }
  c /= n;
  r = 0.0;
  for (int j = 0; j < n; j++)
  {
    double d = Norm(v[j] - c);
    if (d > r)
    {
      r = d;
    }
  }
}



/*!
\brief Overloaded.
\param s Stream.
\param sphere The sphere.
*/
std::ostream& operator<<(std::ostream& s, const Sphere& sphere)
{
  s << "Sphere(" << sphere.c << ", " << sphere.r << ")";
  return s;
}


/*!
\brief Translates a sphere.
\param t Translation vector.
*/
void Sphere::Translate(const Vector& t)
{
  c += t;
}

/*!
\brief Scales a sphere.
\param s Scaling.
*/
void Sphere::Scale(double s)
{
  r *= abs(s);
}

/*!
\brief Wrap from a sphere.
\param vTrans Transformation vector (unitary).
\param v Point to wrap.
*/
void Sphere::Wrap(const Vector& vTrans, const Vector& v, Vector& vTemp) const
{
  Vector w = v - c;
  double d = Norm(w);
  if (d >= r)
  {
    vTemp = v;
  }
  else
  {
    vTemp = v + vTrans * (1 - d / r);
  }
}

