#ifndef SPHEREWRAPPED_H
#define SPHEREWRAPPED_H

#endif // SPHEREWRAPPED_H
#include "sphere.h"

class SphereWrapped
{
protected:
  Sphere s ;// Sphere.
  Vector v ; // Vector of deformation.

public :
  SphereWrapped(){};
  explicit SphereWrapped(Sphere s, const Vector&);

  //Deformation
  void SphereDeformation(Sphere s, const Vector&);
};
