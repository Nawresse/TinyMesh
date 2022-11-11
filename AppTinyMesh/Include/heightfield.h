// Heightfield class

#pragma once
#include <QImage>
#include <QString>

#include "mesh.h"

class HeightField: public Mesh
{
public:
    // Empty
    explicit HeightField();
    // Constructors from image
    explicit HeightField(std::string, double);
    ~HeightField();
    //getteur

};
