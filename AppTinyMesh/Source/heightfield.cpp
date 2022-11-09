#include "heightfield.h"


/*!
\brief Empty
*/
HeightField::HeightField()
{
}



/*!
\brief Empty
*/
HeightField::~HeightField()
{
}



/*!
\brief Creates a terrain mesh from an image.
The object has width*height vertices, width*height normals and  2*(width-1)*(height-1) triangles.
\param image_path The image path.
\param scale The scale of the terrain.
*/
HeightField::HeightField(std::string image_path, double scale)
{
  QImage img(QString::fromStdString(image_path));

  int width = img.width();
  int height = img.height();
  // Reserve space for the triangle array
  varray.reserve(2 * (width - 1) * (height - 1) * 3);
  narray.reserve(2 * (width - 1) * (height - 1) * 3);
  // Reserve space for the vertices and normals
  vertices.reserve(width * height);
  normals.reserve(width * height);

  // Create the vertices
  for (int i=0; i<height; i++)
  {
    for (int j=0; j<width; j++)
    {
      QRgb value_rgb(img.pixel(j, i));
      double value = qGray(value_rgb);
      vertices.push_back(Vector(10*double(j)/(width-1)-5, 10*double(i)/(height-1)-5, value/255.0*scale));
    }
  }

  // Create the normals
  Vector a, b, c, d, n;
  for (int i=0; i<height; i++)
  {
    for (int j=0; j<width; j++)
    {
      if (i == 0)
      {
        a = vertices[i*width+j];
      }
      else
      {
        a = vertices[(i-1)*width+j];
      }
      if (i == height - 1)
      {
        d = vertices[i*width+j];
      }
      else
      {
        d = vertices[(i+1)*width+j];
      }
      if (j == 0)
      {
        b = vertices[i*width+j];
      }
      else
      {
        b = vertices[i*width+j-1];
      }
      if (j == width - 1)
      {
        c = vertices[i*width+j];
      }
      else
      {
        c = vertices[i*width+j+1];
      }
      n = (c - b) / (d - a);
      normals.push_back(n/Norm(n));
    }
  }

  // Create the triangles
  for (int i=0; i<height-1; i++)
  {
    for (int j=0; j<width-1; j++)
    {
      AddSmoothTriangle(i*width+j, i*width+j, i*width+j+1, i*width+j+1, (i+1)*width+j+1, (i+1)*width+j+1);
      AddSmoothTriangle(i*width+j, i*width+j, (i+1)*width+j+1, (i+1)*width+j+1, (i+1)*width+j, (i+1)*width+j);
    }
  }
}
