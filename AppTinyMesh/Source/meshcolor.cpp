#include "meshcolor.h"

/*!
\brief Create an empty mesh.
*/
MeshColor::MeshColor()
{
}

/*!
\brief Constructor from a Mesh with color array and indices.
\param m Base mesh.
\param cols Color array.
\param carr Color indexes, should be the same size as Mesh::varray and Mesh::narray.
*/
MeshColor::MeshColor(const Mesh& m, const std::vector<Color>& cols, const std::vector<int>& carr) : Mesh(m), colors(cols), carray(carr)
{
}

/*!
\brief Constructor from a Mesh.
\param m the base mesh
*/
MeshColor::MeshColor(const Mesh& m) : Mesh(m)
{
	colors.resize(vertices.size(), Color(1.0, 1.0, 1.0));
	carray = varray;
}

/*!
\brief Empty.
*/
MeshColor::~MeshColor()
{
}
/*!
\brief Creates the terrain colors from a mesh.
\param m The mesh.
\param scale The scale of the terrain.
*/
MeshColor::MeshColor(const Mesh& m, double scale, bool flattenSea) : Mesh(m)
{
  // Defining landscape levels
  double oceanLevel = 0.0;
  double coastLevel = 0.3;
  double beachLevel = 0.4;
  double plainLevel = 0.5;
  double forestLevel = 0.7;
  double hillLevel = 0.8;
  double mountainLevel = 1.0;

  // Defining landscape colors
  Color oceanColor = Color(0, 0, 179);
  Color coastColor = Color(0, 204, 255);
  Color beachColor = Color(255, 230, 179);
  Color plainColor = Color(57, 172, 57);
  Color forestColor = Color(0, 77, 0);
  Color hillColor = Color(89, 89, 89);
  Color mountainColor = Color(242, 242, 242);

  // Reserving memory
  colors.reserve(vertices.size());
  carray = varray;

  // Defining the color array
  Color vertexColor;
  for (int i=0; i<vertices.size(); i++)
  {
    double level = vertices[i] * Vector::Z / scale; // Between 0 and 1
    if (level < coastLevel)
    {
      vertexColor = Color::Lerp((level - oceanLevel) / (coastLevel - oceanLevel), oceanColor, coastColor);
      if (flattenSea)
      {
        vertices[i] = Vector(vertices[i]*Vector::X, vertices[i]*Vector::Y, beachLevel * scale);
      }
    }
    else if (level < beachLevel)
    {
      vertexColor = Color::Lerp((level - coastLevel) / (beachLevel - coastLevel), coastColor, beachColor);
      if (flattenSea)
      {
        vertices[i] = Vector(vertices[i]*Vector::X, vertices[i]*Vector::Y, beachLevel * scale);
      }
    }
    else if (level < plainLevel)
    {
      vertexColor = Color::Lerp((level - beachLevel) / (plainLevel - beachLevel), beachColor, plainColor);
    }
    else if (level < forestLevel)
    {
      vertexColor = Color::Lerp((level - plainLevel) / (forestLevel - plainLevel), plainColor, forestColor);
    }
    else if (level < hillLevel)
    {
      vertexColor = Color::Lerp((level - forestLevel) / (hillLevel - forestLevel), forestColor, hillColor);
    }
    else
    {
      vertexColor = Color::Lerp((level - hillLevel) / (mountainLevel - hillLevel), hillColor, mountainColor);
    }
    colors.push_back(vertexColor);
  }
}
