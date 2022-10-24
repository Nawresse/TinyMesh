#include "mesh.h"

/*!
\class Mesh mesh.h

\brief Core triangle mesh class.
*/



/*!
\brief Initialize the mesh to empty.
*/
Mesh::Mesh()
{
}

/*!
\brief Initialize the mesh from a list of vertices and a list of triangles.

Indices must have a size multiple of three (three for triangle vertices and three for triangle normals).

\param vertices List of geometry vertices.
\param indices List of indices wich represent the geometry triangles.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<int>& indices) :vertices(vertices), varray(indices)
{
  normals.resize(vertices.size(), Vector::Z);
}

/*!
\brief Create the mesh.

\param vertices Array of vertices.
\param normals Array of normals.
\param va, na Array of vertex and normal indexes.
*/
Mesh::Mesh(const std::vector<Vector>& vertices, const std::vector<Vector>& normals, const std::vector<int>& va, const std::vector<int>& na) :vertices(vertices), normals(normals), varray(va), narray(na)
{
}

/*!
\brief Reserve memory for arrays.
\param nv,nn,nvi,nvn Number of vertices, normals, vertex indexes and vertex normals.
*/
void Mesh::Reserve(int nv, int nn, int nvi, int nvn)
{
  vertices.reserve(nv);
  normals.reserve(nn);
  varray.reserve(nvi);
  narray.reserve(nvn);
}

/*!
\brief Empty
*/
Mesh::~Mesh()
{
}

/*!
\brief Smooth the normals of the mesh.

This function weights the normals of the faces by their corresponding area.
\sa Triangle::AreaNormal()
*/
void Mesh::SmoothNormals()
{
  // Initialize
  normals.resize(vertices.size(), Vector::Null);

  narray = varray;

  // Accumulate normals
  for (int i = 0; i < varray.size(); i += 3)
  {
    Vector tn = Triangle(vertices[varray.at(i)], vertices[varray.at(i + 1)], vertices[varray.at(i + 2)]).AreaNormal();
    normals[narray[i + 0]] += tn;
    normals[narray[i + 1]] += tn;
    normals[narray[i + 2]] += tn;
  }

  // Normalize
  for (int i = 0; i < normals.size(); i++)
  {
    Normalize(normals[i]);
  }
}

/*!
\brief Add a smooth triangle to the geometry.
\param a, b, c Index of the vertices.
\param na, nb, nc Index of the normals.
*/
void Mesh::AddSmoothTriangle(int a, int na, int b, int nb, int c, int nc)
{
  varray.push_back(a);
  narray.push_back(na);
  varray.push_back(b);
  narray.push_back(nb);
  varray.push_back(c);
  narray.push_back(nc);
}

/*!
\brief Add a triangle to the geometry.
\param a, b, c Index of the vertices.
\param n Index of the normal.
*/
void Mesh::AddTriangle(int a, int b, int c, int n)
{
  varray.push_back(a); // pourquoi ici on push un sommet 1 normal, 1sommet 1normal
  narray.push_back(n);
  varray.push_back(b);
  narray.push_back(n);
  varray.push_back(c);
  narray.push_back(n);
}

/*!
\brief Add a smmoth quadrangle to the geometry.

Creates two smooth triangles abc and acd.

\param a, b, c, d  Index of the vertices.
\param na, nb, nc, nd Index of the normal for all vertices.
*/
void Mesh::AddSmoothQuadrangle(int a, int na, int b, int nb, int c, int nc, int d, int nd)
{
  // First triangle
  AddSmoothTriangle(a, na, b, nb, c, nc);

  // Second triangle
  AddSmoothTriangle(a, na, c, nc, d, nd);
}

/*!
\brief Add a quadrangle to the geometry.

\param a, b, c, d  Index of the vertices and normals.
*/
void Mesh::AddQuadrangle(int a, int b, int c, int d)
{
  AddSmoothQuadrangle(a, a, b, b, c, c, d, d);
}

/*!
\brief Compute the bounding box of the object.
*/
Box Mesh::GetBox() const
{
  if (vertices.size() == 0)
  {
    return Box::Null;
  }
  return Box(vertices);
}

/*!
\brief Creates an axis aligned box.

The object has 8 vertices, 6 normals and 12 triangles.
\param box The box.
*/
Mesh::Mesh(const Box& box)
{
  // Vertices
  vertices.resize(8);

  for (int i = 0; i < 8; i++)
  {
    vertices[i] = box.Vertex(i);
  }

  // Normals
  normals.push_back(Vector(-1, 0, 0));
  normals.push_back(Vector(1, 0, 0));
  normals.push_back(Vector(0, -1, 0));
  normals.push_back(Vector(0, 1, 0));
  normals.push_back(Vector(0, 0, -1));
  normals.push_back(Vector(0, 0, 1));

  // Reserve space for the triangle array
  varray.reserve(12 * 3);
  narray.reserve(12 * 3);

  AddTriangle(0, 2, 1, 4);
  AddTriangle(1, 2, 3, 4);

  AddTriangle(4, 5, 6, 5);
  AddTriangle(5, 7, 6, 5);

  AddTriangle(0, 4, 2, 0);
  AddTriangle(4, 6, 2, 0);

  AddTriangle(1, 3, 5, 1);
  AddTriangle(3, 7, 5, 1);

  AddTriangle(0, 1, 5, 2);
  AddTriangle(0, 5, 4, 2);

  AddTriangle(3, 2, 7, 3);
  AddTriangle(6, 7, 2, 3);
}

/*!
\brief Creates an axis aligned sphere.

The object has (n-2)*n+2 vertices, n*(n-1) normals and 2*n*(n-2) triangles.
\param sphere The sphere.
\param nc The number of circles.
\param npc The number of subdivisions per circle.
*/
Mesh::Mesh(const Sphere& sphere, int nc, int npc)
{
if (nc < 3) // Pour éviter une boucle vide , évider une discrétisation qui n'a pas de sens
 {
    nc = 3;
 }
if (npc <3){

    npc = 3;
}
  double r = sphere.Radius(); // on declare le rayon de la sphere
  Vector c = sphere.Center(); // ici on déclare le centre de la sphere
  // Reserve space for the triangle array
  varray.reserve(2 * npc * (nc - 2) * 3);
  narray.reserve(2 * npc * (nc - 2) * 3);

  // Create the vertices
  vertices.push_back(c + Vector(0, 0, r)); // pôle nord
  for (int i=1; i<nc-1; i++)
  {
    double theta = M_PI * i / (nc-1);
    double rho = r * sin(theta); // projection dans le plan x y
    double z = r * cos(theta); // projection sur l'axe z
    for (int j=0; j<npc; j++)
    {
      double phi = 2 * M_PI * j / npc;
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      vertices.push_back(c + Vector(x, y, z));
    }
  }
  vertices.push_back(c + Vector(0, 0, -r)); // pôle sud

  // Create the triangles and normals
  Triangle t; // Temporary triangle to compute the normals
  for (int j=0; j<npc; j++)
  {
    t = Triangle(vertices[0], vertices[j+1], vertices[(j+1)%npc+1]);
    normals.push_back(t.Normal());
    AddTriangle(0, j+1, (j+1)%npc+1, j);
  }
  for (int i=1; i<nc-2; i++)
  {
    for (int j=0; j<npc; j++)
    {
      t = Triangle(vertices[(i-1)*npc+j+1], vertices[i*npc+j+1], vertices[i*npc+(j+1)%npc+1]);
      normals.push_back(t.Normal());
      AddTriangle((i-1)*npc+j+1, i*npc+j+1, i*npc+(j+1)%npc+1, (i-1)*npc+j+npc);
      AddTriangle((i-1)*npc+j+1, i*npc+(j+1)%npc+1, (i-1)*npc+(j+1)%npc+1, (i-1)*npc+j+npc);
    }
  }
  for (int j=0; j<npc; j++)
  {
    t = Triangle(vertices[(nc-3)*npc+j+1], vertices[(nc-2)*npc+1], vertices[(nc-3)*npc+(j+1)%npc+1]);
    normals.push_back(t.Normal());
    AddTriangle((nc-3)*npc+j+1, (nc-2)*npc+1, (nc-3)*npc+(j+1)%npc+1, (nc-3)*npc+j+npc);
  }

}




/*!
\brief Creates an axis aligned cylinder.

The object has (nh-2)*nr+2 vertices, nr*(nh-1) normals and 2*nr*(nh-2) triangles.
\param cylinder The cylinder.
\param nh The number of vertical subdivisions.
\param nr The number of radial subdivisions.
*/
Mesh::Mesh(const Cylinder& cylinder, int nh, int nr)
{
  if (nr == -1)
  {
    nr = 2 * nh;
  }
  if (nh < 3)
  {
    nh = 3;
  }
  if (nr < 3)
  {
    nr = 3;
  }
  Vector c = cylinder.Center();
  double h = cylinder.Height();
  double r = cylinder.Radius();
  // Reserve space for the triangle array
  varray.reserve(2 * nr * (nh - 2) * 3);
  narray.reserve(2 * nr * (nh - 2) * 3);

  // Create the vertices
  vertices.push_back(c + Vector(0, 0, h/2));
  for (int i=0; i<nh-2; i++)
  {
    double z = h/2 - h * i / (nh-3) ;
    for (int j=0; j<nr; j++)
    {
      double theta = 2 * M_PI * j / nr;
      double x = r * cos(theta);
      double y = r * sin(theta);
      vertices.push_back(c + Vector(x, y, z));
    }
  }
  vertices.push_back(c + Vector(0, 0, -h/2));

  // Create the triangles and normals
  Triangle t; // Temporary triangle to compute the normals
  for (int j=0; j<nr; j++)
  {
    t = Triangle(vertices[0], vertices[j+1], vertices[(j+1)%nr+1]);
    normals.push_back(t.Normal());
    AddTriangle(0, j+1, (j+1)%nr+1, j);
  }
  for (int i=1; i<nh-2; i++)
  {
    for (int j=0; j<nr; j++)
    {
      t = Triangle(vertices[(i-1)*nr+j+1], vertices[i*nr+j+1], vertices[i*nr+(j+1)%nr+1]);
      normals.push_back(t.Normal());
      AddTriangle((i-1)*nr+j+1, i*nr+j+1, i*nr+(j+1)%nr+1, (i-1)*nr+j+nr);
      AddTriangle((i-1)*nr+j+1, i*nr+(j+1)%nr+1, (i-1)*nr+(j+1)%nr+1, (i-1)*nr+j+nr);
    }
  }
  for (int j=0; j<nr; j++)
  {
    t = Triangle(vertices[(nh-3)*nr+j+1], vertices[(nh-2)*nr+1], vertices[(nh-3)*nr+(j+1)%nr+1]);
    normals.push_back(t.Normal());
    AddTriangle((nh-3)*nr+j+1, (nh-2)*nr+1, (nh-3)*nr+(j+1)%nr+1, (nh-3)*nr+j+nr);
  }

}



/*!
\brief Creates an axis aligned capsule.

\param capsule The capsule.
\param nh The number of vertical subdivisions.
\param nr The number of radial subdivisions.
*/
Mesh::Mesh(const Capsule& capsule, int nh, int nr)
{
  if (nr == -1)
  {
    nr = 2 * nh;
  }
  if (nh < 3)
  {
    nh = 3;
  }
  if (nr < 3)
  {
    nr = 3;
  }
  if (nr % 2 == 0)
  {
    nr++;
  }

  Vector c = capsule.Center();
  double h = capsule.Height();
  double r = capsule.Radius();
  // Reserve space for the triangle array
  varray.reserve(2 * nr * (nr + nh - 5) * 3);
  narray.reserve(2 * nr * (nr + nh - 5) * 3);

  // Create the vertices of the first half sphere
  vertices.push_back(c + Vector(0, 0, r + h/2)); // north pole of the first half sphere
  for (int i=1; i<nr/2; i++)
  {
    double theta = M_PI * i / (nr-1);
    double rho = r * sin(theta);
    double z = r * cos(theta);
    for (int j=0; j<nr; j++)
    {
      double phi = 2 * M_PI * j / nr;
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      vertices.push_back(c + Vector(0, 0, h/2) + Vector(x, y, z));
    }
  }
  // Create the vertices of the inner cylinder
  for (int i=0; i<nh-2; i++)
  {
    double z = h/2 - h * i / (nh-3) ;
    for (int j=0; j<nr; j++)
    {
      double theta = 2 * M_PI * j / nr;
      double x = r * cos(theta);
      double y = r * sin(theta);
      vertices.push_back(c + Vector(x, y, z));
    }
  }
  // Create the vertices of the second half sphere
  for (int i=nr/2+1; i<nr-1; i++)
  {
    double theta = M_PI * i / (nr - 1);
    double rho = r * sin(theta);
    double z = r * cos(theta);
    for (int j=0; j<nr; j++)
    {
      double phi = 2 * M_PI * j / nr;
      double x = rho * cos(phi);
      double y = rho * sin(phi);
      vertices.push_back(c + Vector(0, 0, -h/2) + Vector(x, y, z)); // translation
    }
  }
  vertices.push_back(c + Vector(0, 0, -r - h/2)); // south pole of the second half sphere


  // Create the triangles and normals
  Triangle t; // Temporary triangle to compute the normals
  for (int j=0; j<nr; j++)
  {
    t = Triangle(vertices[0], vertices[j+1], vertices[(j+1)%nr+1]);
    normals.push_back(t.Normal());
    AddTriangle(0, j+1, (j+1)%nr+1, j);
  }
  for (int i=1; i<nr+nh-5; i++)
  {
    for (int j=0; j<nr; j++)
    {
      t = Triangle(vertices[(i-1)*nr+j+1], vertices[i*nr+j+1], vertices[i*nr+(j+1)%nr+1]);
      normals.push_back(t.Normal());
      AddTriangle((i-1)*nr+j+1, i*nr+j+1, i*nr+(j+1)%nr+1, (i-1)*nr+j+nr);
      AddTriangle((i-1)*nr+j+1, i*nr+(j+1)%nr+1, (i-1)*nr+(j+1)%nr+1, (i-1)*nr+j+nr);
    }
  }
  for (int j=0; j<nr; j++)
  {
    t = Triangle(vertices[(nr+nh-6)*nr+j+1], vertices[(nr+nh-5)*nr+1], vertices[(nr+nh-6)*nr+(j+1)%nr+1]);
    normals.push_back(t.Normal());
    AddTriangle((nr+nh-6)*nr+j+1, (nr+nh-5)*nr+1, (nr+nh-6)*nr+(j+1)%nr+1, (nr+nh-6)*nr+j+nr);
  }

}


/*!
\brief Creates an axis aligned tore.
The object has n*n vertices, n*n normals and 2*n*n triangles.
\param tore The tore.
\param n The number of subdivisions.
*/
Mesh::Mesh(const Torus& torus, int nc, int npc)
{
  if (nc < 3)
  {
    nc = 3;
  }
  if (npc < 3)
  {
    npc = 3;
  }
  double ri = torus.InnerRadius();
  double ro = torus.OuterRadius();
  Vector c = torus.Center();
  // Reserve space for the triangle array
  varray.reserve(2 * nc * npc* 3);
  narray.reserve(2 * nc * npc * 3);

  // Create the vertices
  for (int i=0; i<nc; i++)
  {
    double u = 2 * M_PI * i / nc;
    double z = ro * sin(u);
    for (int j=0; j<npc; j++)
    {
      double v = 2 * M_PI * j / npc;
      double y = (ri - ro * cos(u)) * sin(v);
      double x = (ri - ro * cos(u)) * cos(v);
      vertices.push_back(c + Vector(x, y, z));
    }
  }

  // Create the triangles and normals
  Triangle t; // Temporary triangle to compute the normals
  for (int i=0; i<nc; i++)
  {
    for (int j=0; j<npc; j++)
    {
      t = Triangle(vertices[i*npc+j], vertices[((i+1)%nc)*npc+j], vertices[((i+1)%nc)*npc+(j+1)%npc]);
      normals.push_back(t.Normal());
      AddTriangle(i*npc+j, ((i+1)%nc)*npc+j, ((i+1)%nc)*npc+(j+1)%npc, i*npc+j); // ajoute un triangle
      AddTriangle(i*npc+j, ((i+1)%nc)*npc+(j+1)%npc, i*npc+(j+1)%npc, i*npc+j);
    }
  }
}



/*!
\brief Scale the mesh.
\param s Scaling factor.
*/
void Mesh::Scale(double s)
{
  // Vertexes
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] *= s;
  }

  if (s < 0.0)
  {
    // Normals
    for (int i = 0; i < normals.size(); i++)
    {
      normals[i] = -normals[i];
    }
  }
}

void Mesh::Translate(const Vector& v)
{
  // Vertexes
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] += v;
  }
}

void Mesh::Rotate(const Vector& axis, double angle)
{
  Matrix rot = Matrix(axis, angle);
  // Vertexes
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] = rot * vertices[i];
  }

  // Normals
  for (int i = 0; i < normals.size(); i++)
  {
    normals[i] = rot * normals[i];
  }
}


void Mesh::Rotate(const Matrix& rot)
{
  // Vertexes
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] = rot * vertices[i];
  }

  // Normals
  for (int i = 0; i < normals.size(); i++)
  {
    normals[i] = rot * normals[i];
  }
}

void Mesh::SphereWrap(const Sphere& sphere)
{
  // Vertexes
  for (int i = 0; i < vertices.size(); i++)
  {
    vertices[i] = sphere.Wrap(vertices[i]);
  }
}

void Mesh::Merge(const Mesh& mesh1, const Mesh& mesh2)
{
  int nv1 = mesh1.vertices.size();
  int nv2 = mesh2.vertices.size();
  int nn1 = mesh1.normals.size();
  int nn2 = mesh2.normals.size();

  // Vertexes
  vertices.reserve(nv1 + nv2);
  vertices.insert(vertices.end(), mesh1.vertices.begin(), mesh1.vertices.end());
  vertices.insert(vertices.end(), mesh2.vertices.begin(), mesh2.vertices.end());

  // Normals
  normals.reserve(nn1 + nn2);
  normals.insert(normals.end(), mesh1.normals.begin(), mesh1.normals.end());
  normals.insert(normals.end(), mesh2.normals.begin(), mesh2.normals.end());

  // Triangles
  varray.reserve(mesh1.varray.size() + mesh2.varray.size());
  varray.insert(varray.end(), mesh1.varray.begin(), mesh1.varray.end());
  narray.reserve(mesh1.narray.size() + mesh2.narray.size());
  narray.insert(narray.end(), mesh1.narray.begin(), mesh1.narray.end());

  for (int i = 0; i < mesh2.varray.size(); i++)
  {
    varray.push_back(mesh2.varray[i] + nv1);
  }
  for (int i = 0; i < mesh2.narray.size(); i++)
  {
    narray.push_back(mesh2.narray[i] + nn1);
  }
}


#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtCore/QRegularExpression>
#include <QtCore/qstring.h>

/*!
\brief Import a mesh from an .obj file.
\param filename File name.
*/
void Mesh::Load(const QString& filename)
{
  vertices.clear();
  normals.clear();
  varray.clear();
  narray.clear();

  QFile data(filename);

  if (!data.open(QFile::ReadOnly))
    return;
  QTextStream in(&data);

  // Set of regular expressions : Vertex, Normal, Triangle
  QRegularExpression rexv("v\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
  QRegularExpression rexn("vn\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)\\s*([-|+|\\s]\\d*\\.\\d+)");
  QRegularExpression rext("f\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)\\s*(\\d*)/\\d*/(\\d*)");
  while (!in.atEnd())
  {
    QString line = in.readLine();
    QRegularExpressionMatch match = rexv.match(line);
    QRegularExpressionMatch matchN = rexn.match(line);
    QRegularExpressionMatch matchT = rext.match(line);
    if (match.hasMatch())//rexv.indexIn(line, 0) > -1)
    {
      Vector q = Vector(match.captured(1).toDouble(), match.captured(2).toDouble(), match.captured(3).toDouble()); vertices.push_back(q);
    }
    else if (matchN.hasMatch())//rexn.indexIn(line, 0) > -1)
    {
      Vector q = Vector(matchN.captured(1).toDouble(), matchN.captured(2).toDouble(), matchN.captured(3).toDouble());  normals.push_back(q);
    }
    else if (matchT.hasMatch())//rext.indexIn(line, 0) > -1)
    {
      varray.push_back(matchT.captured(1).toInt() - 1);
      varray.push_back(matchT.captured(3).toInt() - 1);
      varray.push_back(matchT.captured(5).toInt() - 1);
      narray.push_back(matchT.captured(2).toInt() - 1);
      narray.push_back(matchT.captured(4).toInt() - 1);
      narray.push_back(matchT.captured(6).toInt() - 1);
    }
  }
  data.close();
}

/*!
\brief Save the mesh in .obj format, with vertices and normals.
\param url Filename.
\param meshName %Mesh name in .obj file.
*/
void Mesh::SaveObj(const QString& url, const QString& meshName) const
{
  QFile data(url);
  if (!data.open(QFile::WriteOnly))
    return;
  QTextStream out(&data);
  out << "g " << meshName << Qt::endl;
  for (int i = 0; i < vertices.size(); i++)
    out << "v " << vertices.at(i)[0] << " " << vertices.at(i)[1] << " " << vertices.at(i)[2] << QString('\n');
  for (int i = 0; i < normals.size(); i++)
    out << "vn " << normals.at(i)[0] << " " << normals.at(i)[1] << " " << normals.at(i)[2] << QString('\n');
  for (int i = 0; i < varray.size(); i += 3)
  {
    out << "f " << varray.at(i) + 1 << "//" << narray.at(i) + 1 << " "
      << varray.at(i + 1) + 1 << "//" << narray.at(i + 1) + 1 << " "
      << varray.at(i + 2) + 1 << "//" << narray.at(i + 2) + 1 << " "
      << "\n";
  }
  out.flush();
  data.close();
}

