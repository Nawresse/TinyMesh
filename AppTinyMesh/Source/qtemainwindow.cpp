#include "qte.h"
#include "implicits.h"
#include "ui_interface.h"
#include "heightfield.h"
#include <chrono>

MainWindow::MainWindow() : QMainWindow(), uiw(new Ui::Assets)
{
	// Chargement de l'interface
    uiw->setupUi(this);

	// Chargement du GLWidget
	meshWidget = new MeshWidget;
	QGridLayout* GLlayout = new QGridLayout;
	GLlayout->addWidget(meshWidget, 0, 0);
	GLlayout->setContentsMargins(0, 0, 0, 0);
    uiw->widget_GL->setLayout(GLlayout);

	// Creation des connect
	CreateActions();

	meshWidget->SetCamera(Camera(Vector(10, 0, 0), Vector(0.0, 0.0, 0.0)));
}


MainWindow::~MainWindow()
{
	delete meshWidget;
}

void MainWindow::CreateActions()
{
	// Buttons
    connect(uiw->boxMesh, SIGNAL(clicked()), this, SLOT(BoxMeshExample()));

    connect(uiw->sphereMesh, SIGNAL(clicked()), this, SLOT(SphereMesh()));
    connect(uiw->cylinderMesh, SIGNAL(clicked()), this, SLOT(CylinderMesh()));
    connect(uiw->capsuleMesh, SIGNAL(clicked()), this, SLOT(CapsuleMesh()));
    connect(uiw->torusMesh, SIGNAL(clicked()), this, SLOT(TorusMesh()));
    connect(uiw->wrappedMesh, SIGNAL(clicked()), this, SLOT(WrappedMesh()));
    connect(uiw->mergedMesh, SIGNAL(clicked()), this, SLOT(MergedMesh()));

    connect(uiw->terrainMesh, SIGNAL(clicked()), this, SLOT(TerrainMesh()));
    connect(uiw->sphereImplicit, SIGNAL(clicked()), this, SLOT(SphereImplicitExample()));
    connect(uiw->resetcameraButton, SIGNAL(clicked()), this, SLOT(ResetCamera()));

    connect(uiw->wireframe, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_1, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_2, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));

	// Widget edition
	connect(meshWidget, SIGNAL(_signalEditSceneLeft(const Ray&)), this, SLOT(editingSceneLeft(const Ray&)));
	connect(meshWidget, SIGNAL(_signalEditSceneRight(const Ray&)), this, SLOT(editingSceneRight(const Ray&)));
}

void MainWindow::editingSceneLeft(const Ray&)
{
}

void MainWindow::editingSceneRight(const Ray&)
{
}

void MainWindow::BoxMeshExample()
{
	Mesh boxMesh = Mesh(Box(1.0));
    meshColor = MeshColor(boxMesh);
	UpdateGeometry();
}

void MainWindow::SphereMesh()
{
    auto start = std::chrono::high_resolution_clock::now();
    Sphere sphere = Sphere(Vector(1.0),1.0);
    Mesh sphereMesh = Mesh(sphere,31,31);
    meshColor = MeshColor(sphereMesh);
    UpdateGeometry();
    auto stop =std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout<<"The generation time of the sphere is equal to : "<<duration.count()<<"ms\n======"<<std::endl;
}
void MainWindow::CylinderMesh()
{
    auto start = std::chrono::high_resolution_clock::now();
    Mesh cylinderMesh = Mesh(Cylinder(Vector::Null));
    meshColor = MeshColor(cylinderMesh);
    UpdateGeometry();
    auto stop =std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout<<"The generation time of the cylinder is equal to : "<<duration.count()<<"ms\n======"<<std::endl;
}
void MainWindow::CapsuleMesh()
{
    auto start = std::chrono::high_resolution_clock::now();
    Mesh capsuleMesh = Mesh(Capsule(Vector::Null));
    meshColor = MeshColor(capsuleMesh);
    UpdateGeometry();
    auto stop =std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout<<"The generation time of the capsul is equal to : "<<duration.count()<<"ms\n======"<<std::endl;
}
void MainWindow::TorusMesh()
{
    auto start = std::chrono::high_resolution_clock::now();
    Torus torus = Torus(1.0,0.5);
    Mesh torusMesh = Mesh(torus,31,31);
    meshColor = MeshColor(torusMesh);
    UpdateGeometry();
    auto stop =std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout<<"The generation time of the torus is equal to : "<<duration.count()<<"ms\n======"<<std::endl;

}

void MainWindow::MergedMesh()
{
    auto start = std::chrono::high_resolution_clock::now();
    Torus torus = Torus(Vector(0.0) , 2.0, 0.1);
    Mesh torusMesh = Mesh(torus,20,61);

    Sphere sphere = Sphere(1.5);
    Mesh sphereMesh = Mesh(sphere, 31);

    Mesh mergedMesh;
    mergedMesh.Merge(sphereMesh, torusMesh);

    std::vector<Color> cols;
    cols.resize(mergedMesh.Vertexes());
    for (int i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(mergedMesh, cols, mergedMesh.VertexIndexes());
    UpdateGeometry();
    auto stop =std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout<<"The generation time is equal to : "<<duration.count()<<"ms\n======"<<std::endl;


}

void MainWindow::WrappedMesh()
{
    auto start = std::chrono::high_resolution_clock::now();
    Sphere sphere = Sphere(1);
    sphere.Translate(Vector(0.5, 0.0, 0.0));

    Capsule capsule = Capsule(Vector(0.0) , 2.0, 1.0);
    Mesh capsuleMesh = Mesh(capsule, 31);
    capsuleMesh.SphereWrap(sphere,Vector(-1.0));

    std::vector<Color> cols;
    cols.resize(capsuleMesh.Vertexes());
    for (int i = 0; i < cols.size(); i++)
        cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

    meshColor = MeshColor(capsuleMesh, cols, capsuleMesh.VertexIndexes());
    UpdateGeometry();
    auto stop =std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout<<"The generation time is equal to : "<<duration.count()<<"ms\n======"<<std::endl;


}

void MainWindow::SphereImplicitExample()
{
  AnalyticScalarField implicit;

  Mesh implicitMesh;
  implicit.Polygonize(31, implicitMesh, Box(2.0));

  std::vector<Color> cols;
  cols.resize(implicitMesh.Vertexes());
  for (size_t i = 0; i < cols.size(); i++)
    cols[i] = Color(0.8, 0.8, 0.8);

  meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
  UpdateGeometry();
}


void MainWindow::UpdateGeometry()
{
	meshWidget->ClearAll();
	meshWidget->AddMesh("BoxMesh", meshColor);

    uiw->lineEdit->setText(QString::number(meshColor.Vertexes()));
    uiw->lineEdit_2->setText(QString::number(meshColor.Triangles()));

	UpdateMaterial();
}

void MainWindow::UpdateMaterial()
{
    meshWidget->UseWireframeGlobal(uiw->wireframe->isChecked());

    if (uiw->radioShadingButton_1->isChecked())
		meshWidget->SetMaterialGlobal(MeshMaterial::Normal);
	else
		meshWidget->SetMaterialGlobal(MeshMaterial::Color);
}

void MainWindow::ResetCamera()
{
	meshWidget->SetCamera(Camera(Vector(-10.0), Vector(0.0)));
}
void MainWindow::TerrainMesh()
{
    auto start = std::chrono::high_resolution_clock::now();
    std::string image_path = "C:\\Users\\Nawresse\\TinyMesh\\mt-taranaki.png";

    HeightField terrainMesh;
    terrainMesh = HeightField(image_path, 5.0);
    meshColor = MeshColor(terrainMesh, 5.0, true);
    UpdateGeometry();
    auto stop =std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout<<"The generation time is equal to : "<<duration.count()<<"ms\n======"<<std::endl;

}
