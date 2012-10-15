#include "mainwindow.h"
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderer.h>
#include <QVTKInteractor.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkGlyph3D.h>
#include <vtkPoints.h>
#include <vtkPointSet.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkRegularPolygonSource.h>    
#include "MatrixCode.h"
#include "ExpLogTransform.h"

using namespace std;

typedef vtkSmartPointer<vtkPolyData> vtkPolyDataPointer;
typedef vtkSmartPointer<vtkPoints> vtkPointsPointer;
typedef vtkSmartPointer<vtkPointSet> vtkPointSetPointer;
typedef vtkSmartPointer<vtkSphereSource> vtkSphereSourcePointer;
typedef vtkSmartPointer<vtkGlyph3D> vtkGlyph3DPointer;
typedef vtkSmartPointer<vtkPolyDataMapper> vtkPolyDataMapperPointer;
typedef vtkSmartPointer<vtkActor> vtkActorPointer;
typedef vtkSmartPointer<vtkProperty> vtkPropertyPointer;

vtkPolyDataPointer CreateCircle(PNSBase::VectorType normal, double phi) {
    normal.print("Normal: ");
    normal /= arma::norm(normal, 2);
    double r = sin(phi) + .1;

    PNSBase::VectorType center(3);
    center = cos(phi) * normal;
    center.print("Center: ");
    vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
    polygonSource->SetNumberOfSides(64);
    polygonSource->SetRadius(r);
    polygonSource->SetCenter(center.colptr(0));
    polygonSource->SetNormal(normal[0],normal[1],normal[2]);
    polygonSource->Update();
    return vtkPolyDataPointer(polygonSource->GetOutput());
}


vtkPolyDataPointer CreateCircle(double nx, double ny, double nz, double cx, double cy, double cz, double r) {
    vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
    polygonSource->SetNumberOfSides(64);
    polygonSource->SetRadius(r);
    polygonSource->SetCenter(cx,cy,cz);
    polygonSource->SetNormal(nx,ny,nz);
    polygonSource->Update();
    return vtkPolyDataPointer(polygonSource->GetOutput());
}

vtkPolyDataPointer CreateSinglePoint(PNSBase::VectorType& point, double radius = 0.05) {
    vtkPolyDataPointer pointSet = vtkPolyDataPointer::New();
    vtkPointsPointer points = vtkPointsPointer::New();
    points->SetNumberOfPoints(1);
    points->SetPoint(0, point[0], point[1], point[2]);
    pointSet->SetPoints(points);
    vtkGlyph3DPointer glyph = vtkGlyph3DPointer::New();
    vtkSphereSourcePointer sphereSource = vtkSphereSourcePointer::New();
    sphereSource->SetThetaResolution(64);
    sphereSource->SetPhiResolution(64);
    sphereSource->SetRadius(radius);
    glyph->SetSourceConnection(sphereSource->GetOutputPort());
    glyph->SetInputData(pointSet);
    glyph->Update();
    return vtkPolyDataPointer(glyph->GetOutput());
}

vtkPolyDataPointer CreatePoints(PNSBase::MatrixType& randomPoints, double radius) {
    vtkPolyDataPointer pointSet = vtkPolyDataPointer::New();
    vtkPointsPointer points = vtkPointsPointer::New();
    points->SetNumberOfPoints(randomPoints.n_cols);
    for (int i = 0; i < randomPoints.n_cols; i++) {
        points->SetPoint(i, randomPoints.at(0,i), randomPoints.at(1,i), randomPoints.at(2,i));
    }
    pointSet->SetPoints(points);
    vtkGlyph3DPointer glyph = vtkGlyph3DPointer::New();
    vtkSphereSourcePointer sphereSource = vtkSphereSourcePointer::New();
    sphereSource->SetThetaResolution(64);
    sphereSource->SetPhiResolution(64);
    sphereSource->SetRadius(radius);
    glyph->SetSourceConnection(sphereSource->GetOutputPort());
    glyph->SetInputData(pointSet);
    glyph->Update();
    return vtkPolyDataPointer(glyph->GetOutput());
}

MainWindow::MainWindow(QWidget* parent): QMainWindow(parent) {
    ui.setupUi(this);
    m_Renderer = vtkRenderer::New();
    vtkGenericOpenGLRenderWindow* renWin = ui.qvtkwidget->GetRenderWindow();
    renWin->AddRenderer(m_Renderer);
    m_Interactor = ui.qvtkwidget->GetInteractor();
    m_Interactor->Initialize();
}

MainWindow::~MainWindow() {
    
}

void MainWindow::RemovePolyData(std::string name) {
	vtkProp* prop = m_PropMap[name];
	m_Renderer->RemoveActor(prop);
	m_PropMap.erase(name);
}

void MainWindow::RemoveAllPolyData() {
	m_Renderer->RemoveAllViewProps();
    m_Renderer->ResetCamera();
    m_Interactor->Render();
    m_PropMap.clear();
}

void MainWindow::AddPolyData(std::string name, vtkPolyData* poly, float r, float g, float b, float opacity) {
    vtkPolyDataMapperPointer mapper = vtkPolyDataMapperPointer::New();
    mapper->SetInputData(poly);
    vtkPropertyPointer props = vtkPropertyPointer::New();
    props->SetColor(r, g, b);
    props->SetOpacity(opacity);

    vtkActorPointer actor = vtkActorPointer::New();
    actor->SetMapper(mapper);
    actor->SetProperty(props);


    if (m_PropMap.find(name) != m_PropMap.end()) {
    	m_Renderer->RemoveActor(m_PropMap[name]);
    	m_PropMap.erase(name);
    }

    vtkProp* newActor = actor.GetPointer();
    m_PropMap.insert(NamedProp(name, newActor));

    m_Renderer->AddActor(actor);
    m_Renderer->ResetCamera();
    m_Interactor->Render();

}

void MainWindow::on_resetButton_clicked() {
	RemoveAllPolyData();

    vtkSphereSourcePointer sphereSource = vtkSphereSourcePointer::New();
    sphereSource->SetThetaResolution(64);
    sphereSource->SetPhiResolution(64);
    sphereSource->SetRadius(1);
    sphereSource->Update();
    vtkPolyDataPointer poly(sphereSource->GetOutput());
    vtkPolyDataMapperPointer mapper = vtkPolyDataMapperPointer::New();
    mapper->SetInputData(poly);
    vtkActorPointer actor = vtkActorPointer::New();
    actor->SetMapper(mapper);

    AddPolyData("zAxis", CreateCircle(0, 0, 1, 0, 0, 0, 1.01), 0, 0, 1);
    AddPolyData("yAxis", CreateCircle(0, 1, 0, 0, 0, 0, 1.01), 0, 1, 0);
    AddPolyData("xAxis", CreateCircle(1, 0, 0, 0, 0, 0, 1.01), 1, 1, 0);
    PNSBase::VectorType xPole("1 0 0"), yPole("0 1 0"), zPole("0 0 1");
    AddPolyData("zPole", CreateSinglePoint(zPole), 0, 0, 1);
    AddPolyData("yPole", CreateSinglePoint(yPole), 0, 1, 0);
    AddPolyData("xPole", CreateSinglePoint(xPole), 1, 1, 0);

    m_Renderer->AddActor(actor);
    m_Renderer->ResetCamera();
    m_Interactor->Render();
}

void MainWindow::on_addButton_clicked() {
    PNSBase::VectorType normal(3);
    normal.randu();
    normal /= arma::norm(normal, 2);
    
    PNSBase::CreateSphereRandoms(100, M_PI_4, M_PI / 18.f, normal, m_Math.m_Data);
    vtkPolyDataPointer pointsSpheres = CreatePoints(m_Math.m_Data, 0.01);
    AddPolyData("Samples", pointsSpheres, 1, 0, 0);
}

void MainWindow::on_runButton_clicked() {
    PNSBase::VectorType northPole("0 0 1");
    m_Math.startOptimization(northPole, M_PI_2, ui.tauSpinBox->value());
    vtkPolyDataPointer newCircle = CreateCircle(m_Math.m_Normal, m_Math.m_Phi);
    cout << "R: " << m_Math.m_Phi << endl;

    AddPolyData("principalArc", newCircle, 0, 1, 1);
    AddPolyData("centerPoint", CreateSinglePoint(m_Math.m_Normal, 0.05), 0, 1, 1);
    AddPolyData("centerPointAtTangent", CreateSinglePoint(m_Math.m_CenterAtTangent, 0.05), 1, 0, 0);
}

void MainWindow::on_testButton_clicked() {
    ExpLogTransform expLog;
    PNSBase::VectorType p("0 0 1");
    expLog.SetTangentPoint(p);

    PNSBase::MatrixType tangentSpace, sphereSpace;

    expLog.TransformLog(m_Math.m_Data, tangentSpace);
    expLog.TransformExp(tangentSpace, sphereSpace);

//    expLog.PrintMatrix(cout, m_Math.m_Data);
    vtkPolyDataPointer transformedSamples = CreatePoints(tangentSpace, 0.02);
    vtkPolyDataPointer sphereSamples = CreatePoints(sphereSpace, 0.05);
    AddPolyData("logTransformedSamples", transformedSamples, 1, .7, .7);
    AddPolyData("expTransformedSamples", sphereSamples, .3, .3, 1, .5);
}
