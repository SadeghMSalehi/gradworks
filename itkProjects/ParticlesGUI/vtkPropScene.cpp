//
//  vtkPropScene.cpp
//  laplacePDE
//
//  Created by Joohwi Lee on 11/15/12.
//
//

#include "vtkPropScene.h"

#include "vtkProp.h"
#include "vtkProperty.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyData.h"
#include "vtkActor.h"
#include "vtkSphereSource.h"
#include "vtkGlyph3D.h"
#include "vtkSmartPointer.h"
#include "vtkProperty.h"
#include "vtkPolyDataReader.h"
#include "vtkRenderer.h"
#include "vtkDiskSource.h"
#include "vtkPlaneSource.h"
#include "vtkTriangleFilter.h"

vtkPropScene::vtkPropScene() {
    m_LastActor = NULL;
    m_Renderer = NULL;
}

vtkPropScene::~vtkPropScene() {

}

vtkPolyData* vtkPropScene::LoadPolyData(const char *filename) {
    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    reader->SetFileName(filename);
    reader->Update();
    vtkPolyData* poly = reader->GetOutput();
    return poly;
}

void vtkPropScene::Render() {
    m_Renderer->Render();
}

vtkPolyData* vtkPropScene::CreateDisk(int innerRadius, int outerRadius) {
    vtkDiskSource* diskSource = vtkDiskSource::New();
    diskSource->SetOuterRadius(outerRadius);
    diskSource->SetInnerRadius(innerRadius);
    diskSource->SetCircumferentialResolution(36);
    diskSource->SetRadialResolution(18);
    diskSource->Update();
    return diskSource->GetOutput();
}

vtkPolyData* vtkPropScene::CreatePlane(int xr, int yr, double cx, double cy) {
    vtkPlaneSource* planeSource = vtkPlaneSource::New();
    planeSource->SetResolution(10, 10);
    planeSource->SetCenter(cx, cy, 0);
    planeSource->SetNormal(0, 0, 1);
    planeSource->Update();
//    vtkTriangleFilter* filter = vtkTriangleFilter::New();
//    filter->SetInputConnection(planeSource->GetOutputPort());
//    filter->PassLinesOff();
//    filter->PassVertsOff();
//    filter->Update();
    cout << "Setting done" << endl;
    return planeSource->GetOutput();
}

vtkActor* vtkPropScene::AddPolyData(std::string name, vtkPolyData *poly) {
    vtkActor* actor = NewActor(poly);
    AddProp(name, actor);
    return actor;
}

vtkActor* vtkPropScene::NewActor(vtkPolyData* poly) {
    vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
    #if VTK_MAJOR_VERSION >= 6
    mapper->SetInputData(poly);
    #else
    mapper->SetInput(poly);
    #endif

    vtkProperty* props = vtkProperty::New();
    vtkActor* actor = vtkActor::New();
    actor->SetMapper(mapper);
    actor->SetProperty(props);
    m_LastActor = actor;
    return actor;
}

void vtkPropScene::AddProp(std::string name, vtkProp* prop) {
    if (m_PropMap.find(name) != m_PropMap.end()) {
    	m_PropMap.erase(name);
    	m_Renderer->RemoveActor(m_PropMap[name]);
    }
    m_PropMap.insert(NamedProp(name, prop));
    m_Renderer->AddViewProp(prop);
}

vtkProp* vtkPropScene::FindProp(std::string name) {
    PropMapType::iterator it = m_PropMap.find(name);
    if (it != m_PropMap.end()) {
        return it->second;
    }
    return NULL;
}

vtkProp* vtkPropScene::RemoveProp(std::string name) {
    vtkProp* prop = FindProp(name);
    if (prop != NULL) {
        m_PropMap.erase(name);
        m_Renderer->RemoveViewProp(prop);
        prop->Delete();
    }
    return prop;
}

vtkPolyData* vtkPropScene::FindPolyData(std::string name) {
    vtkActor* actor = FindActor(name);
    if (actor == NULL) {
        return NULL;
    }
    return dynamic_cast<vtkPolyData*>(actor->GetMapper()->GetInput());
}
vtkActor* vtkPropScene::FindActor(std::string name) {
    vtkProp* prop = FindProp(name);
    vtkActor* actor = dynamic_cast<vtkActor*>(prop);
    if (actor != NULL) {
        m_LastActor = actor;
        return actor;
    }
    return NULL;
}

void vtkPropScene::SetColor(float r, float g, float b) {
    if (m_LastActor != NULL) {
        m_LastActor->GetProperty()->SetColor(r, g, b);
    }
}

void vtkPropScene::SetOpacity(float o) {
    if (m_LastActor != NULL) {
        m_LastActor->GetProperty()->SetOpacity(o);
    }
}

void vtkPropScene::SetRepresentation(int n) {
    if (m_LastActor != NULL) {
        m_LastActor->GetProperty()->SetRepresentation(n);
    }
}

void vtkPropScene::ModifyLastActor() {
    if (m_LastActor != NULL) {
        m_LastActor->GetMapper()->Modified();
    }
}

void vtkPropScene::Clear() {
    m_PropMap.clear();
    m_Renderer->RemoveAllViewProps();
}
