//
//  piVTK.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/16/13.
//
//

#include "piVTK.h"
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
#include "vtkPolyDataWriter.h"
#include "vtkRenderer.h"
#include "vtkDiskSource.h"
#include "vtkPlaneSource.h"
#include "vtkTriangleFilter.h"
#include "vtkSphereSource.h"

namespace pivtk {
    PolyDataPointer CreateSphere(int phiRes, int thetaRes) {
        __vtk(SphereSource);
        SphereSourcePointer source = SphereSourcePointer::New();
        source->SetPhiResolution(phiRes);
        source->SetThetaResolution(thetaRes);
        source->LatLongTessellationOn();
        source->Update();
        return PolyDataPointer(source->GetOutput());
    }
    
    PolyDataPointer ConstructPathLines(pi::ParticleSetSeries& snapshot) {
        PointsPointer points = PointsPointer::New();
        CellArrayPointer cells = CellArrayPointer::New();
        int maxTime = snapshot.timeSeries.size();
        const int nPoints = snapshot.maxIdx + 1;
        points->SetNumberOfPoints(nPoints * maxTime);
        int idx = 0;
        for (int j = 0; j < nPoints; j++) {
            PolyLinePointer polyLine = PolyLinePointer::New();
            polyLine->GetPointIds()->SetNumberOfIds(maxTime);
            for (int t = 0; t < maxTime; t++) {
                polyLine->GetPointIds()->SetId(t, idx);
                points->SetPoint(idx, snapshot.timeSeries[t][j].x);
                ++idx;
            }
            cells->InsertNextCell(polyLine);
        }

        PolyDataPointer polyData = PolyDataPointer::New();
        polyData->SetPoints(points);
        polyData->SetLines(cells);

        return polyData;
    }

    PropScene::PropScene() {
        m_LastActor = NULL;
        m_Renderer = NULL;
    }

    PropScene::~PropScene() {

    }

    vtkPolyData* PropScene::LoadPolyData(const char *filename) {
        vtkPolyDataReader* reader = vtkPolyDataReader::New();
        reader->SetFileName(filename);
        reader->Update();
        vtkPolyData* poly = reader->GetOutput();
        return poly;
    }

    void PropScene::Render() {
        m_Renderer->Render();
    }

    vtkPolyData* PropScene::CreateDisk(int innerRadius, int outerRadius) {
        vtkDiskSource* diskSource = vtkDiskSource::New();
        diskSource->SetOuterRadius(outerRadius);
        diskSource->SetInnerRadius(innerRadius);
        diskSource->SetCircumferentialResolution(36);
        diskSource->SetRadialResolution(18);
        diskSource->Update();
        return diskSource->GetOutput();
    }

    vtkPolyData* PropScene::CreatePlane(int xr, int yr, double cx, double cy) {
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

    vtkActor* PropScene::AddPolyData(std::string name, vtkPolyData *poly) {
        vtkActor* actor = NewActor(poly);
        AddProp(name, actor);
        return actor;
    }

    vtkActor* PropScene::NewActor(vtkPolyData* poly) {
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

    void PropScene::AddProp(std::string name, vtkProp* prop) {
        if (m_PropMap.find(name) != m_PropMap.end()) {
            vtkProp* foundProp = m_PropMap[name];
            m_Renderer->RemoveActor(foundProp);
            m_PropMap.erase(name);
        }
        m_PropMap.insert(NamedProp(name, prop));
        m_Renderer->AddViewProp(prop);
    }

    vtkProp* PropScene::FindProp(std::string name) {
        PropMapType::iterator it = m_PropMap.find(name);
        if (it != m_PropMap.end()) {
            return it->second;
        }
        return NULL;
    }

    vtkProp* PropScene::RemoveProp(std::string name) {
        vtkProp* prop = FindProp(name);
        if (prop != NULL) {
            m_PropMap.erase(name);
            m_Renderer->RemoveViewProp(prop);
            prop->Delete();
        }
        return prop;
    }

    vtkPolyData* PropScene::FindPolyData(std::string name) {
        vtkActor* actor = FindActor(name);
        if (actor == NULL) {
            return NULL;
        }
        return dynamic_cast<vtkPolyData*>(actor->GetMapper()->GetInput());
    }
    vtkActor* PropScene::FindActor(std::string name) {
        vtkProp* prop = FindProp(name);
        vtkActor* actor = dynamic_cast<vtkActor*>(prop);
        if (actor != NULL) {
            m_LastActor = actor;
            return actor;
        }
        return NULL;
    }

    void PropScene::SetColor(float r, float g, float b) {
        if (m_LastActor != NULL) {
            m_LastActor->GetProperty()->SetColor(r, g, b);
        }
    }

    void PropScene::SetOpacity(float o) {
        if (m_LastActor != NULL) {
            m_LastActor->GetProperty()->SetOpacity(o);
        }
    }
    
    void PropScene::SetRepresentation(int n) {
        if (m_LastActor != NULL) {
            m_LastActor->GetProperty()->SetRepresentation(n);
        }
    }
    
    void PropScene::ModifyLastActor() {
        if (m_LastActor != NULL) {
            m_LastActor->GetMapper()->Modified();
        }
    }
    
    void PropScene::Clear() {
        m_PropMap.clear();
        m_Renderer->RemoveAllViewProps();
    }


    void vtk_write_polydata(const char* f, vtkPolyData* p) {
        vtkPolyDataWriter* w = vtkPolyDataWriter::New();
        w->SetInput(p);
        w->SetFileName(f);
        w->Write();
        w->Delete();
    }
}