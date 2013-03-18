//
//  piVTK.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/16/13.
//
//

#ifndef __ParticleGuidedRegistration__piVTK__
#define __ParticleGuidedRegistration__piVTK__

#include <map>
#include <iostream>
//#include "piParticleTrace.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkIdList.h"
#include "vtkCellArray.h"
#include "vtkPolyLine.h"
#include "vtkRenderer.h"
#include "vtkProp.h"
#include "vtkActor.h"

#define __vtk(x) typedef vtkSmartPointer<vtk##x> x##Pointer
namespace pivtk {
    __vtk(CellArray);
    __vtk(Points);
    __vtk(PolyData);
    __vtk(PolyLine);


    PolyDataPointer CreateSphere(int phiRes, int thetaRes);
    void vtk_write_polydata(const char* f, vtkPolyData* p);
//    PolyDataPointer ConstructPathLines(pi::ParticleSetSeries& snapshot);

    class PropScene {
    public:
        typedef std::map<std::string, vtkProp*> PropMapType;
        typedef std::pair<std::string, vtkProp*> NamedProp;

        PropScene();
        ~PropScene();

        void SetRenderer(vtkRenderer* renderer) { m_Renderer = renderer; }
        vtkRenderer* GetRenderer() { return m_Renderer; }
        void Render();

        vtkPolyData* LoadPolyData(const char* filename);
        vtkActor* AddPolyData(std::string name, vtkPolyData* poly);
        vtkPolyData* FindPolyData(std::string name);
        vtkActor* NewActor(vtkPolyData* poly);


        // scene related function
        void AddProp(std::string name, vtkProp* prop);
        vtkActor* FindActor(std::string name);
        vtkProp* FindProp(std::string name);
        vtkProp* RemoveProp(std::string name);
        void ModifyLastActor();

        vtkPolyData* CreateDisk(int innerRadius, int outerRadius);
        vtkPolyData* CreatePlane(int xr, int yr, double cx, double cy);


        void SetColor(float r, float g, float b);
        void SetOpacity(float o);
        void SetRepresentation(int n);
        
        
        void Clear();
        
    private:
        PropMapType m_PropMap;
        vtkRenderer* m_Renderer;
        vtkActor* m_LastActor;
    };
}

#endif /* defined(__ParticleGuidedRegistration__piVTK__) */
