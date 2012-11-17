//
//  vtkPropScene.h
//  laplacePDE
//
//  Created by Joohwi Lee on 11/15/12.
//
//

#ifndef __laplacePDE__vtkPropScene__
#define __laplacePDE__vtkPropScene__

#include <iostream>
#include <string>
#include <map>

class vtkPolyData;
class vtkProp;
class vtkActor;
class vtkRenderer;

class vtkPropScene {
public:
    typedef std::map<std::string, vtkProp*> PropMapType;
    typedef std::pair<std::string, vtkProp*> NamedProp;

    vtkPropScene();
    virtual ~vtkPropScene();
    
    void SetRenderer(vtkRenderer* renderer) { m_Renderer = renderer; }
    vtkRenderer* GetRenderer() { return m_Renderer; }
    void Render();

    vtkPolyData* LoadPolyData(const char* filename);
    vtkActor* AddPolyData(std::string name, vtkPolyData* poly);
    vtkPolyData* FindPolyData(std::string name);
    
    vtkActor* NewActor(vtkPolyData* poly);
    void AddProp(std::string name, vtkProp* prop);
    vtkActor* FindActor(std::string name);
    vtkProp* FindProp(std::string name);
    vtkProp* RemoveProp(std::string name);
    void ModifyLastActor();

    void SetColor(float r, float g, float b);
    void SetOpacity(float o);
    void SetRepresentation(int n);

    
    void Clear();

private:
    PropMapType m_PropMap;
    vtkRenderer* m_Renderer;
    vtkActor* m_LastActor;
};

#endif /* defined(__laplacePDE__vtkPropScene__) */
