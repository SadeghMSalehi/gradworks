//
//  particleAni.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/14/13.
//
//

#ifndef __ParticlesGUI__particleAni__
#define __ParticlesGUI__particleAni__

#include <iostream>

#include "QMainWindow"
#include "ui_particleAni.h"

class vtkRenderer;
class vtkPoints;
class vtkUnstructuredGrid;

class AniWindow : public QMainWindow {
    Q_OBJECT

public:
    AniWindow(QWidget* parent = NULL);
    ~AniWindow();
    void CreateParticles();

public slots:
    void on_actionOpen_Trace_triggered();
    void on_actionOpen_System_triggered();
    void on_actionForward_triggered();
    void on_actionBackward_triggered();
    void on_actionFirst_triggered();
    void on_actionLast_triggered();
    void on_actionLabel_Registration_triggered();
    void on_actionMark_At_Image_triggered();
    void on_actionCreate_Pathline_triggered();
    void on_timeStepSpin_valueChanged(int value);
    void on_timeSteps_valueChanged(int value);


    // Tools
    void on_actionReset_Camera_triggered();
    void on_actionImage_Particle_View_triggered(bool checked);
    void on_subjects_currentIndexChanged(int n);
    void on_glyphRadius_valueChanged(double r);
    
private:
    void ShowTraceParticles();

private:
    Ui::AniWindow ui;
    vtkRenderer* m_Renderer;
    vtkPoints* m_Points;
    vtkUnstructuredGrid* m_Grid;
};
#endif /* defined(__ParticlesGUI__particleAni__) */
