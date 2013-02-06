#ifndef __particle_viewer_window_h__
#define __particle_viewer_window_h__

#include "ui_particleViewerWindow.h"
#include "QGraphicsScene"

class ParticleViewerWindow : public QMainWindow {
    Q_OBJECT
public:
    ParticleViewerWindow(QWidget* parent = NULL);
    virtual ~ParticleViewerWindow();

    void load(char* file);

public slots:
    void on_action_Open_triggered();
    void on_action_Close_triggered();
    void updateScene();


private:
    void createGrid();
    
    Ui::MainWindow ui;
    QGraphicsScene m_Scene;
};

#endif