//
//  piSimul2.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/23/13.
//
//

#ifndef __ParticleGuidedRegistration__piSimul2__
#define __ParticleGuidedRegistration__piSimul2__

#include <iostream>
#include "QMainWindow"
#include "QTimer"

#include "ui_simul2d.h"
#include "piqSimulCore.h"

namespace pi {
    class Simul2: public QMainWindow {
        Q_OBJECT
    public:
        Simul2(QWidget* parent = NULL);
        void centerToDesktop();
        
    public slots:
        void on_applyButton_clicked(bool value);
        void on_runStepButton_clicked();

        //////////////////////////////////////////////////
        // VISUALIZATIONS
        void on_zoom_sliderMoved(int);

        void on_loadTrace_clicked();
        void on_saveTrace_clicked();
        void on_traceSteps_valueChanged(int n);

        void on_actionViewOrientation_triggered();
        void on_actionShowWarped_triggered();
        void on_actionPrint_triggered();
        void on_actionTest_triggered();
        
        void tick();

    private:
        QColor getColor(int i);
        
        Ui_Simul2D ui;
        piq::SimulCore core;
        
        QGraphicsScene m_scene;
        QTimer m_timer;

    };

}
#endif /* defined(__ParticleGuidedRegistration__piSimul2__) */
