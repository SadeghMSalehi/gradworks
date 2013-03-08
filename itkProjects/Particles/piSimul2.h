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

namespace pi {
    class Simul2: public QMainWindow {
        Q_OBJECT
    public:
        Simul2(QWidget* parent = NULL);
        void setupParticles();
        void updateParticles();
        void removeParticles();

    public slots:
        void on_applyButton_clicked(bool value);
        void on_runStepButton_clicked();
        void on_showImage_toggled(bool value);


        //////////////////////////////////////////////////
        // VISUALIZATIONS
        void on_loadMoving_clicked();
        void on_loadFixed_clicked();
        void on_loadTrace_clicked();
        void on_showWarped01_clicked();
        void on_showWarped10_clicked();
        void on_printAsPDF_clicked();
        void on_showMoving_toggled(bool checked);
        void on_showFixed_toggled(bool checked);
        void on_traceSteps_valueChanged(int n);
        
        void tick();

    private:
        QColor getColor(int i);
        
        Ui_Simul2D ui;
        QGraphicsScene m_scene;
        QTimer m_timer;
    };

}
#endif /* defined(__ParticleGuidedRegistration__piSimul2__) */
