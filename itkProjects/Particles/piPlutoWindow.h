//
//  piPlutoWindow.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/2/13.
//
//

#ifndef __ParticleGuidedRegistration__piPlutoWindow__
#define __ParticleGuidedRegistration__piPlutoWindow__

#include <iostream>

#include <QMainWindow>
#include "ui_plutowindow.h"

namespace pi {
    class PlutoWindow: public QMainWindow {
        Q_OBJECT
    public:
        PlutoWindow(QWidget* parent = NULL);
        virtual ~PlutoWindow();
        
        void centerToDesktop();
        void connectSignals();
        
    public slots:
        void on_actionOpen_triggered();
        
    private:
        Ui_PlutoMain _ui;
    };
}

#endif /* defined(__ParticleGuidedRegistration__piPlutoWindow__) */
