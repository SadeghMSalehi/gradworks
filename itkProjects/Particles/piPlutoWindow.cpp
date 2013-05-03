//
//  piPlutoWindow.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/2/13.
//
//

#include "piPlutoWindow.h"
#include <QDesktopWidget>

namespace pi {
    PlutoWindow::PlutoWindow(QWidget* parent): QMainWindow(parent) {
        _ui.setupUi(this);
    }
    
    PlutoWindow::~PlutoWindow() {
        
    }
    
    void PlutoWindow::centerToDesktop() {
        QDesktopWidget *desktop = QApplication::desktop();
        
        int screenWidth, width;
        int screenHeight, height;
        int x, y;
        QSize windowSize;
        
        screenWidth = desktop->width();
        screenHeight = desktop->height();
        
        windowSize = size();
        width = windowSize.width();
        height = windowSize.height();
        
        x = (screenWidth - width) / 2;
        y = (screenHeight - height) / 2;
        y -= 50;
        
        move(x, y);
        resize(windowSize.width(), windowSize.height());
    }
    
    void PlutoWindow::connectSignals() {
        
    }
    
    void PlutoWindow::on_actionOpen_triggered() {
        
    }
}