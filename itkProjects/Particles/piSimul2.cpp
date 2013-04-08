//
//  piSimul2.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/23/13.
//
//

#include "piSimul2.h"
#include "piParticleCore.h"
#include "piParticleSystemSolver.h"
#include "piParticleTrace.h"
#include "piImageSlice.h"

#include <QDesktopWidget>
#include <QGraphicsItem>
#include <QGraphicsPixmapItem>

#include "QParticlesGraphicsItem.h"
#include "itkImageIO.h"
#include "piTimer.h"
#include "qutils.h"



using namespace std;

namespace pi {
    itkcmds::itkImageIO<LabelImage> imageIO;

    ParticleSystemSolver main;
    ParticleSystem& system = main.m_System;
    ImageContext& images = main.m_ImageContext;

    string configcache = "/tmp/psim.txt";
    
    Simul2::Simul2(QWidget* parent): core(this) {
        ui.setupUi(this);
        centerToDesktop();

        // signal-slot connection
        QObject::connect(&m_timer, SIGNAL(timeout()), this, SLOT(tick()));

        // load config cache
        ifstream in(configcache.c_str());
        if (in) {
            std::stringstream buffer;
            buffer << in.rdbuf();
            ui.config->setPlainText(buffer.str().c_str());
        }

        core.setUi(&ui);
        ui.costPlot->hide();
    }


    void Simul2::centerToDesktop() {
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


    void Simul2::on_actionTest_triggered() {
        core.openImage1("/NIRAL/work/joohwi/nadia/SliceImages/C31/C31_E04_slice_128.nii.gz");
        core.openImage2("/NIRAL/work/joohwi/nadia/SliceImages/C33/C33_E04_slice_128.nii.gz");
        core.openLabel1("/NIRAL/work/joohwi/nadia/SliceImages/C31/C31_E04_label_128.nii.gz");
        core.openLabel2("/NIRAL/work/joohwi/nadia/SliceImages/C33/C33_E04_label_128.nii.gz");
    }

    void Simul2::on_actionPrint_triggered() {
        cout << main.m_Options << endl;
        main.m_System.GetInitialSubject().WriteParticlePositions(cout);
    }

    void Simul2::on_loadImage1_clicked() {
        core.openImage1(__fileManager.openFile(QFileManager::Image, this, "Open Gray Image1"));
    }
    void Simul2::on_loadImage2_clicked() {
        core.openImage2(__fileManager.openFile(QFileManager::Image, this, "Open Gray Image2"));
    }
    void Simul2::on_loadLabel1_clicked() {
        core.openLabel1(__fileManager.openFile(QFileManager::Image, this, "Open Label1"));
    }
    void Simul2::on_loadLabel2_clicked() {
        core.openLabel2(__fileManager.openFile(QFileManager::Image, this, "Open Label2"));
    }


    void Simul2::on_applyButton_clicked(bool value) {
        ui.splitter->setOrientation(Qt::Vertical);

        stringstream configstream;
        string config = ui.config->toPlainText().toStdString();
        configstream.str(config);

        main.LoadConfig(configstream);
        main.Preprocessing();
        main.Setup();

        core.setParticleSolver(&main);
        core.updateParticles();

        stringstream cachestream;
        main.SaveConfig("/tmp/psim.txt");
        cachestream << main.m_Options << endl;
        ui.config->setPlainText(QString::fromStdString(cachestream.str()));
        
        ui.statusbar->showMessage(QString("%1 particles loaded!").arg(system.GetNumberOfParticles()));
    }

    void Simul2::on_runStepButton_clicked() {
        // prepare simultation
        main.t0 = main.m_Options.GetRealVectorValue("TimeRange", 0);
        main.t1 = main.m_Options.GetRealVectorValue("TimeRange", 2);
        main.dt = main.m_Options.GetRealVectorValue("TimeRange", 1);
        main.t = 0;

        main.verbose = ui.verboseOutput->isChecked();
        m_timer.start(100);
    }
    
    void Simul2::tick() {
        if (main.t > main.t1) {
            main.t = 0;
            m_timer.stop();
            return;
        }
        // procede 1 secs
        DataReal nextT = main.t + 1;
        while (main.t < nextT) {
            main.RunStep();
        }
        core.updateParticles();
    }


    QColor Simul2::getColor(int i) {
        switch (i) {
            case 0:
                return Qt::green;
            case 1:
                return Qt::yellow;
            case 2:
                return Qt::cyan;
            case 3:
                return Qt::blue;
            default:
                return Qt::white;
        }
    }

    void Simul2::on_traceSteps_valueChanged(int n) {
        /*
        if (n >= trace.system[0].timeSeries.size() || n < 0) {
            return;
        }
        ui.lcdNumber->display(n);
        removeParticles();
        for (int i = 0; i < trace.system.size(); i++) {
            for (int j = 0; j < trace.system[i].timeSeries[n].size(); j++) {
                Particle& p = trace.system[i].timeSeries[n][j];
                QColor pointColor = getColor(i);
                pointColor.setAlpha(60);
                QGraphicsEllipseItem* item = m_scene.addEllipse(p.x[0]-.5, p.x[1]-.5, 1, 1, QPen(pointColor), QBrush(pointColor, Qt::SolidPattern));
                ellipseItems.push_back(item);
            }
        }
         */
    }

    void Simul2::on_showWarped01_clicked() {

    }

    void Simul2::on_showWarped10_clicked() {

    }

    void Simul2::on_printAsPDF_clicked() {

    }
}