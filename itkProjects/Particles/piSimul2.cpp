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
#include "QGraphicsRectWidget.h"
#include "itkImageIO.h"
#include "piTimer.h"
#include "qutils.h"



using namespace std;

namespace pi {
    ParticleSystemSolver main;
    ParticleSystem& system = main.m_System;
    ImageContext& images = main.m_ImageContext;

    ParticleTrace trace;

    string configcache = "/tmp/psim.txt";
    
    Simul2::Simul2(QWidget* parent): core(this) {
        ui.setupUi(this);

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

        centerToDesktop();
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


    void Simul2::on_applyButton_clicked(bool value) {
        stringstream configstream;
        string config = ui.config->toPlainText().toStdString();
        configstream.str(config);

        cout << configstream.str() << endl;

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
        main.t0 = main.m_Options.GetRealVectorValue("TimeRange:", 0);
        main.t1 = main.m_Options.GetRealVectorValue("TimeRange:", 2);
        main.dt = main.m_Options.GetRealVectorValue("TimeRange:", 1);
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

    void Simul2::on_zoom_sliderMoved(int val) {
        QTransform transform;
        transform.scale(val / 100.0f,  val / 100.0f);
        ui.graphicsView->setTransform(transform);
        ui.graphicsView2->setTransform(transform);
    }


    void Simul2::on_loadTrace_clicked() {
        string traceFile = main.m_Options.GetString("RunTrace:");
        if (traceFile == "") {
            QString file = __fileManager.openFile(QFileManager::Image, this, "Open Particle Trace File");
            if (file.isEmpty()) {
                return;
            }
            traceFile = file.toStdString();
        }

        trace.Clear();
        ifstream fin(traceFile.c_str());
        trace.Read(fin);

        ui.traceSteps->setValue(0);
        ui.traceSteps->setMaximum(trace.system[0].timeSeries.size());
    }

    void Simul2::on_saveTrace_clicked() {
        string traceFile = main.m_Options.GetString("RunTrace:");
        if (traceFile == "") {
            return;
        }
        ofstream traceOut(traceFile.c_str());
        main.trace.Write(traceOut);
        traceOut.close();
    }

    void Simul2::on_traceSteps_valueChanged(int n) {
        if (n >= trace.system[0].timeSeries.size() || n < 0) {
            return;
        }
        ui.lcdNumber->display(n);
        core.updateParticles(0, trace.system[0].timeSeries[n]);
        core.updateParticles(1, trace.system[1].timeSeries[n]);
    }


    void Simul2::on_actionViewOrientation_triggered() {
        if (ui.viewSplitter->orientation() == Qt::Vertical) {
            ui.viewSplitter->setOrientation(Qt::Horizontal);
        } else {
            ui.viewSplitter->setOrientation(Qt::Vertical);
        }
    }

    void Simul2::on_actionShowWarped_triggered() {
//        main.m_System.WarpImage(0, 1);
        RealImage::Pointer warpedImage1 = main.WarpImage(0);
        core.showAuxImage(0, warpedImage1);
        RealImage::Pointer warpedImage2 = main.WarpImage(1);
        core.showAuxImage(1, warpedImage2);
    }

    void Simul2::on_actionPrint_triggered() {
        main.PrintPoints();
    }
}