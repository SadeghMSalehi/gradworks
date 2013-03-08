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
#include "QGraphicsItem"
#include "QGraphicsPixmapItem"
#include "itkImageIO.h"



using namespace std;

namespace pi {
    itkcmds::itkImageIO<LabelImage> imageIO;

    ParticleSystemSolver main;
    ParticleSystem& system = main.m_System;
    ImageContext& images = main.m_ImageContext;
    LabelVector& labels = images.m_LabelImages;
    ParticleTrace trace;

    // storage for particle drawing
    QTransform viewTransform;
    ImageSlice<LabelImage> viewSlice;
    QGraphicsPixmapItem* pixmapItem;
    vector<QGraphicsEllipseItem*> ellipseItems;

    ImageSlice<LabelImage> sourceSlice;
    ImageSlice<LabelImage> targetSlice;
    
    string configcache = "psim.txt";
    
    Simul2::Simul2(QWidget* parent) {
        ui.setupUi(this);
//        ui.toolBar->addWidget(ui.showImage);

        ui.graphicsView->setScene(&m_scene);
        ui.graphicsView->setBackgroundBrush(QBrush(Qt::black));

        viewTransform.scale(4,4);
        ui.graphicsView->setTransform(viewTransform);

        // signal-slot connection
        QObject::connect(&m_timer, SIGNAL(timeout()), this, SLOT(tick()));

        // load config cache
        ifstream in(configcache.c_str());
        if (in) {
            std::stringstream buffer;
            buffer << in.rdbuf();
            ui.config->setPlainText(buffer.str().c_str());
        }
    }

    void Simul2::on_showImage_toggled(bool check) {
        if (check) {
            // constraint view
            viewSlice.alpha = 60;
            viewSlice.SetImage(labels[0]);
            viewSlice.pixmapCache = m_scene.addPixmap(viewSlice.GetPixmap());
            viewSlice.pixmapCache->setZValue(-1);
        } else {
            m_scene.removeItem((QGraphicsItem*) viewSlice.pixmapCache);
        }
    }

    void Simul2::on_applyButton_clicked(bool value) {

        stringstream configstream;
        string config = ui.config->toPlainText().toStdString();
        configstream.str(config);
        main.LoadConfig(configstream);

        // together
        setupParticles();
        updateParticles();

        ui.showImage->setChecked(true);

        // save current configuration
        ofstream cfg(configcache.c_str());
        cfg << config;
        cfg.close();

        ui.statusbar->showMessage(QString("%1 particles loaded!").arg(system.GetNumberOfParticles()));
    }

    void Simul2::on_runStepButton_clicked() {
        // prepare simultation
        main.t0 = 0;
        main.t1 = 10;
        main.dt = 0.05;
        main.t = 0;

        main.verbose = false;
        main.Setup();
        m_timer.start(200);
    }
    
    void Simul2::tick() {
        if (main.t > main.t1) {
            main.t = 0;
            m_timer.stop();
            return;
        }
        // procede 1 secs
        DataReal nextT = main.t + 2;
        while (main.t < nextT) {
            main.RunStep();
        }
        updateParticles();
    }

    void Simul2::setupParticles() {
        if (system.size() < 1 || labels.size() < 1) {
            return;
        }
        if (labels[0].IsNull()) {
            return;
        }
        system[0].InitializeRandomPoints(labels[0]);
        for (int i = 0; i < system.size(); i++) {
            for (int j = 0; j < system[i].size(); j++) {
                Particle& p = system[i][j];
                p.label = 1;
//                p.label = p.x[1] > 80;
//                p.label = j%2;
            }
        }
        system[0].friendImage = viewSlice.GetImage();
    }

    void Simul2::updateParticles() {
        removeParticles();
        for (int i = 0; i < system.size(); i++) {
            for (int j = 0; j < system[i].size(); j++) {
                Particle& p = system[i][j];
                QColor pointColor = getColor(p.label);
                QGraphicsEllipseItem* item = m_scene.addEllipse(p.x[0]-.5, p.x[1]-.5, 1, 1, QPen(pointColor), QBrush(pointColor, Qt::SolidPattern));
                ellipseItems.push_back(item);
            }
        }
    }

    void Simul2::removeParticles() {
        for (int i = 0; i < ellipseItems.size(); i++) {
            m_scene.removeItem((QGraphicsItem*)ellipseItems[i]);
        }
        ellipseItems.clear();
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


    void Simul2::on_loadMoving_clicked() {
        LabelImage::Pointer moving = imageIO.ReadImageT("/NIRAL/work/joohwi/data/synth/fixed2d.nii.gz");
        sourceSlice.alpha = 60;
        sourceSlice.SetImage(moving);
    }

    void Simul2::on_loadFixed_clicked() {
        LabelImage::Pointer fixed = imageIO.ReadImageT("/NIRAL/work/joohwi/data/synth/moving2d.nii.gz");
        targetSlice.alpha = 60;
        targetSlice.SetImage(fixed);
    }

    void Simul2::on_loadTrace_clicked() {
        ifstream is("/NIRAL/work/joohwi/data/synth/run_trace_f2dm2d.txt");
        trace.Read(is);
        ui.traceSteps->setMaximum(trace.system[0].timeSeries.size()-1);
        ui.lcdNumber->display(0);
    }

    void Simul2::on_traceSteps_valueChanged(int n) {
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
    }

    void Simul2::on_showWarped01_clicked() {

    }

    void Simul2::on_showWarped10_clicked() {

    }

    void Simul2::on_printAsPDF_clicked() {

    }

    void Simul2::on_showMoving_toggled(bool checked) {
        if (checked) {
            // constraint view
            sourceSlice.alpha = 60;
            sourceSlice.pixmapCache = m_scene.addPixmap(sourceSlice.GetPixmap());
            sourceSlice.pixmapCache->setZValue(-1);
        } else {
            m_scene.removeItem((QGraphicsItem*) sourceSlice.pixmapCache);
        }
    }

    void Simul2::on_showFixed_toggled(bool checked) {
        if (checked) {
            // constraint view
            targetSlice.alpha = 60;
            targetSlice.pixmapCache = m_scene.addPixmap(targetSlice.GetPixmap());
            targetSlice.pixmapCache->setZValue(-1);
        } else {
            m_scene.removeItem((QGraphicsItem*) targetSlice.pixmapCache);
        }
    }
}