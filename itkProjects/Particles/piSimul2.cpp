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
#include "piImageSlice.h"
#include "QGraphicsItem"
#include "QGraphicsPixmapItem"

using namespace std;

namespace pi {
    ParticleSystemSolver main;
    ParticleSystem& system = main.m_System;
    ImageContext& images = main.m_ImageContext;
    LabelVector& labels = images.m_LabelImages;


    // storage for particle drawing
    QTransform viewTransform;
    ImageSlice<LabelImage> viewSlice;
    QGraphicsPixmapItem* pixmapItem;
    vector<QGraphicsEllipseItem*> ellipseItems;
    
    string configcache = "psim.txt";
    
    Simul2::Simul2(QWidget* parent) {
        ui.setupUi(this);
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

    void Simul2::on_actionAddLabel_triggered() {

    }

    void Simul2::on_actionNewParticles_triggered() {

    }

    void Simul2::on_applyButton_clicked(bool value) {

        stringstream configstream;
        string config = ui.config->toPlainText().toStdString();
        configstream.str(config);
        main.LoadConfig(configstream);

        // together
        setupParticles();
        updateParticles();

        // constraint view
        viewSlice.alpha = 60;
        viewSlice.SetLabel(labels[0]);
        m_scene.removeItem((QGraphicsItem*) viewSlice.pixmapCache);
        viewSlice.pixmapCache = m_scene.addPixmap(viewSlice.GetPixmap());
        viewSlice.pixmapCache->setZValue(-1);

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
//                p.label = p.x[1] > 80;
                p.label = j%2;
            }
        }
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
}