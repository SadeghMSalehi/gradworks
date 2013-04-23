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
#include <QMovie>
#include <QtConcurrentRun>

#include "QParticlesGraphicsItem.h"
#include "QGraphicsRectWidget.h"
#include "piTimer.h"
#include "qutils.h"

#include "piImagePatch.h"
#include "piImageEntropyComputer.h"
#include "piParticleBSpline.h"
#include "piImageRegistration.h"
#include "itkSIFTImageFilter.h"
#include "piImageProcessing.h"

using namespace std;

namespace pi {
    ParticleSystemSolver main;
    ParticleSystem& system = main.m_System;
    ImageContext& images = main.m_ImageContext;
    ImageIO<RealImage> __imageIO;

    ParticleTrace trace;

    string configcache = "/tmp/psim.txt";
    
    Simul2::Simul2(QWidget* parent): core(this) {
        ui.setupUi(this);

        static QMovie* movie = new QMovie(":/Icons/Images/loading.gif");
        ui.loadingAnimation->setMovie(movie);

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

    void Simul2::threadedRun() {
        while (main.t < main.t1 && main.continueToRun) {
            main.RunStep();
        }
    }

    void Simul2::on_runStepButton_clicked() {
        if (!ui.runStepButton->isChecked()) {
            stopLoadingAnimation();
        } else {
            startLoadingAnimation();
        }
    }


    void Simul2::startLoadingAnimation() {
        if (m_future.isRunning()) {
            return;
        }
        stringstream configstream;
        string config = ui.config->toPlainText().toStdString();
        configstream.str(config);

        cout << configstream.str() << endl;

        main.LoadParameters(configstream);
        main.SetupParameters();
        main.continueToRun = true;
        main.verbose = ui.verboseOutput->isChecked();


        // prepare simultation
        main.t0 = main.m_Options.GetRealVectorValue("TimeRange:", 0);
        main.t1 = main.m_Options.GetRealVectorValue("TimeRange:", 2);
        main.dt = main.m_Options.GetRealVectorValue("TimeRange:", 1);
        main.t = 0;

        // UI status
        ui.runStepButton->setText("Stop");
        ui.runStepButton->setChecked(true);
        ui.loadingAnimation->movie()->start();

        m_timer.start(500);
        m_future = QtConcurrent::run(this, &Simul2::threadedRun);
        m_futureWatcher.setFuture(m_future);

        // Connect signals
        connect(&m_futureWatcher, SIGNAL(finished()), SLOT(stopLoadingAnimation()));
        connect(&m_futureWatcher, SIGNAL(canceled()), SLOT(stopLoadingAnimation()));
    }

    void Simul2::stopLoadingAnimation() {
        ui.loadingAnimation->movie()->stop();
        ui.runStepButton->setText("Step");
        ui.runStepButton->setChecked(false);
        m_timer.stop();
        if (m_future.isRunning()) {
            main.continueToRun = false;
        }
        disconnect(&m_futureWatcher, SIGNAL(finished()), this, SLOT(stopLoadingAnimation()));
        disconnect(&m_futureWatcher, SIGNAL(canceled()), this, SLOT(stopLoadingAnimation()));
        m_future = QFuture<void>();
    }

    void Simul2::tick() {
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

    void Simul2::on_actionBsplineWarp_triggered() {
        ParticleBSpline bspline;
        ImageIO<LabelImage> io;
        LabelImage::Pointer labelImg = io.CastImageFromS<RealImage>(core.getImage(0));
        bspline.SetReferenceImage(labelImg);
        bspline.EstimateTransform(main.m_System[0], main.m_System[1]);
        RealImage::Pointer warpedImage = bspline.WarpImage(core.getImage(1));
        core.showAuxImage(0, warpedImage);
    }

    void Simul2::on_actionImageBsplineWarp_triggered() {
        RealImage::Pointer warpedImage = bsplineRegistration(core.getImage(1), core.getImage(0));
        core.showAuxImage(1, warpedImage);
    }

    void Simul2::on_actionIntensityGradient_triggered() {
        core.showAuxImage(0, main.m_System[0].kappaImage, true);
    }

    void Simul2::on_actionShowWarped_triggered() {
//        main.m_System.WarpImage(0, 1);
        RealImage::Pointer warpedImage1 = main.WarpImage(0);
        core.showAuxImage(0, warpedImage1);
        RealImage::Pointer warpedImage2 = main.WarpImage(1);
        core.showAuxImage(1, warpedImage2);
    }

    void Simul2::on_actionShowCostPlot_toggled(bool toggle) {
        if (toggle) {
            ui.costPlot->show();
        } else {
            ui.costPlot->hide();
        }
    }

    void Simul2::on_actionPrint_triggered() {
        main.PrintPoints();
    }

    void Simul2::on_actionTest_triggered() {
        core.showEntropyMeasurementImage();
    }

    void Simul2::on_actionThreadTest_toggled(bool toggle) {
        startLoadingAnimation();
    }

    void Simul2::on_actionGenerateSIFTImage1_triggered() {
        ImageProcessing proc;
        GradientImage::Pointer gradImg = proc.ComputeGaussianGradient(core.getImage(0), 1);
        itk::SIFTImageFilter::Pointer filter = itk::SIFTImageFilter::New();
        filter->SetInput(gradImg);
        try {
            filter->Update();
        } catch (itk::ExceptionObject& ex) {
            cout << ex << endl;
        }
        itk::SIFTImage::Pointer siftImage = filter->GetOutput();
        ImageIO<itk::SIFTImage> io;
        io.WriteImage("/tmpfs/sift1.nrrd", siftImage);
        
        itk::SIFTImagePCAComputer pcaComp;
        pcaComp.computePCA(siftImage);
        RealImage::Pointer pc1 = pcaComp.computePCImage(siftImage, 0);
        RealImage::Pointer pc2 = pcaComp.computePCImage(siftImage, 1);
        RealImage::Pointer pc3 = pcaComp.computePCImage(siftImage, 2);

        core.showCompositeImage(pc1, pc2, pc3);
    }

    void Simul2::on_actionGenerateSIFTImage2_triggered() {
        ImageProcessing proc;
        GradientImage::Pointer gradImg = proc.ComputeGaussianGradient(core.getImage(1), 1);
        itk::SIFTImageFilter::Pointer filter = itk::SIFTImageFilter::New();
        filter->SetInput(gradImg);
        try {
            filter->Update();
        } catch (itk::ExceptionObject& ex) {
            cout << ex << endl;
        }
        itk::SIFTImage::Pointer siftImage = filter->GetOutput();
        ImageIO<itk::SIFTImage> io;
        io.WriteImage("/tmpfs/sift2.nii.gz", siftImage);
    }

}