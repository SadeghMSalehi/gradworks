//
//  piGroupSimul.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/28/13.
//
//

#include <sstream>

#include "piGroupSimul.h"
#include "piParticleSystemSolver.h"
#include "piParticleBSpline.h"
#include "piImageIO.h"
#include "piImageProcessing.h"
#include <QMovie>

using namespace std;

namespace pi {

    pi::ImageIO<RealImage> __realIO;
    pi::ImageIO<LabelImage> __labelIO;

    GroupSimul::GroupSimul() {

    }

    GroupSimul::~GroupSimul() {

    }

    void GroupSimul::setupUi(Ui_Simul2D *ui) {
        this->_ui = ui;
        this->_ui->graphicsView2->hide();
        this->_ui->graphicsView->setScene(&_scene);

        connectSignals();
    }

    void GroupSimul::connectSignals() {
        connect(_ui->actionGrid, SIGNAL(toggled(bool)), SLOT(showGrid(bool)));
        connect(_ui->actionShowParticles, SIGNAL(toggled(bool)), SLOT(showParticles(bool)));
        connect(_ui->actionDistanceMap, SIGNAL(triggered()), SLOT(createDistanceMap()));

    }

    void GroupSimul::setupParticles(pi::ParticleSystemSolver *solver) {
        this->_solver = solver;
    }

    void GroupSimul::loadConfig() {
        stringstream configstream;
        string config = _ui->config->toPlainText().toStdString();
        configstream.str(config);

        cout << configstream.str() << endl;

        _solver->LoadConfig(configstream);
        _solver->Preprocessing();
        _solver->Setup();

        _nSubjs = _solver->m_System.GetNumberOfSubjects();
        _nParticles = _solver->m_System.GetNumberOfParticles();

        stringstream cachestream;
        _solver->SaveConfig("/tmp/psim_group.txt");
        cachestream << _solver->m_Options;
        _ui->config->setPlainText(QString::fromStdString(cachestream.str()));
        _ui->statusbar->showMessage(QString("%1 particles loaded!").arg(_solver->m_System.GetNumberOfParticles()));

        initImages();
        initParticles();
    }

    void GroupSimul::initImages() {
        _images.resize(_nSubjs);
        _warpedImages.resize(_nSubjs);
        _labels.resize(_nSubjs);
        _imageItems.resize(_nSubjs);
        _secondImageItems.resize(_nSubjs);
        _particleGroups.resize(_nSubjs);

        int hPos = 0;
        int wPos = 0;
        for (int i = 0; i < _nSubjs; i++) {
            ParticleSubject& subj = _solver->m_System[i];
            // source images
            _images[i] = subj.GetImage(0);
            _labels[i] = subj.GetLabel();

            // image items
            _imageItems[i] = new QGraphicsRealImageItem();
            _imageItems[i]->setFlip(QGraphicsRealImageItem::UpDown);

            // compute visualization
            if (i == 0) {
                _imageItems[i]->setImage<RealImage>(_images[i], true);
            } else {
                _imageItems[i]->setImage<RealImage>(_images[i], false);
                _imageItems[i]->setRange(_imageItems[0]->getRange(0), _imageItems[0]->getRange(1));
            }
            _imageItems[i]->refresh();
            _imageItems[i]->setPos(0, hPos);


            wPos = _imageItems[0]->boundingRect().width() + 3;
            _secondImageItems[i] = new QGraphicsRealImageItem();
            _secondImageItems[i]->setFlip(QGraphicsRealImageItem::UpDown);
            _secondImageItems[i]->setPos(wPos, hPos);
            _secondImageItems[i]->hide();



            // add item
            _scene.addItem(_imageItems[i]);
            _scene.addItem(_secondImageItems[i]);

            // particle showing
            _particleGroups[i].setListener(this);
            _particleGroups[i].hideParticles(true);
            _particleGroups[i].setScene(&_scene);
            _particleGroups[i].setParentItem(_imageItems[i]);
            _particleGroups[i].createParticles(&(_solver->m_System[i]));

            hPos += (_imageItems[i]->boundingRect().height() + 3);
        }

        _ui->graphicsView->fitInView(0, 0, wPos, hPos, Qt::KeepAspectRatio);
    }

    void GroupSimul::initParticles() {
        
    }

    void GroupSimul::spreadParticles() {
        _solver->SpreadParticles();
        showParticles(true);
    }

    void GroupSimul::showGrid(bool on) {
        for (int i = 0; i < _nSubjs; i++) {
            if (_imageItems[i] != NULL) {
                _imageItems[i]->showGrid(on);
                _secondImageItems[i]->showGrid(on);
            }
        }
    }

    void GroupSimul::showParticles(bool on) {
        for (int i = 0; i < _nSubjs; i++) {
            _particleGroups[i].hideParticles(!on);
            _particleGroups[i].updateParticles();
        }
    }

    void GroupSimul::computeWarpedImages() {
        const int n = _solver->m_System.GetNumberOfSubjects();
        ParticleSubject& meanSubj = _solver->m_System[0];

        for (int i = 0; i < n; i++) {
            ParticleBSpline bspline;
            bspline.SetReferenceImage(_labels[i]);
            bspline.EstimateTransform(meanSubj, _solver->m_System[i]);
            _warpedImages[i] = bspline.WarpImage(_images[i]);
            FieldTransformType::Pointer transform = bspline.GetTransform();

            // show warpedImage
            _secondImageItems[i]->setImage<RealImage>(_warpedImages[i]);
            _secondImageItems[i]->setRange(_imageItems[0]->getRange(0), _imageItems[0]->getRange(1));
            _secondImageItems[i]->generateUserGrids<FieldTransformType>(transform);
            _secondImageItems[i]->refresh();
            _secondImageItems[i]->show();

            __realIO.WriteImage(QString("/tmpfs/WarpedImage_%1.nii.gz").arg(i).toStdString(), _warpedImages[i]);
        }
    }

    void GroupSimul::startRun() {
        if (m_future.isRunning()) {
            return;
        }
        stringstream configstream;
        string config = _ui->config->toPlainText().toStdString();
        configstream.str(config);

        cout << configstream.str() << endl;

        _solver->LoadParameters(configstream);
        _solver->SetupParameters();
        _solver->continueToRun = true;
        _solver->verbose = _ui->verboseOutput->isChecked();


        // prepare simultation
        _solver->t0 = _solver->m_Options.GetRealVectorValue("TimeRange:", 0);
        _solver->t1 = _solver->m_Options.GetRealVectorValue("TimeRange:", 2);
        _solver->dt = _solver->m_Options.GetRealVectorValue("TimeRange:", 1);
        _solver->t = 0;

        // UI status
        _ui->runStepButton->setText("Stop");
        _ui->runStepButton->setChecked(true);
        _ui->loadingAnimation->movie()->start();

        m_timer.start(500, this);
        m_future = QtConcurrent::run(this, &GroupSimul::threadedRun);
        m_futureWatcher.setFuture(m_future);

        // Connect signals
        connect(&m_futureWatcher, SIGNAL(finished()), SLOT(stopLoadingAnimation()));
        connect(&m_futureWatcher, SIGNAL(canceled()), SLOT(stopLoadingAnimation()));
    }


    void GroupSimul::threadedRun() {
        while (_solver->t < _solver->t1 && _solver->continueToRun) {
            _solver->RunStep();
        }
    }

    void GroupSimul::timerEvent(QTimerEvent* event) {
        showParticles(_ui->actionShowParticles->isChecked());
    }

    void GroupSimul::stopRun() {
        _ui->loadingAnimation->movie()->stop();
        _ui->runStepButton->setText("Step");
        _ui->runStepButton->setChecked(false);
        m_timer.stop();
        if (m_future.isRunning()) {
            _solver->continueToRun = false;
        }
        disconnect(&m_futureWatcher, SIGNAL(finished()), this, SLOT(stopRun()));
        disconnect(&m_futureWatcher, SIGNAL(canceled()), this, SLOT(stopRun()));
        m_future = QFuture<void>();
    }

    void GroupSimul::mousePressed(QGraphicsEllipseEventItem *sender, QGraphicsSceneMouseEvent *event) {
        event->accept();
    }

    void GroupSimul::mouseReleased(QGraphicsEllipseEventItem *sender, QGraphicsSceneMouseEvent *event) {
        int particleId = sender->data(0).value<int>();
        for (int i = 0; i < _particleGroups.size(); i++) {
            _particleGroups[i].selectParticle(particleId);
        }
        
        IntensityForce& force = _solver->intensityForce;
        ParticleAttribute& attr = force.GetAttribute(0, particleId);

        cout << attr.F[0] << "," << attr.F[1] << endl;
    }

    void GroupSimul::createDistanceMap() {
        LabelImage::Pointer image = _solver->m_System[0].GetLabel();
        ImageProcessing proc;
        VectorImage::Pointer dmap = proc.ComputeDistanceMap(image);
        RealImage::Pointer dmapMag = proc.ComputeMagnitudeMap(dmap);
        __realIO.WriteImage("/tmpfs/dmap.nii.gz", dmapMag);


    }
}