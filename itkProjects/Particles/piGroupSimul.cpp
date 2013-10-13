//
//  piGroupSimul.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/28/13.
//
//

#include <sstream>
#include <QMovie>

#include "piGroupSimul.h"
#include "piImageRegistration.h"
#include "piParticleSystemSolver.h"
#include "piParticleBSpline.h"
#include "piParticleTrainer.h"
#include "piImageIO.h"
#include "piImageProcessing.h"
#include "QGraphicsDirectionItem.h"

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
        //this->_ui->graphicsView2->hide();
        this->_ui->graphicsView->setScene(&_scene);
        this->_ui->graphicsView2->setScene(&_scene2);

        connectSignals();

        _particleWindow = new QGraphicsRealImageItem();
    }

    void GroupSimul::connectSignals() {
        connect(_ui->actionGrid, SIGNAL(toggled(bool)), SLOT(showGrid(bool)));
        connect(_ui->actionShowParticles, SIGNAL(toggled(bool)), SLOT(showParticles(bool)));
        connect(_ui->actionDistanceMap, SIGNAL(triggered()), SLOT(createDistanceMap()));
        connect(_ui->actionShowAttributes, SIGNAL(triggered()), SLOT(showAttributes()));
        connect(_ui->actionShowWarped, SIGNAL(triggered()), SLOT(showWarpedImages()));
        connect(_ui->actionImageBsplineWarp, SIGNAL(triggered()), SLOT(computeImageBspline()));
        connect(_ui->actionSaveTransform, SIGNAL(triggered()), SLOT(saveTransform()));
        connect(_ui->actionParticleBsplineWarp2, SIGNAL(triggered()), SLOT(computeParticleWarp()));

        connect(_ui->actionLevel0, SIGNAL(triggered()), SLOT(setResolutionLevel()));
        connect(_ui->actionLevel1, SIGNAL(triggered()), SLOT(setResolutionLevel()));
        connect(_ui->actionLevel2, SIGNAL(triggered()), SLOT(setResolutionLevel()));
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


        // load images
        initImages();

        // sample particles from brain mask
        initParticles();

        // test trainer
        testTrainer();
    }

    void GroupSimul::testTrainer() {
        ParticleTrainer trainer(&_solver->m_System);
        IntVector counts;
        counts.push_back(80);
        counts.push_back(50);
        counts.push_back(50);
        counts.push_back(50);

        trainer.setNumberOfPointSets(2);
        trainer.setNumberOfParticles(counts);
        trainer.sampleParticles();


        // phase #1
        _solver->SpreadMultiParticles();


        trainer.initialClosestCorrespondences();
        trainer.removeOutliers();
        _solver->m_System.RemoveDisabledParticles();

        // phase #2
        _solver->SpreadMultiParticles();


        for (int s = 0; s < _solver->m_System.GetNumberOfSubjects(); s++) {
            trainer.computeDelaunayTriangulation(s, 1, s);
            NeighborWeightMatrix& neighbors = trainer.getPointSets(s);
            
            // index for display
            IntIndex x, y;
            QPainterPath pathItem;
            ParticleSubject& subj = _solver->m_System[s];
            for (int i = 0; i < neighbors.size(); i++) {
                int sourceId = neighbors[i].first;
                Particle& sourceParticle = subj[sourceId];
                ParticleWeightVector& targets = neighbors[i].second;
                subj.ComputeIndexX(sourceParticle, x);
                for (int j = 0; j < targets.size(); j++) {
                    pathItem.moveTo(x[0], x[1]);
                    int targetId = targets[j].first;
                    Particle& targetParticle = subj[targetId];
                    subj.ComputeIndexX(targetParticle, y);
                    pathItem.lineTo(y[0], y[1]);
                }
            }

            QColor penColor = Qt::yellow;
            penColor.setAlphaF(0.1);
            QPen pen(penColor, 1);

            QGraphicsPathItem* graphicsPath = new QGraphicsPathItem();
            graphicsPath->setPath(pathItem);
            graphicsPath->setPen(pen);
            graphicsPath->setParentItem(_imageItems[s]);
        }


        // establish correspondences
        //
//        trainer.establishCorrespondences();
//        trainer.removeOutliers();

        int nSubjs = _solver->m_System.GetNumberOfSubjects();
        ParticleSubject& subj0 = _solver->m_System[0];
        for (int i = 0; i < nSubjs; i++) {
            ParticleSubject& subj = _solver->m_System[i];
            VNLVector laplaceSignature(subj.GetNumberOfPoints());
            for (int j = 0; j < laplaceSignature.size(); j++) {
                laplaceSignature[j] = subj[j].correspondenceScore;
            }
            _particleGroups[i].setScalars(laplaceSignature);
            _particleGroups[i].useScalars(true);
            _particleGroups[i].setReferenceSubject(&subj0);
            _particleGroups[i].updateParticles();

            cout << "Laplace: " << laplaceSignature << endl;
        }

    }

    void GroupSimul::setResolutionLevel() {
        if (_ui->actionLevel0->isChecked()) {
            _solver->m_System.SetCurrentResolutionLevel(0);
        } else if (_ui->actionLevel1->isChecked()) {
            _solver->m_System.SetCurrentResolutionLevel(1);
        } else if (_ui->actionLevel2->isChecked()) {
            _solver->m_System.SetCurrentResolutionLevel(2);
        }
    }

    void GroupSimul::initImages() {
        _images.resize(_nSubjs);
        _warpedImages.resize(_nSubjs);
        _labels.resize(_nSubjs);
        _imageItems.resize(_nSubjs);
        _secondImageItems.resize(_nSubjs);
        _particleGroups.resize(_nSubjs);
        _secondParticleGroups.resize(_nSubjs);

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

            _secondParticleGroups[i].useParticleX(false);
            _secondParticleGroups[i].setListener(this);
            _secondParticleGroups[i].hideParticles(true);
            _secondParticleGroups[i].setScene(&_scene);
            _secondParticleGroups[i].setParentItem(_secondImageItems[i]);
            _secondParticleGroups[i].createParticles(&(_solver->m_System[i]));

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

//            _particleGroups[i].useScalars(true);
//            _particleGroups[i].setScalars(_solver->m_System.NormalizedImageEnergy);

            _secondParticleGroups[i].hideParticles(!on);
            _secondParticleGroups[i].updateParticles();
        }
        showAttribute(_particleGroups[0].getSelectedParticleId());
//        cout << "Image Energy: " << _solver->m_System.ImageEnergy.sum() << endl;
    }

    void GroupSimul::computeWarpedImages() {
        stringstream configstream;
        string config = _ui->config->toPlainText().toStdString();
        configstream.str(config);

        Options options;
        configstream >> options;

        const int n = _solver->m_System.GetNumberOfSubjects();
        for (int i = 0; i < n; i++) {
            // show warpedImage
            if (_solver->m_System[i].m_WarpedImage.IsNotNull()) {
                cout << "showing the cached warped image" << endl;
                _secondImageItems[i]->setImage<RealImage>(_solver->m_System[i].m_WarpedImage);
            }
            
            _secondImageItems[i]->setRange(_imageItems[0]->getRange(0), _imageItems[0]->getRange(1));
            _secondImageItems[i]->refresh();
            _secondImageItems[i]->show();

            _secondParticleGroups[i].updateParticles();
            _secondParticleGroups[i].hideParticles(false);
            
            __realIO.WriteImage(QString("/tmpfs/WarpedImage_%1.nii.gz").arg(i).toStdString(), _warpedImages[i]);
        }
    }

    // compute warped image using current particle display which may be different from the one used during registration
    void GroupSimul::computeParticleWarp() {
        stringstream configstream;
        string config = _ui->config->toPlainText().toStdString();
        configstream.str(config);

        Options options;
        configstream >> options;

        cout << "ControlPointSpacing: " << options.GetInt("ControlPointSpacing:", 8);

        const int n = _solver->m_System.GetNumberOfSubjects();
        _solver->m_System.ComputeXMeanSubject();
        ParticleSubject& meanSubj = _solver->m_System.GetMeanSubject();

        for (int i = 0; i < n; i++) {
            // warp to the first
            ParticleBSpline bspline;
            bspline.SetControlPointSpacing(options.GetInt("ControlPointSpacing:", 8));
            bspline.SetReferenceImage(_labels[i]);
            bspline.EstimateTransform(_solver->m_System[0], _solver->m_System[i]);

            FieldTransformType::Pointer transform = bspline.GetTransform();
            _warpedImages[i] = bspline.WarpImage(_solver->m_System[i].GetImage(0));

            LabelImage::Pointer labelOutput = bspline.WarpLabel(_solver->m_System[i].GetLabel());

            string labelFile = _solver->m_Options.GetStringVectorValue("LabelImages:", i);
            string filename = QString::fromStdString(labelFile).replace("_label", "_warpedLabelTo0").toStdString();
            __labelIO.WriteImage(filename, labelOutput);

            // warp to the mean
            {
                ParticleBSpline bspline;
                bspline.SetControlPointSpacing(options.GetInt("ControlPointSpacing:", 8));
                bspline.SetReferenceImage(_labels[i]);
                bspline.EstimateTransform(meanSubj, _solver->m_System[i]);
                FieldTransformType::Pointer transform = bspline.GetTransform();
                LabelImage::Pointer labelOutput = bspline.WarpLabel(_solver->m_System[i].GetLabel());
                string labelFile = _solver->m_Options.GetStringVectorValue("LabelImages:", i);
                string filename = QString::fromStdString(labelFile).replace("_label", "_warpedLabelToMean").toStdString();
                __labelIO.WriteImage(filename, labelOutput);
            }



            _secondImageItems[i]->setImage<RealImage>(_warpedImages[i]);
            _secondImageItems[i]->generateUserGrids<FieldTransformType>(transform);

            _secondImageItems[i]->setRange(_imageItems[0]->getRange(0), _imageItems[0]->getRange(1));
            _secondImageItems[i]->refresh();
            _secondImageItems[i]->show();

            _secondParticleGroups[i].updateParticles();
            _secondParticleGroups[i].hideParticles(false);
        }
    }

    void GroupSimul::computeImageBspline() {
        RealImage::Pointer image = bsplineRegistration(_solver->m_System[1].GetImage(0), _solver->m_System[0].GetImage(0));
        if (image.IsNotNull()) {
            _secondImageItems[0]->setImage<RealImage>(image);
            _secondImageItems[0]->setRange(_imageItems[0]->getRange(0), _imageItems[0]->getRange(1));
            _secondImageItems[0]->refresh();
            _secondImageItems[0]->show();
        }
    }

    void GroupSimul::showWarpedImages() {
        
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
        _solver->verbose = true;//_ui->verboseOutput->isChecked();


        // prepare simultation
        _solver->t0 = _solver->m_Options.GetRealVectorValue("TimeRange:", 0);
        _solver->t1 = _solver->m_Options.GetRealVectorValue("TimeRange:", 2);
        _solver->dt = _solver->m_Options.GetRealVectorValue("TimeRange:", 1);
        _solver->t = 0;

        // UI status
        _ui->runStepButton->setText("Stop");
        _ui->runStepButton->setChecked(true);
        _ui->loadingAnimation->movie()->start();

        setResolutionLevel();

        m_timer.start(500, this);
        m_future = QtConcurrent::run(this, &GroupSimul::threadedRun);
        m_futureWatcher.setFuture(m_future);
    }


    void GroupSimul::threadedRun() {
        _solver->RunLoopBegin();
        while (_solver->t < _solver->t1 && _solver->continueToRun) {
            _solver->RunStep();
        }
        _solver->RunLoopEnd();
        stopRun();
    }

    void GroupSimul::timerEvent(QTimerEvent* event) {
        showParticles(_ui->actionShowParticles->isChecked());
        if (!_solver->continueToRun) {
            m_timer.stop();
        }
    }

    void GroupSimul::stopRun() {
        _ui->loadingAnimation->movie()->stop();
        _ui->runStepButton->setText("Step");
        _ui->runStepButton->setChecked(false);
        _solver->continueToRun = false;
        disconnect(&m_futureWatcher, SIGNAL(finished()), this, SLOT(stopRun()));
        disconnect(&m_futureWatcher, SIGNAL(canceled()), this, SLOT(stopRun()));
        m_future = QFuture<void>();
    }

    void GroupSimul::mousePressed(QGraphicsEllipseEventItem *sender, QGraphicsSceneMouseEvent *event) {
        event->accept();

        int particleId = sender->data(0).value<int>();

        if (_solver->m_System.NormalizedImageEnergy.size() > particleId) {
            cout << "Particle #" << particleId << " Energy: " << _solver->m_System.NormalizedImageEnergy[particleId] << endl;
        }
    }

    void GroupSimul::mouseReleased(QGraphicsEllipseEventItem *sender, QGraphicsSceneMouseEvent *event) {
        int particleId = sender->data(0).value<int>();
        for (int i = 0; i < _particleGroups.size(); i++) {
            _particleGroups[i].selectParticle(particleId);
        }

        showAttribute(particleId);

        QPointF pos = sender->mapToParent(event->pos().x(), event->pos().y());

        RealImage::PointType idx, point;
        idx[0] = pos.x();
        idx[1] = pos.y();
        
        ParticleSystem& system = _solver->m_System;
        for (int i = 0; i < system.GetNumberOfSubjects(); i++) {
            if (_particleGroups[i].getItem(particleId) != sender) {
                continue;
            }
            system[i].ContinuousIndexToPoint(idx, point);
            fordim (k) {
                system[i][particleId].x[k] = point[k];
            }
            _particleGroups[i].updateParticles();
        }
    }

    void GroupSimul::createDistanceMap() {
        LabelImage::Pointer image = _solver->m_System[0].GetLabel();
        ImageProcessing proc;
        VectorImage::Pointer dmap = proc.ComputeDistanceMap(image);
        RealImage::Pointer dmapMag = proc.ComputeMagnitudeMap(dmap);
        __realIO.WriteImage("/tmpfs/dmap.nii.gz", dmapMag);


    }

    void GroupSimul::showAttribute(int particleId) {
        if (particleId < 0) {
            return;
        }
        _scene2.clear();
        const int w = ATTR_SIZE;
        int j = 0;
        cout << "F: " << endl;
        for (int i = 0; i < _solver->m_System.GetNumberOfSubjects(); i++) {
            QGraphicsRealImageItem* item = new QGraphicsRealImageItem();
            ParticleAttribute* attr = _solver->intensityForce.GetAttribute(i, particleId);
            if (attr) {
                item->setImage(attr->x, w, w);
//                item->setRange(_imageItems[0]->getRange(0), _imageItems[0]->getRange(1));
                item->computeRange();
                item->setFlip(QGraphicsRealImageItem::UpDown);
                item->refresh();
                item->setPos((w+3)*i, (w+3)*j);

                _scene2.addItem(item);
                _ui->graphicsView2->fitInView(0, 0, (w+3)*(i+1), (w+3)*(j+1), Qt::KeepAspectRatio);
//                _ui->graphicsView2->fitInView(0, 0, (w+3), w, Qt::KeepAspectRatio);
                cout << attr->F[0] << "," << attr->F[1] << "; ";
            }
            cout << _solver->m_System[i][particleId].correspondenceScore << "; ";
        }
        cout << endl;
    }

    void GroupSimul::showAttributes() {
        _scene2.clear();
        const int w = ATTR_SIZE;
        for (int j = 0; j < _solver->m_System.GetNumberOfParticles(); j++) {
            for (int i = 0; i < _solver->m_System.GetNumberOfSubjects(); i++) {
                QGraphicsRealImageItem* item = new QGraphicsRealImageItem();
                ParticleAttribute* attr = _solver->intensityForce.GetAttribute(i, j);
                if (attr) {
                    item->setImage(attr->x, w, w);
                    item->setRange(_imageItems[0]->getRange(0), _imageItems[0]->getRange(1));
                    item->setFlip(QGraphicsRealImageItem::UpDown);
                    item->refresh();
                    item->setPos((w+3)*j, (w+3)*i);

                    QGraphicsDirectionItem* dirItem = new QGraphicsDirectionItem(item);
                    dirItem->setPen(QPen(Qt::yellow));
                    dirItem->setDirection(QPointF(attr->F[0], attr->F[1]));

                    _scene2.addItem(item);
                    _ui->graphicsView2->fitInView(0, 0, (w+3)*(j+1), (w+3)*(i+1), Qt::KeepAspectRatio);

                    cout << "F = [" << attr->F[0] << "," << attr->F[1] << "]" << endl;
                }
            }
        }
    }

    void GroupSimul::saveTransform() {
        ParticleSystem& system = _solver->m_System;
        for (int i = 0; i < system.GetNumberOfSubjects(); i++) {
            if (system[i].m_InverseDeformableTransform.IsNotNull()) {
                string labelFile = _solver->m_Options.GetStringVectorValue("LabelImages:", i);
                string filename = QString::fromStdString(labelFile).replace("_label", "_warpedLabel").toStdString();
                
                LabelImage::Pointer label = system[i].WarpLabelToMeanSpace();
                __labelIO.WriteImage(filename, label);
            }
        }
    }
}
