//
//  mainwindow.cpp
//  laplacePDE
//
//  Created by Joohwi Lee on 11/13/12.
//
//

#include "mainwindow.h"
#include "QFileDialog"
#include "vtkRenderer.h"
#include "vtkGenericOpenGLRenderWindow.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "QVTKInteractor.h"
#include "myParticleAlgorithm.h"
#include <vtkBoxMuellerRandomSequence.h>
#include <vtkMinimalStandardRandomSequence.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include "myImageParticlesAlgorithm.h"
#include "itkARGBColorFunction.h"

static int g_AnimFrame = 0;
static bool g_showTraceParticles = false;

ParametersVectorType g_Params;
ImageParticlesAlgorithm::Pointer g_imageParticlesAlgo;

MainWindow::MainWindow(QWidget* parent): m_ParticleColors(this), m_Props(this) {
    ui.setupUi(this);

    m_Renderer = vtkRenderer::New();
    m_PropScene.SetRenderer(m_Renderer);

    vtkGenericOpenGLRenderWindow* renWin = ui.qvtkwidget->GetRenderWindow();
    renWin->AddRenderer(m_Renderer);
    m_Interactor = ui.qvtkwidget->GetInteractor();

    ui.graphicsView->setScene(&m_scene);
    ui.graphicsView->setBackgroundBrush(QBrush(Qt::black, Qt::SolidPattern));

    QObject::connect(&m_Timer, SIGNAL(timeout()), this, SLOT(on_animationTimeout()));
    QObject::connect(ui.sliceIndex, SIGNAL(sliderMoved(int)), this, SLOT(updateScene()));
    QObject::connect(ui.showXY, SIGNAL(toggled(bool)), this, SLOT(chooseSlice()));
    QObject::connect(ui.showYZ, SIGNAL(toggled(bool)), this, SLOT(chooseSlice()));
    QObject::connect(ui.showZX, SIGNAL(toggled(bool)), this, SLOT(chooseSlice()));
    QObject::connect(ui.showSource, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showTarget, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showSourceLabel, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showTargetLabel, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.zoomRatio, SIGNAL(valueChanged(double)), this, SLOT(updateScene()));
    QObject::connect(ui.labelOpacity, SIGNAL(valueChanged(int)), this, SLOT(updateScene()));
    QObject::connect(&m_ParticleColors, SIGNAL(triggered(QAction*)), this, SLOT(updateScene()));

	ui.graphicsView->setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform);
    m_ParticleColors.addAction(ui.actionParticleBlack);
    m_ParticleColors.addAction(ui.actionParticleWhite);
    m_ParticleColors.addAction(ui.actionParticleRed);
    m_ParticleColors.addAction(ui.actionParticleGreen);
    m_ParticleColors.addAction(ui.actionParticleBlue);
    m_ParticleColors.addAction(ui.actionParticleHSV);
}

MainWindow::~MainWindow() {

}

void MainWindow::ReadyToExperiments() {
    LoadImage("/data/00.T2.nrrd", 0);
    LoadLabel("/data/00.Label.nrrd", 0);
    on_actionRandomParticlesInit_triggered();
}

void MainWindow::LoadImage(QString fileName, int idx) {
    ImageContainer::Pointer image;
    if (m_ImageList.size() > idx) {
        image = m_ImageList[idx];
    } else {
        image = ImageContainer::New();
        for (int j = 0; j < idx; j++) {
            if (image.IsNull()) {
                m_ImageList.push_back(ImageContainer::Pointer());
            }
        }
        m_ImageList.push_back(image);
    }
    image->LoadImage(fileName.toUtf8().data());
    if (image->HasImage() && !image->HasLabel()) {
        ui.sliceIndex->setMaximum(m_ImageList[idx]->GetSize()[0]);
        ui.sliceIndex->setValue(m_ImageList[idx]->GetSliceIndex()[0]);
        updateScene();
    }
    ui.tabWidget->setCurrentWidget(ui.imageTab);
}

void MainWindow::LoadLabel(QString fileName, int idx) {
    ImageContainer::Pointer image;
    if (m_ImageList.size() > idx) {
        image = m_ImageList[idx];
    } else {
        image = ImageContainer::New();
        for (int j = 0; j < idx; j++) {
            if (image.IsNull()) {
                m_ImageList.push_back(ImageContainer::Pointer());
            }
        }
        m_ImageList.push_back(image);
    }
    image->LoadLabel(fileName.toUtf8().data());
    if (!image->HasImage() && image->HasLabel()) {
        ui.sliceIndex->setMaximum(m_ImageList[idx]->GetSize()[0]);
        ui.sliceIndex->setValue(m_ImageList[idx]->GetSliceIndex()[0]);
    }
    updateScene();
    ui.tabWidget->setCurrentWidget(ui.imageTab);
}

void MainWindow::on_actionOpenSource_triggered() {
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Volumes (*.nrrd *.nii *.gipl.gz)"));
    if (fileName.isNull()) {
        return;
    }
    LoadImage(fileName, 0);
}

void MainWindow::on_actionOpenSourceLabel_triggered() {
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Volumes (*.nrrd *.nii *.gipl.gz)"));
    if (fileName.isNull()) {
        return;
    }
    ui.showSourceLabel->setChecked(true);
    LoadLabel(fileName, 0);
}


void MainWindow::on_actionOpenTarget_triggered() {
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Volumes (*.nrrd *.nii *.gipl.gz)"));
    if (fileName.isNull()) {
        return;
    }
    LoadImage(fileName, 1);
}


void MainWindow::on_actionOpenSurface_triggered() {
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Volumes (*.vtk)"));
    if (fileName.isNull()) {
        return;
    }

    vtkPolyData* poly = m_PropScene.LoadPolyData(fileName.toUtf8().data());
    m_PropScene.AddPolyData("mainSurface", poly);
    m_PropScene.SetColor(1, 0, 0);
    m_PropScene.SetRepresentation(0);

    m_Renderer->ResetCamera();
    m_Interactor->Render();

    ParticleAlgorithm::Pointer algo = ParticleAlgorithm::New();
    algo->SetPropertyAccess(m_Props);
    algo->SetInitialParticles(poly);
    algo->RunOptimization();
    g_Params = algo->GetParameterTrace();
    ////
    //    vtkPolyData* poly2 = poly->NewInstance();
    //    poly2->DeepCopy(poly);
    //
    //    vtkPointSet* result = algo->GetResultPoints();
    //    poly2->SetPoints(result->GetPoints());
}

void MainWindow::chooseSlice() {
    int sliceView = 0;
    if (ui.showXY->isChecked()) {
        sliceView = 0;
    } else if (ui.showYZ->isChecked()) {
        sliceView = 1;
    } else if (ui.showZX->isChecked()) {
        sliceView = 2;
    }
    for (int i = 0; i < m_ImageList.size(); i++) {
        m_ImageList[i]->SetSliceDir(sliceView);
    }
    updateScene();
    g_imageParticlesAlgo = ImageParticlesAlgorithm::Pointer(NULL);
}

void MainWindow::updateScene() {
    m_scene.clear();


    QTransform viewTransform;
    viewTransform.scale(m_Props.GetDouble("zoomRatio", 1), m_Props.GetDouble("zoomRatio", 1));
    ui.graphicsView->setTransform(viewTransform);

    int dim = GetCurrentView();
    int image = GetCurrentImage();

    if (IsImageAvailable(image) || IsLabelAvailable(image)) {
        ui.sliceIndex->setMaximum(m_ImageList[image]->GetSize()[dim]);
        m_ImageList[image]->SetSliceIndex(dim, ui.sliceIndex->value());
    }
    if (IsImageAvailable(image)) {
        QGraphicsPixmapItem* item = m_scene.addPixmap(m_ImageList[image]->GetPixmap(dim));
        ui.graphicsView->centerOn((const QGraphicsItem*) item);

    }
    if (ui.showSourceLabel->isChecked()) {
        if (IsLabelAvailable(0)) {
            m_ImageList[0]->SetLabelAlpha(m_Props.GetInt("labelOpacity", 128));
            m_scene.addPixmap(m_ImageList[0]->GetLabelPixmap(dim));
        }
    }
    if (ui.showTargetLabel->isChecked()) {
        if (IsLabelAvailable(1)) {
            m_ImageList[1]->SetLabelAlpha(m_Props.GetInt("labelOpacity", 128));
            m_scene.addPixmap(m_ImageList[1]->GetLabelPixmap(dim));
        }
    }
    if (g_imageParticlesAlgo.IsNotNull() && ui.actionShowParticles->isChecked()) {
        const OptimizerParametersType* particles = NULL;
        if (g_showTraceParticles) {
            particles = g_imageParticlesAlgo->GetTraceParameters(g_AnimFrame);
        } else {
            particles = &(g_imageParticlesAlgo->GetCurrentParams());
        }
        if (particles != NULL) {
//            const int nSubj = g_imageParticlesAlgo->GetNumberOfSubjects();
            const int nPoints = g_imageParticlesAlgo->GetNumberOfPoints();
            const int nParams = g_imageParticlesAlgo->GetNumberOfParams();
            const int nDims = 2;

            itk::Function::HSVColormapFunction<double, itk::RGBAPixel<unsigned char> >::Pointer colorFunc = itk::Function::HSVColormapFunction<double, itk::RGBAPixel<unsigned char> >::New();
            colorFunc->SetMinimumInputValue(0);
            colorFunc->SetMaximumInputValue(nPoints);

            for (int i = 0; i < nPoints; i++) {
                double x = particles->GetElement(image * nParams + nDims * i);
                double y = particles->GetElement(image * nParams + nDims * i + 1);

                QColor pointColor = Qt::black;
                if (ui.actionParticleRed->isChecked()) {
                    pointColor = Qt::red;
                } else if (ui.actionParticleGreen->isChecked()) {
                    pointColor = Qt::green;
                } else if (ui.actionParticleBlue->isChecked()) {
                    pointColor = Qt::blue;
                } else if (ui.actionParticleWhite->isChecked()) {
                    pointColor = Qt::white;
                } else if (ui.actionParticleHSV->isChecked()) {
                    itk::RGBAPixel<unsigned char> rgb = colorFunc->operator()(i);
                    pointColor = QColor(rgb[0], rgb[1], rgb[2]);
                }

                m_scene.addEllipse(x, y, 1, 1, QPen(pointColor),
                                   QBrush(pointColor, Qt::SolidPattern));
            }
        }
    }
}


void MainWindow::on_actionAnimation_triggered() {
    g_AnimFrame = 0;
    m_Timer.start(m_Props.GetInt("animationTimeout", 100));
}

void MainWindow::on_animationTimeout() {
    if (ui.tabWidget->currentWidget() == ui.imageTab) {
        // animation should play image particels if imageTab is showing
        g_showTraceParticles = true;
        if (g_imageParticlesAlgo.IsNotNull() && g_AnimFrame >= g_imageParticlesAlgo->GetNumberOfTraces()) {
            g_showTraceParticles = false;
            m_Timer.stop();
            return;
        }
        // the particle movement is handeled in updateScene() with g_showTraceParticles flag
        // this is not an optimal implementation because this will update every images together
        updateScene();
        ui.statusbar->showMessage(QString("%1 frame played...").arg(g_AnimFrame));        
        g_AnimFrame++;
    } else if (ui.tabWidget->currentWidget() == ui.modelTab) {
        vtkPolyData* poly = m_PropScene.FindPolyData("mainSurface");
        if (poly == NULL) {
            m_Timer.stop();
            return;
        }

        if (g_AnimFrame < (int) g_Params.size()) {
            OptimizerParametersType param = g_Params[g_AnimFrame];
            if (param.GetSize() != 3 * poly->GetNumberOfPoints()) {
                m_Timer.stop();
                return;
            }

            int nPoints = poly->GetNumberOfPoints();
            for (int i = 0; i < nPoints; i++) {
                poly->GetPoints()->SetPoint(i, param[3*i], param[3*i+1], param[3*i+2]);
            }
            m_PropScene.ModifyLastActor();

            ui.statusbar->showMessage(QString("%1 frame showing").arg(++g_AnimFrame));
            m_Interactor->Render();
        }
    }
}

void MainWindow::on_actionRandomParticlesInit_triggered() {
    g_imageParticlesAlgo = ImageParticlesAlgorithm::New();
    g_imageParticlesAlgo->SetPropertyAccess(m_Props);
    g_imageParticlesAlgo->SetViewingDimension(GetCurrentView());
    g_imageParticlesAlgo->SetImageList(&m_ImageList);
    g_imageParticlesAlgo->CreateRandomInitialPoints(m_Props.GetInt("numberOfPoints", 100));
    updateScene();
}


void MainWindow::on_actionRunImageParticles_triggered() {
    if (IsImageAvailable(0)) {
        if (g_imageParticlesAlgo.IsNotNull()) {
            g_imageParticlesAlgo->SetImageList(&m_ImageList);
            // g_imageParticlesAlgo->SetPropertyAccess(m_Props);
            g_imageParticlesAlgo->RunOptimization();
        }
    }
    updateScene();
}