//
//  mainwindow.cpp
//  laplacePDE
//
//  Created by Joohwi Lee on 11/13/12.
//
//

#include "mainwindow.h"
#include "QFileDialog"
#include "QMessageBox"
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
#include "bitset"
#include "armadillo"

static int g_AnimFrame = 0;
static bool g_showTraceParticles = false;

surface::ParametersVectorType g_Params;
ImageParticlesAlgorithm::Pointer g_imageParticlesAlgo;
surface::ParticleAlgorithm::Pointer g_ParticleAlgo;
my::ImplicitSurfaceConstraint g_constraint;


MainWindow::MainWindow(QWidget* parent): m_ParticleColors(this), m_Props(this), m_CompareDialog(this), m_BSplineVisDialog(this) {
    ui.setupUi(this);

    ui.toolBar->addWidget(ui.zoomLabel);
    ui.toolBar->addWidget(ui.zoomRatio);
    
    m_Renderer = vtkRenderer::New();
    m_PropScene.SetRenderer(m_Renderer);

    vtkGenericOpenGLRenderWindow* renWin = ui.qvtkwidget->GetRenderWindow();
    renWin->AddRenderer(m_Renderer);
    m_Interactor = ui.qvtkwidget->GetInteractor();

    ui.graphicsView->setScene(&m_scene);
    ui.graphicsView->setBackgroundBrush(QBrush(Qt::black, Qt::SolidPattern));


    QObject::connect(&m_Timer, SIGNAL(timeout()), this, SLOT(on_animationTimeout()));
    QObject::connect(ui.sliceIndex, SIGNAL(sliderMoved(int)), this, SLOT(chooseSlice()));
    QObject::connect(ui.showXY, SIGNAL(toggled(bool)), this, SLOT(chooseSlice()));
    QObject::connect(ui.showYZ, SIGNAL(toggled(bool)), this, SLOT(chooseSlice()));
    QObject::connect(ui.showZX, SIGNAL(toggled(bool)), this, SLOT(chooseSlice()));
    QObject::connect(ui.grayImages, SIGNAL(currentIndexChanged(int)), this, SLOT(selectImage(int)));
    QObject::connect(ui.labelImages, SIGNAL(currentIndexChanged(int)), this, SLOT(selectLabel(int)));

    QObject::connect(ui.actionShowCorrespondence, SIGNAL(triggered()), this, SLOT(updateScene()));
    QObject::connect(ui.actionShowDisplacement, SIGNAL(triggered()), this, SLOT(updateScene()));

    QObject::connect(ui.zoomRatio, SIGNAL(valueChanged(double)), this, SLOT(updateScene()));
    QObject::connect(ui.labelOpacity, SIGNAL(valueChanged(int)), this, SLOT(updateScene()));
    QObject::connect(&m_ParticleColors, SIGNAL(triggered(QAction*)), this, SLOT(updateScene()));
    QObject::connect(ui.actionShowParticles, SIGNAL(triggered()), this, SLOT(updateScene()));
    QObject::connect(ui.actionShowBSplineControlPoints, SIGNAL(triggered()), this, SLOT(updateScene()));
    QObject::connect(ui.showGray, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showLabel, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.showDerived, SIGNAL(toggled(bool)), this, SLOT(updateScene()));
    QObject::connect(ui.particlesInitialization, SIGNAL(clicked()), this, SLOT(on_actionRandomParticlesInit_triggered()));
    QObject::connect(ui.selectedPoints, SIGNAL(textChanged()), this, SLOT(updateScene()));
    
	ui.graphicsView->setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform);
    m_ParticleColors.addAction(ui.actionParticleBlack);
    m_ParticleColors.addAction(ui.actionParticleWhite);
    m_ParticleColors.addAction(ui.actionParticleRed);
    m_ParticleColors.addAction(ui.actionParticleGreen);
    m_ParticleColors.addAction(ui.actionParticleBlue);
    m_ParticleColors.addAction(ui.actionParticleHSV);

    ui.costPlot->setColor(QColor(0xf5, 0xf3, 0xff));
    ui.costPlot->addGraph();
}

MainWindow::~MainWindow() {

}

void MainWindow::showEvent(QShowEvent* event) {
    on_actionOpenBSplineVis_triggered();
}

void MainWindow::on_actionOpenCompareWindow_triggered() {
    m_CompareDialog.show();
}

void MainWindow::on_actionOpenBSplineVis_triggered() {
    m_BSplineVisDialog.SetImageParticlesAlgorithm(g_imageParticlesAlgo);
    m_BSplineVisDialog.show();
}


void MainWindow::on_actionLoadParticles_triggered() {
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), "/data/Particles", tr("Armadillo matrix (*.txt)"));

    if (fileName.isNull()) {
        return;
    }
    arma::vec armParams;
    if (armParams.quiet_load(fileName.toStdString(), arma::arma_ascii)) {
        VNLVector params(armParams.memptr(), armParams.size());
        g_imageParticlesAlgo->SetCurrentParams(params);
        cout << "Particles loaded: " << armParams.size() << endl;
        updateScene();
    } else {
        QMessageBox::warning(this, "Particle Load Fail", fileName + " not available");
    }
}

void MainWindow::on_actionSaveParticles_triggered() {
    QString fileName = QFileDialog::getSaveFileName(this, tr("Open File"), "/data/Particles", tr("Armadillo matrix (*.txt)"));
    
    if (fileName.isNull()) {
        return;
    }
    
    VNLVector params(g_imageParticlesAlgo->GetCurrentParams());
    arma::vec armParams(params.data_block(), params.size(), false);
    if (!armParams.quiet_save(fileName.toStdString(), arma::arma_ascii)) {
        QMessageBox::warning(this, "Particle Save Fail", fileName + " not available");
    }
}

void MainWindow::on_actionRunOnIntersection_triggered() {
    g_imageParticlesAlgo->CreateUniformInitialization();
    updateScene();
}


void MainWindow::on_derivedImages_currentIndexChanged(int n) {
    if (n > 0) {
        ui.showDerived->setCheckState(Qt::Checked);
    } else {
        ui.showDerived->setCheckState(Qt::Unchecked);
    }
    updateScene();
}

void MainWindow::EventRaised(int eventId, int eventCode, const void* src, void* data) {
    if (eventId == 0xADDCE) {
        ui.derivedImages->clear();
        ui.derivedImages->addItem("-----");
        ImageContainer::StringList derivedNames;
        ImageContainer::GetDerivedViewNames(derivedNames);
        for (int i = 0; i < derivedNames.size(); i++) {
            ui.derivedImages->addItem(derivedNames[i].c_str());
        }
    } else if (eventId == 0xADDCEC) {
        double* doubleData = (double*) data;
        ui.costPlot->graph()->addData(doubleData[0], doubleData[1]);
        ui.costPlot->rescaleAxes();
        ui.costPlot->xAxis->setLabel("iteration");
        ui.costPlot->yAxis->setLabel("cost");
        ui.costPlot->replot();
    }
}

void MainWindow::ReadyToExperiments() {
//    LoadImage("/data/Particles/00.T2.half.nrrd");
//    LoadLabel("/data/Particles/00.Label.half.nrrd");
//    LoadImage("/data/Particles/16.T2.half.nrrd");
//    LoadLabel("/data/Particles/16.Label.half.nrrd");
//    LoadImage("/data/Particles/cs1_tex.nrrd");
//    LoadLabel("/data/Particles/cs1.nrrd");
//    LoadImage("/data/Particles/cs2_tex.nrrd");
//    LoadLabel("/data/Particles/cs2.nrrd");
//    LoadSurface("/data/Particles/00.vtk");
//    LoadImage("/data/Particles/Image001.nrrd");
//    LoadLabel("/data/Particles/Image001.nrrd");
//    LoadImage("/data/Particles/Image002.nrrd");
//    LoadLabel("/data/Particles/Image002.nrrd");
    LoadImage("/data/Particles/image_circle.nrrd");
    LoadLabel("/data/Particles/image_circle_label.nrrd");
    LoadImage("/data/Particles/image_square.nrrd");
    LoadLabel("/data/Particles/image_square_label.nrrd");
//    LoadImage("/data/Particles/image_black_gray.nrrd");
//    LoadLabel("/data/Particles/image_center_mask.nrrd");

    on_actionRandomParticlesInit_triggered();
    g_constraint.SetImageList(&m_ImageList);
}

ImageParticlesAlgorithm::Pointer MainWindow::GetImageParticlesAlgorithm() {
    return g_imageParticlesAlgo;
}

void MainWindow::selectImage(int idx) {
    updateScene();
}

void MainWindow::selectLabel(int idx) {
    updateScene();
}

void MainWindow::LoadImage(QString fileName) {
    ImageContainer::Pointer image(NULL);
    for (int i = 0; i < (int) m_ImageList.size(); i++) {
        if (!m_ImageList[i]->HasImage()) {
            image = m_ImageList[i];
        }
    }
//    cout << "Loading 1: " << fileName.toStdString() << endl;
    if (image.IsNull()) {
//        cout << "Image is null" << endl;
        image = ImageContainer::New();
//        cout << "ImageContainer::New()" << endl;
        image->SetEventCallback(this);
//        cout << "SetEventCallback(this);" << endl;
        m_ImageList.push_back(image);
    }
//    cout << "Loading 2: " << fileName.toStdString() << endl;
    image->LoadImage(fileName.toUtf8().data());

    if (!image->HasLabel()) {
        ui.sliceIndex->setMaximum(image->GetSize()[0]);
        ui.sliceIndex->setValue(image->GetSliceIndex()[0] > 1 ? image->GetSliceIndex()[0] / 2 : 0);
    }
//    cout << "Loading 3: " << fileName.toStdString() << endl;
    updateScene();
    ui.tabWidget->setCurrentWidget(ui.imageTab);
    ui.toolBox->setCurrentWidget(ui.imageSettings);
    ui.grayImages->addItem(fileName);
    g_imageParticlesAlgo = ImageParticlesAlgorithm::Pointer(NULL);
}

/**
 * Load label image and add to ImageList
 * Handles label selection combobox
 */
void MainWindow::LoadLabel(QString fileName) {
    ImageContainer::Pointer image(NULL);
    for (int i = 0; i < (int) m_ImageList.size(); i++) {
        if (!m_ImageList[i]->HasLabel()) {
            image = m_ImageList[i];
            break;
        }
    }
    if (image.IsNull()) {
        image = ImageContainer::New();
        image->SetEventCallback(this);
        m_ImageList.push_back(image);
    }

    image->LoadLabel(fileName.toUtf8().data());
    if (!image->HasImage()) {
        ui.sliceIndex->setMaximum(image->GetSize()[0]);
        ui.sliceIndex->setValue(image->GetSliceIndex()[0]);
    }
    updateScene();
    ui.tabWidget->setCurrentWidget(ui.imageTab);
    ui.toolBox->setCurrentWidget(ui.imageSettings);
    ui.labelImages->addItem(fileName);
    g_imageParticlesAlgo = ImageParticlesAlgorithm::Pointer(NULL);
}

void MainWindow::LoadSurface(QString fileName) {
    vtkPolyData* poly = m_PropScene.LoadPolyData(fileName.toUtf8().data());
    m_PropScene.AddPolyData("mainSurface", poly);
    m_PropScene.SetColor(1, 0, 0);
    m_PropScene.SetRepresentation(0);

    m_Renderer->ResetCamera();
    m_Interactor->Render();
    

}


void MainWindow::on_actionAddImage_triggered() {
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Volumes (*.nrrd *.nii *.gipl.gz)"));
    if (fileName.isNull()) {
        return;
    }
    LoadImage(fileName);
}

void MainWindow::on_actionAddLabel_triggered() {
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Volumes (*.nrrd *.nii *.gipl.gz)"));
    if (fileName.isNull()) {
        return;
    }
    LoadLabel(fileName);
}


void MainWindow::on_actionOpenSurface_triggered() {
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/tmpfs/data", tr("Volumes (*.vtk)"));
    if (fileName.isNull()) {
        return;
    }

    LoadSurface(fileName);
//    ParticleAlgorithm::Pointer algo = ParticleAlgorithm::New();
//    algo->SetPropertyAccess(m_Props);
//    algo->SetInitialParticles(poly);
//    algo->RunOptimization();
//    g_Params = algo->GetParameterTrace();
    ////
    //    vtkPolyData* poly2 = poly->NewInstance();
    //    poly2->DeepCopy(poly);
    //
    //    vtkPointSet* result = algo->GetResultPoints();
    //    poly2->SetPoints(result->GetPoints());
}

void MainWindow::on_actionSurfaceSmoothing_triggered() {
    // surface smoothing assumes the first image and label is matched with the current surface model
    // apply the surface entropy minimization smoothing with labelmap restriction
    g_ParticleAlgo = surface::ParticleAlgorithm::New();
    g_ParticleAlgo->SetPropertyAccess(m_Props);
    g_ParticleAlgo->SetImageList(&m_ImageList);
    g_ParticleAlgo->SetModelScene(&m_PropScene);
    g_ParticleAlgo->RunOptimization();
}

void MainWindow::on_actionSurfaceSmoothingContinue_triggered() {
    if (g_ParticleAlgo.IsNull()) {
        return;
    }
}

void MainWindow::on_actionTest_triggered() {
    if (ui.tabWidget->currentWidget() == ui.modelTab) {
        vtkPolyData* object = m_PropScene.CreatePlane(10, 10, 0, 0);
        object->Print(cout);
        m_PropScene.AddPolyData("object", object);
        m_PropScene.SetColor(0, 0, 1);
        m_PropScene.SetRepresentation(1);


        m_Renderer->ResetCamera();
        m_Interactor->Render();
    } else {


        OptimizerParametersType p;
        p.SetSize(100);
        for (int i = 0; i < p.GetSize(); i++) {
            p[i] = i;
        }
        vnl_matrix<double> mat(10, 10);
        mat.set(p.data_block());
        cout << mat << endl;
    }
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
    // manage current view status
    int dim = GetCurrentView();
    int image = GetCurrentImage();
    int labelIdx = ui.labelImages->currentIndex();
    ImageContainer::SetCurrentView(dim);
    ImageContainer::SetCurrentImage(image);
    ImageContainer::SetCurrentLabel(labelIdx);
    // g_constraint.SetImageList(&m_ImageList);
    g_imageParticlesAlgo = ImageParticlesAlgorithm::Pointer(NULL);
}

void MainWindow::updateScene() {
    // manage current view status
    int dim = GetCurrentView();
    int image = GetCurrentImage();
    int labelIdx = ui.labelImages->currentIndex();

    const int nImages = m_ImageList.size();

    ImageContainer::SetCurrentView(dim);
    ImageContainer::SetCurrentImage(image);
    ImageContainer::SetCurrentLabel(labelIdx);

    // remove previous drawings
    m_scene.clear();
    QTransform viewTransform;

    // flip y-axis if viewing plane is XY or YZ
    if (dim == 0 || dim == 1) {
        viewTransform.scale(m_Props.GetDouble("zoomRatio", 1), -m_Props.GetDouble("zoomRatio", 1));

    } else {
        viewTransform.scale(m_Props.GetDouble("zoomRatio", 1), m_Props.GetDouble("zoomRatio", 1));

    }
    ui.graphicsView->setTransform(viewTransform);

    // 

    // Gray image rendering
    if (ui.derivedImages->currentIndex() > 0 && ui.showDerived->isChecked()) {
        
        m_scene.addPixmap(ImageContainer::GetDerivedViewPixmap(ui.derivedImages->currentText().toUtf8().data()));
    } else {
        if (m_ImageList.size() > image && ui.showGray->isChecked()) {
            ui.sliceIndex->setMaximum(m_ImageList[image]->GetSize()[dim]);
            for (int i = 0; i < m_ImageList.size(); i++) {
                m_ImageList[i]->SetSliceIndex(dim, ui.sliceIndex->value());
            }
            QGraphicsPixmapItem* item = m_scene.addPixmap(m_ImageList[image]->GetPixmap(dim));
            ui.graphicsView->centerOn((const QGraphicsItem*) item);
        }

    }

    // Label image rendering
    if (m_ImageList.size() > labelIdx && ui.showLabel->isChecked()) {
        m_ImageList[labelIdx]->SetLabelAlpha(m_Props.GetInt("labelOpacity", 128));
        m_scene.addPixmap(m_ImageList[labelIdx]->GetLabelPixmap(dim));
    }

    std::bitset<1000> selectedPointIds;
    QStringList selectedPoints = ui.selectedPoints->toPlainText().split(" ", QString::SkipEmptyParts);
    bool showSelectedPoints = selectedPoints.size() > 0;
    for (int i = 0; i < selectedPoints.size(); i++) {
        if (selectedPoints[i].toInt() < selectedPointIds.size()) {
            selectedPointIds.set(selectedPoints[i].toInt(), true);
        }
    }

    // check if particles can be rendered
    if (g_imageParticlesAlgo.IsNotNull() && g_imageParticlesAlgo->IsCurrentSliceAndView(dim, ui.sliceIndex->value())) {
        const VNLVector* particles = NULL;
        if (g_showTraceParticles) {
            particles = g_imageParticlesAlgo->GetTraceParameters(g_AnimFrame);
        } else {
            particles = &(g_imageParticlesAlgo->GetCurrentParams());
        }

        if (ui.actionShowParticles->isChecked()) {
            if (particles != NULL) {
                // const int nSubj = g_imageParticlesAlgo->GetNumberOfSubjects();
                const int nPoints = g_imageParticlesAlgo->GetNumberOfPoints();
                const int nParams = g_imageParticlesAlgo->GetNumberOfParams();
                const int nDims = 2;

                // show correspondence
                if (ui.actionShowCorrespondence->isChecked()) {
                    itk::Function::HSVColormapFunction<double, itk::RGBAPixel<unsigned char> >::Pointer colorFunc = itk::Function::HSVColormapFunction<double, itk::RGBAPixel<unsigned char> >::New();
                    colorFunc->SetMinimumInputValue(0);
                    colorFunc->SetMaximumInputValue(nPoints);

                    for (int n = 0; n < nImages; n++) {
                        if (n > 0) {
                            VNLMatrixRef src(nPoints, 2, (double*) &(particles->operator[](0)));
                            VNLMatrixRef dst(nParams, 2, (double*) &(particles->operator[](n*nParams)));

                            for (int i = 0; i < nPoints; i++) {
                                m_scene.addLine(src[i][0], src[i][1], dst[i][0], dst[i][1], QPen(Qt::green));
                            }
                        }
                        for (int i = 0; i < nPoints; i++) {
                            double x = particles->get(n * nParams + nDims * i);
                            double y = particles->get(n * nParams + nDims * i + 1);

                            // better if we have random color palette
                            QColor pointColor = Qt::black;
                            if (n == 0) {
                                pointColor = Qt::red;
                            } else if (n == 1) {
                                pointColor = Qt::blue;
                            } else if (n == 2) {
                                pointColor = Qt::yellow;
                            }
                            m_scene.addEllipse(x-.5, y-.5, 1, 1, QPen(pointColor), QBrush(pointColor, Qt::SolidPattern));
                        }
                    }
                } else if (showSelectedPoints) {
                    itk::Function::HSVColormapFunction<double, itk::RGBAPixel<unsigned char> >::Pointer colorFunc = itk::Function::HSVColormapFunction<double, itk::RGBAPixel<unsigned char> >::New();
                    colorFunc->SetMinimumInputValue(0);
                    colorFunc->SetMaximumInputValue(m_ImageList.size());

                    for (int i = 0; i < nPoints; i++) {
                        if (!selectedPointIds[i]) {
                            continue;
                        }
                        VNLMatrix xyVector(m_ImageList.size(), nDims);
                        for (int n = 0; n < m_ImageList.size(); n++) {
                            double x = particles->get(n * nParams + nDims * i);
                            double y = particles->get(n * nParams + nDims * i + 1);
                            xyVector[n][0] = x;
                            xyVector[n][1] = y;

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
                                itk::RGBAPixel<unsigned char> rgb = colorFunc->operator()(n);
                                pointColor = QColor(rgb[0], rgb[1], rgb[2]);
                            }

                            m_scene.addEllipse(::round(x), ::round(y), 1, 1, QPen(pointColor), QBrush(pointColor, Qt::SolidPattern));
                            if (n > 0) {
                                m_scene.addLine(xyVector[0][0], xyVector[0][1], x, y, QPen(Qt::yellow));
                            }
                        }
                    }
                } else {
                    itk::Function::HSVColormapFunction<double, itk::RGBAPixel<unsigned char> >::Pointer colorFunc = itk::Function::HSVColormapFunction<double, itk::RGBAPixel<unsigned char> >::New();
                    colorFunc->SetMinimumInputValue(0);
                    colorFunc->SetMaximumInputValue(nPoints);

                    for (int i = 0; i < nPoints; i++) {
                        double x = particles->get(image * nParams + nDims * i);
                        double y = particles->get(image * nParams + nDims * i + 1);

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

                        m_scene.addEllipse(x-.5, y-.5, 1, 1, QPen(pointColor), QBrush(pointColor, Qt::SolidPattern));
                    }
                }
            }
        }

        if (ui.actionShowDisplacement->isChecked()) {
            if (!g_imageParticlesAlgo->GetDisplacementField().IsNull()) {
                SliceType::Pointer sliceImage = m_ImageList[image]->GetSlice();
                SliceIteratorType iter(sliceImage, sliceImage->GetBufferedRegion());
                DisplacementFieldType::Pointer field = g_imageParticlesAlgo->GetDisplacementField();
                for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
                    SliceType::IndexType idx = iter.GetIndex();
                    VectorType offset = field->GetPixel(idx);
                    VectorType imageOffset;
                    m_ImageList[image]->TransformToImageVector(offset, imageOffset);
                    VectorType src;
                    src[0] = idx[0];
                    src[1] = idx[1];
                    VectorType dst = src + imageOffset;
                    if (imageOffset[0] == 0 && imageOffset[1] == 0) {
                        continue;
                    }
                    cout << "Destination: " << dst << endl;
                    m_scene.addLine(src[0], src[1], dst[0], dst[1], QPen(Qt::yellow));
                }
            }
        }


        // bspline control points
        if (ui.actionShowBSplineControlPoints->isChecked()) {
            DisplacementFieldType::Pointer controlPoints = g_imageParticlesAlgo->GetBSplineRegistration()->GetControlPoints();
            if (controlPoints.IsNotNull()) {
                FieldIteratorType cpsIter(controlPoints, controlPoints->GetBufferedRegion());
                for (cpsIter.GoToBegin(); !cpsIter.IsAtEnd(); ++cpsIter) {
                    DisplacementFieldType::IndexType aIndex = cpsIter.GetIndex();
                    ContinuousIndexType aContIndex;
                    DisplacementFieldType::PointType aPoint;
                    VectorType aOffset = cpsIter.Get();
                    aContIndex[0] = aIndex[0] + aOffset[0];
                    aContIndex[1] = aIndex[1] + aOffset[1];
                    controlPoints->TransformContinuousIndexToPhysicalPoint(aContIndex, aPoint);
                    m_scene.addEllipse(aPoint[0] - .5, aPoint[1] - .5, 1, 1, QPen(Qt::white), QBrush(Qt::white, Qt::SolidPattern));
                }
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
        g_AnimFrame += m_Props.GetInt("animationInterleave", 10);
    } else if (ui.tabWidget->currentWidget() == ui.modelTab) {
        int imageIdx = ui.grayImages->currentIndex();
        if (imageIdx < 0) {
            m_Timer.stop();
            return;
        }

        vtkPolyData* poly = m_PropScene.FindPolyData("plane");
        if (poly == NULL || g_imageParticlesAlgo.IsNull()) {
            m_Timer.stop();
            return;
        }

//        if (g_AnimFrame < (int) g_Params.size()) {
//            OptimizerParametersType param = g_Params[g_AnimFrame];
//            if (param.GetSize() != 3 * poly->GetNumberOfPoints()) {
//                m_Timer.stop();
//                return;
//            }
        if (g_AnimFrame < (int) g_imageParticlesAlgo->GetNumberOfTraces()) {
            const VNLVector *param = g_imageParticlesAlgo->GetTraceParameters(g_AnimFrame);
            int nPoints = poly->GetNumberOfPoints();
            int nOffset = imageIdx * (nPoints * 2);
            for (int i = 0; i < nPoints; i++) {
                poly->GetPoints()->SetPoint(i, param->operator[](nOffset+2*i), param->operator[](nOffset+2*i+1), 0);
            }
            m_PropScene.ModifyLastActor();
            ui.statusbar->showMessage(QString("%1 frame showing").arg(g_AnimFrame));
            g_AnimFrame += m_Props.GetInt("animationInterleave", 1);
            m_Renderer->ResetCamera();
            m_Interactor->Render();
        }
    }
}

void MainWindow::on_actionRandomParticlesInit_triggered() {
    ui.toolBox->setCurrentWidget(ui.optimizerSettings);
    g_imageParticlesAlgo = ImageParticlesAlgorithm::New();
    g_imageParticlesAlgo->SetPropertyAccess(m_Props);
    g_imageParticlesAlgo->SetCurrentSliceAndView(GetCurrentView(), ui.sliceIndex->value());
    g_imageParticlesAlgo->SetImageList(&m_ImageList);
    g_imageParticlesAlgo->SetEventCallback(this);
    g_imageParticlesAlgo->CreateRandomInitialPoints(m_Props.GetInt("numberOfPoints", 100));
    g_constraint.SetImageList(&m_ImageList);
    updateScene();
}


void MainWindow::on_actionRunImageParticles_triggered() {
    // run optimization via itk-optimizers or ODE system
    ui.toolBox->setCurrentWidget(ui.optimizerSettings);
    ui.costPlot->graph()->clearData();

    // prepare image particles algorithm
    if (IsImageAvailable(0)) {
        if (g_imageParticlesAlgo.IsNotNull()) {
            ImageContainer::ClearDerivedViews();
            g_imageParticlesAlgo->SetImageList(&m_ImageList);
            g_imageParticlesAlgo->SetPropertyAccess(m_Props);

            // for viewing purpose
            g_constraint.SetImageList(&m_ImageList);
        }
    }

    if (g_imageParticlesAlgo.IsNotNull()) {
        if (ui.actionUseODESolver->isChecked()) {
            g_imageParticlesAlgo->RunODE();

        } else {
            g_imageParticlesAlgo->RunOptimization();
        }
    }
    updateScene();
}

void MainWindow::on_actionContinue_triggered() {
    ui.toolBox->setCurrentWidget(ui.optimizerSettings);
    if (g_imageParticlesAlgo.IsNotNull()) {
        if (ui.actionUseODESolver->isChecked()) {
            g_imageParticlesAlgo->RunODE();
            updateScene();
        } else {
            g_imageParticlesAlgo->SetPropertyAccess(m_Props);
            g_imageParticlesAlgo->ContinueOptimization();
            updateScene();
        }
    }
}


// apply TPS transform using current particle positions
void MainWindow::on_actionTPS_triggered() {
    // must after particles run
    if (g_imageParticlesAlgo.IsNull()) {
        return;
    }
    g_imageParticlesAlgo->ApplyTPSorEBSTransform(0);
}

// apply TPS transform using current particle positions
void MainWindow::on_actionEBS_triggered() {
    // must after particles run
    if (g_imageParticlesAlgo.IsNull()) {
        return;
    }
    g_imageParticlesAlgo->ApplyTPSorEBSTransform(1);
}

void MainWindow::on_actionR2LogRTPS_triggered() {
    // must after particles run
    if (g_imageParticlesAlgo.IsNull()) {
        return;
    }
    g_imageParticlesAlgo->ApplyTPSorEBSTransform(2);
}

void MainWindow::on_actionBSpline_triggered() {
    if (g_imageParticlesAlgo.IsNull()) {
        return;
    }
    g_imageParticlesAlgo->ApplyBSplineTransform();
}

void MainWindow::on_graphicsView_mousePressed(QMouseEvent* event) {
    QPoint o = event->pos();
    QPointF p = ui.graphicsView->mapToScene(o);

    cout << p.x() << "," << p.y() << endl;

    if (ui.derivedImages->currentIndex() > 0) {
        // g_imageParticlesAlgo->ProbeDerivedImage(ui.derivedImages->currentText().toUtf8().data());

    }

    if (GetCurrentImage() > -1) {
        my::ImplicitSurfaceConstraint::ContinuousIndexType idx;
        idx[0] = p.x();
        idx[1] = p.y();
        SliceType::IndexType idx2;
        idx2[0] = p.x(); idx2[1] = p.y();

        if (!g_constraint.IsInsideRegion(GetCurrentImage(), idx)) {
            return;
        }

        cout << "Distance: " << g_constraint.GetDistance(GetCurrentImage(), idx) << endl;
        cout << "Inside Offset: " << g_constraint.GetInsideOffset(GetCurrentImage(), idx2) << endl;
        cout << "Outside Offset: " << g_constraint.GetOutsideOffset(GetCurrentImage(), idx2) << endl;

        if (g_constraint.GetDistance(GetCurrentImage(), idx) < 0) {
            my::ImplicitSurfaceConstraint::DistanceVectorType offset = g_constraint.GetInsideOffset(GetCurrentImage(), idx2);
            m_scene.addLine(idx[0], idx[1], idx[0] + offset[0], idx[1] + offset[1], QPen(Qt::white));
        } else {
            my::ImplicitSurfaceConstraint::DistanceVectorType offset = g_constraint.GetOutsideOffset(GetCurrentImage(), idx2);
            m_scene.addLine(idx[0], idx[1], idx[0] + offset[0], idx[1] + offset[1], QPen(Qt::white));
        }

        my::ImplicitSurfaceConstraint::GradientPixelType gx = g_constraint.GetGradient(GetCurrentImage(), idx);
        cout << "Gradient: " << gx << endl;

        if (g_imageParticlesAlgo.IsNotNull()) {
            g_imageParticlesAlgo->OnClick(p.x(), p.y(), GetCurrentImage());
        }
    }
    
//    int xyRes = 10; //::round(::sqrt(m_Props.GetInt("numberOfPoints", 100)));
//    vtkPolyData* plane = m_PropScene.CreatePlane(xyRes, xyRes, p.x(), p.y());
//    plane->Print(cout);
//
//    m_PropScene.AddPolyData("plane", plane);
//    m_PropScene.SetColor(0, 0, 1);
//    m_PropScene.SetRepresentation(1);
//    m_Renderer->ResetCamera();
//    m_Interactor->Render();
//
//    cout << "Render a plane ..." << endl;

//    g_imageParticlesAlgo = ImageParticlesAlgorithm::New();
//    g_imageParticlesAlgo->SetPropertyAccess(m_Props);
//    g_imageParticlesAlgo->SetViewingDimension(GetCurrentView());
//    g_imageParticlesAlgo->SetImageList(&m_ImageList);
//    g_imageParticlesAlgo->SetEventCallback(this);
//    // g_imageParticlesAlgo->CreateRandomInitialPoints(m_Props.GetInt("numberOfPoints", 100));
//    g_imageParticlesAlgo->CreateInitialPoints(plane->GetPoints());
//    g_imageParticlesAlgo->SetSliceMarker(ui.sliceIndex->value());
    // updateScene();
}
