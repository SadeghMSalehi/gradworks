#include "mainwindow.h"
#include "QGraphicsView"
#include "QGraphicsScene"
#include "QFileDialog"
#include "QImage"
#include "QDebug"
#include "QGraphicsPixmapItem"
#include "QMessageBox"
#include "vector"
#include "algorithm"

//#include "SurfaceEntropyNLP.h"
#include "imageParticleTypes.h"

std::vector<BitmapType::Pointer> g_bitmapList;
std::vector<BitmapType::Pointer> g_gradmagBitmapList;
std::vector<BitmapType::Pointer> g_distanceMapBitmapList;
std::vector<BitmapType::Pointer> g_maskBitmapList;
std::vector<BitmapType::Pointer> g_cannyMaskBitmapList;

std::vector<ImageType::Pointer> g_imageList;
std::vector<GradientImageType::Pointer> g_gradientImageList;
std::vector<ImageType::Pointer> g_gradmagImageList;
std::vector<ImageType::Pointer> g_distanceMapList;
std::vector<DistanceVectorImageType::Pointer> g_distanceVectorList;
std::vector<ImageType::Pointer> g_boundaryMapList;
std::vector<ImageType::Pointer> g_cannyMaskList;

ListOfPointVectorType g_phantomParticles;
ListOfPointVectorType g_pointList;
ListOfParametersType g_pointHistory;

typedef SurfaceEntropyCostFunction<POINT_DIMENSIONS> CostFunctionType;

unsigned int g_pointHistoryIdx;
unsigned int g_numberOfPoints;
unsigned int g_numberOfIterations;

std::vector<double> g_costHistory;


class OptimizerProgress: public itk::Command {
private:
    ListOfParametersType* m_ParametersHistory;
    int m_Counter;

public:
	/** Standard class typedefs. */
	typedef OptimizerProgress Self;
	typedef itk::Command Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;

	/** Run-time type information (and related methods). */
	itkTypeMacro(OptimizerProgress, itk::Command)
	;

	itkNewMacro(OptimizerProgress)
	;


    void SetParticlesHistory(ListOfParametersType* parametersHistory) {
        if (parametersHistory == NULL) {
            return;
        }
        m_ParametersHistory = parametersHistory;
        m_Counter = 0;
    }

	/** Abstract method that defines the action to be taken by the command. */
	virtual void Execute(itk::Object *caller, const itk::EventObject & event) {
		this->Execute((const itk::Object*) caller, event);
	}

	/** Abstract method that defines the action to be taken by the command.
	 * This variant is expected to be used when requests comes from a
	 * const Object */
	virtual void Execute(const itk::Object *caller,
                         const itk::EventObject & event) {
		const OptimizerType* realCaller = dynamic_cast<const OptimizerType*>(caller);
		if (realCaller == NULL) {
			return;
		}
        if (++m_Counter % 1 == 0) {
            if (dynamic_cast<const LBFGSOptimizerType*>(caller) != NULL) {
                const LBFGSOptimizerType* opti = dynamic_cast<const LBFGSOptimizerType*>(caller);
                g_costHistory.push_back(opti->GetCachedValue());
            } else if (dynamic_cast<const FRPROptimizerType*>(caller) != NULL) {
                const FRPROptimizerType* opti = dynamic_cast<const FRPROptimizerType*>(caller);
                g_costHistory.push_back(opti->GetValue());
            } else if (dynamic_cast<const GDOptimizerType*>(caller) != NULL) {
                const GDOptimizerType* opti = dynamic_cast<const GDOptimizerType*>(caller);
                g_costHistory.push_back(opti->GetValue());
            }
            if (realCaller->GetCurrentPosition().GetSize() > 0) {
                m_ParametersHistory->push_back(realCaller->GetCurrentPosition());
            }

        }
	}

protected:
	OptimizerProgress() {
	}
	virtual ~OptimizerProgress() {
	}
private:
	OptimizerProgress(const Self &);        //purposely not implemented
	void operator=(const Self &); //purposely not implemented
};


static void createPhantomParticles(ImageType::Pointer edgeImg) {
    PointVectorType phantoms;
    itk::ImageRegionConstIteratorWithIndex<ImageType> iter(edgeImg, edgeImg->GetBufferedRegion());
    for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
        if (iter.Get() > 0) {
            phantoms.push_back(iter.GetIndex()[0]);
            phantoms.push_back(iter.GetIndex()[1]);
        }
    }
    g_phantomParticles.push_back(phantoms);
}

BitmapType::Pointer loadImage(QString f) {
	BitmapType::Pointer rgbImage;
	try {
		itkcmds::itkImageIO<ImageType> io;
		ImageType::Pointer img = io.ReadImageT(f.toAscii().data());
        g_imageList.push_back(img);

        ScalarToRGBFilter::Pointer rgbFilter = ScalarToRGBFilter::New();
		rgbFilter->SetInput(img);
		rgbFilter->SetNumberOfThreads(8);
        rgbFilter->SetAlphaValue(128);
		rgbFilter->Update();
		rgbImage = rgbFilter->GetOutput();
		g_bitmapList.push_back(rgbImage);

        GradientImageFilter::Pointer gradFilter = GradientImageFilter::New();
        gradFilter->SetInput(img);
        gradFilter->SetSigma(3);
        gradFilter->Update();
        g_gradientImageList.push_back(gradFilter->GetOutput());

        VectorMagnitudeImageFilter::Pointer magFilter = VectorMagnitudeImageFilter::New();
        magFilter->SetInput(g_gradientImageList.back());
        magFilter->Update();
        g_gradmagImageList.push_back(magFilter->GetOutput());

		ScalarToRGBFilter::Pointer rgbMagFilter = ScalarToRGBFilter::New();
		rgbMagFilter->SetInput(magFilter->GetOutput());
		rgbMagFilter->SetNumberOfThreads(1);
		rgbMagFilter->Update();
        BitmapType::Pointer rgbGradMagImage = rgbMagFilter->GetOutput();
        g_gradmagBitmapList.push_back(rgbGradMagImage);

	} catch (itk::ExceptionObject& ex) {
		cout << ex.what() << endl;
	}
	return rgbImage;
}


void loadMask(QString f) {
    itkcmds::itkImageIO<ImageType> io;
    ImageType::Pointer img = io.ReadImageT(f.toAscii().data());
    g_boundaryMapList.push_back(img);

    ScalarToRGBFilter::Pointer rgbFilter1 = ScalarToRGBFilter::New();
    rgbFilter1->SetInput(img);
    rgbFilter1->SetAlphaValue(128);
    rgbFilter1->Update();
    BitmapType::Pointer maskBitmap = rgbFilter1->GetOutput();
    g_maskBitmapList.push_back(maskBitmap);

    DistanceMapFilter::Pointer distmapFilter = DistanceMapFilter::New();
    distmapFilter->SetInput(img);
    distmapFilter->Update();
    ImageType::Pointer distImg = distmapFilter->GetOutput();
    g_distanceMapList.push_back(distmapFilter->GetOutput());
    
    DistanceVectorImageType::Pointer distVector = distmapFilter->GetVectorDistanceMap();
    g_distanceVectorList.push_back(distVector);

    ScalarToRGBFilter::Pointer rgbFilter = ScalarToRGBFilter::New();
    rgbFilter->SetInput(distImg);
    rgbFilter->SetNumberOfThreads(8);
    rgbFilter->Update();
    BitmapType::Pointer rgbImage = rgbFilter->GetOutput();
    g_distanceMapBitmapList.push_back(rgbImage);

    EdgeDetectionFilterType::Pointer edgeFilter = EdgeDetectionFilterType::New();
    edgeFilter->SetInput(img);
    edgeFilter->Update();
    ImageType::Pointer edgeImg = edgeFilter->GetOutput();

    ScalarToRGBFilter::Pointer edgeRGBFilter = ScalarToRGBFilter::New();
    edgeRGBFilter->SetInput(edgeImg);
    edgeRGBFilter->SetNumberOfThreads(8);
    edgeRGBFilter->Update();
    BitmapType::Pointer edgeRGB = edgeRGBFilter->GetOutput();

    g_cannyMaskList.push_back(edgeImg);
    g_cannyMaskBitmapList.push_back(edgeRGB);

    createPhantomParticles(edgeImg);
}

QPixmap getPixmap(BitmapType::Pointer bitmap) {
	BitmapType::SizeType bitmapSz = bitmap->GetBufferedRegion().GetSize();
	QImage qImg = QImage((unsigned char*) bitmap->GetBufferPointer(),
                         bitmapSz[0], bitmapSz[1], QImage::Format_ARGB32);
	return QPixmap::fromImage(qImg);
}

void drawGrid(QGraphicsScene* gs, double x0, double y0, double xSpacing,
              double ySpacing) {
	BitmapType::SizeType bmpSz = g_bitmapList[0]->GetBufferedRegion().GetSize();
	for (double x = x0; x < bmpSz[0]; x += xSpacing) {
		gs->addLine(x, 0, x, bmpSz[1], QPen(QColor(255, 255, 0)));
		cout << x << endl;
	}
	for (double y = y0; y < bmpSz[1]; y += ySpacing) {
		gs->addLine(0, y, bmpSz[0], y, QPen(QColor(255, 255, 0)));
		cout << y << endl;
	}
}


MainWindow::MainWindow(QWidget *parent) :
QMainWindow(parent) {
	ui.setupUi(this);

    m_particleColors = new QActionGroup(this);
    m_particleColors->addAction(ui.actionParticleBlack);
    m_particleColors->addAction(ui.actionParticleRed);
    m_particleColors->addAction(ui.actionParticleGreen);
    m_particleColors->addAction(ui.actionParticleBlue);
    m_particleColors->addAction(ui.actionParticleWhite);
    m_particleColors->addAction(ui.actionParticleHSV);
    ui.actionParticleHSV->setChecked(true);

    //    QSize gvSize = ui.graphicsView->size();
    //    ui.graphicsView->setSceneRect(0, 0, gvSize.width(), gvSize.height());
    //    qDebug() << "Scene Rect: " << ui.graphicsView->sceneRect();
	ui.graphicsView->setScene(&gs);
	ui.graphicsView->setRenderHints(
                                    QPainter::Antialiasing | QPainter::SmoothPixmapTransform);
	ui.graphicsView->scale(2, 2);
    ui.graphicsView->setBackgroundBrush(QBrush(Qt::black, Qt::SolidPattern));
    //	addImage(tr("/base/imageParticles/data/12.T2.slice.nrrd"));
    //	addImage(tr("/base/imageParticles/data/13.T2.slice.nrrd"));

#ifdef __APPLE__
    addImage(tr("/data/2D/circle.jpg"));
    loadMask(tr("/data/2D/circle-boundary.png"));
    addImage(tr("/data/2D/circle.jpg"));
    loadMask(tr("/data/2D/random-boundary.jpg"));
#else
    addImage(tr("/base/imageParticles/data/circle.jpg"));
    addImage(tr("/base/imageParticles/data/circle.jpg"));
    loadMask(tr("/base/imageParticles/data/circle-boundary.png"));
    loadMask(tr("/base/imageParticles/data/random-boundary.jpg"));
#endif

    QObject::connect(&m_Timer, SIGNAL(timeout()), this, SLOT(particleAnimationTimeout()));
    QObject::connect(ui.actionShowImage, SIGNAL(triggered()), this, SLOT(updateScene()));
    QObject::connect(ui.actionShowShapeMask, SIGNAL(triggered()), this, SLOT(updateScene()));
    QObject::connect(ui.actionShowShapeDistanceMap, SIGNAL(triggered()), this, SLOT(updateScene()));
    QObject::connect(m_particleColors, SIGNAL(triggered(QAction*)), this, SLOT(updateScene()));

    ui.listWidget->setCurrentRow(0);
}

MainWindow::~MainWindow() {
}

void MainWindow::on_actionOpen_triggered() {
	QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                    "/data", tr("Images (*.nrrd)"));
	if (fileName != NULL) {
		addImage(fileName);
	}
}

void MainWindow::on_actionClose_triggered() {
	this->close();
}

void MainWindow::on_actionZoomIn_triggered() {
	ui.graphicsView->scale(2, 2);
}

void MainWindow::on_actionZoomOut_triggered() {
	ui.graphicsView->scale(.5, .5);
}

/**
 * Deploy particles over the volume
 */
void MainWindow::on_actionDeploy_triggered() {
    g_pointList.clear();
    g_numberOfPoints = ui.numberOfParticles->value();
    const unsigned nVars = g_numberOfPoints * POINT_DIMENSIONS;
    for (unsigned j = 0; j < g_boundaryMapList.size(); j++) {
        std::vector<ImageType::IndexType> insideIndex;
        itk::ImageRegionConstIteratorWithIndex<ImageType> iter(g_boundaryMapList[j], g_boundaryMapList[j]->GetBufferedRegion());
        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            if (iter.Get() > 0) {
                insideIndex.push_back(iter.GetIndex());
            }
        }
        std::random_shuffle(insideIndex.begin(), insideIndex.end());
        if (insideIndex.size() > g_numberOfPoints) {
            PointVectorType initialPoints;
            initialPoints.reserve(nVars);
            for (unsigned i = 0; i < g_numberOfPoints; i++) {
                for (unsigned j = 0; j < POINT_DIMENSIONS; j++) {
                    int k = i;
                    initialPoints.push_back((float) insideIndex[k][j]);
                }
            }
            g_pointList.push_back(initialPoints);
        }
    }

	updateScene();
}

void MainWindow::on_listWidget_itemClicked(QListWidgetItem *item) {
}

void MainWindow::on_listWidget_currentRowChanged(int currentRow) {
	updateScene();
}

/**
 * set up cost function to compute optimal particle positions
 * parameters must be set up due to pre-processing
 */
void MainWindow::on_actionRun_triggered() {
    g_numberOfIterations = ui.numberOfIterations->value();
    g_costHistory.clear();
    g_pointHistory.clear();
    g_pointHistoryIdx = 0;

    m_CostFunc = CostFunctionType::New();
    m_CostFunc->SetSigma(ui.sigma->value());
    m_CostFunc->SetCutoffDistance(ui.cutoffDistance->value());
    m_CostFunc->SetPhantomCutoffDistance(ui.cutoffDistancePhantom->value());
    m_CostFunc->SetMaxKappa(ui.maxKappa->value());
    m_CostFunc->SetUseAdaptiveSampling(ui.actionAdaptiveSampling->isChecked());
    for (unsigned i = 0; i < g_imageList.size() && i < g_pointList.size(); i++) {
        m_CostFunc->AddSubjects(g_pointList[i], g_gradmagImageList[i], g_boundaryMapList[i]);
    }

    runOptimization();
}

void MainWindow::on_actionContinueOptimization_triggered() {
	g_numberOfIterations = ui.numberOfIterations->value();
    if (m_CostFunc.IsNull()) {
        QMessageBox::information(this, tr("Error"), tr("Cost function is not set up!"));
        return;
    }
    runOptimization();
}


OptimizerParameters optimizeByGD(CostFunctionType::Pointer costFunc, OptimizerParameters initialParams, OptimizerProgress::Pointer progress) {
	GDOptimizerType::Pointer opti = GDOptimizerType::New();
	opti->AddObserver(itk::IterationEvent(), progress.GetPointer());
	opti->SetCostFunction(costFunc);
	opti->SetInitialPosition(initialParams);
    opti->SetNumberOfIterations(g_numberOfIterations);
    cout << "Starting optimization ..." << endl;
	opti->StartOptimization();
	cout << "Optimization done ..." << endl;
    costFunc->Print();
    return opti->GetCurrentPosition();
}

OptimizerParameters optimizeByLBFGS(CostFunctionType::Pointer costFunc, OptimizerParameters initialParams, OptimizerProgress::Pointer progress) {
	LBFGSOptimizerType::Pointer opti = LBFGSOptimizerType::New();
	opti->AddObserver(itk::IterationEvent(), progress.GetPointer());
	opti->SetCostFunction(costFunc);
	opti->SetInitialPosition(initialParams);
	opti->SetMinimize(true);

    cout << "Starting optimization ..." << endl;
	opti->StartOptimization();
	cout << "Optimization done ..." << endl;
    costFunc->Print();
    return opti->GetCurrentPosition();
}

OptimizerParameters optimizeByFRPR(CostFunctionType::Pointer costFunc, OptimizerParameters initialParams, OptimizerProgress::Pointer progress) {
	FRPROptimizerType::Pointer opti = FRPROptimizerType::New();
	opti->AddObserver(itk::IterationEvent(), progress.GetPointer());
	opti->SetCostFunction(costFunc);
	opti->SetInitialPosition(initialParams);
    opti->SetMaximumIteration(g_numberOfIterations);
    opti->SetMaximumLineIteration(10);
    opti->SetStepLength(.25);
    opti->SetUseUnitLengthGradient(true);
    cout << "Starting optimization ..." << endl;
	opti->StartOptimization();
	cout << "Optimization done ..." << endl;
    costFunc->Print();
    return opti->GetCurrentPosition();
}

void MainWindow::on_optiCG_toggled(bool toggled) {
}

void MainWindow::on_optiLBFGS_toggled(bool toggled) {

}

void MainWindow::on_optiGD_toggled(bool toggled) {

}



void MainWindow::runOptimization() {
	if (m_CostFunc.IsNull() || g_pointList.size() == 0) {
		return;
	}

    OptimizerParameters initialParams;
    ConvertListOfPointVectorsToParameters(g_pointList, initialParams);

	OptimizerProgress::Pointer progress = OptimizerProgress::New();
    progress->SetParticlesHistory(&g_pointHistory);
    
    OptimizerParameters result;
    if (ui.optiCG->isChecked()) {
        result = optimizeByFRPR(m_CostFunc, initialParams, progress);
    } else if (ui.optiGD->isChecked()) {
        result = optimizeByGD(m_CostFunc, initialParams, progress);
    } else if (ui.optiLBFGS->isChecked()) {
        result = optimizeByLBFGS(m_CostFunc, initialParams, progress);
    }

    // draw cost function
    ConvertParametersToListOfPointVectors(result, g_pointList.size(), g_pointList[0].size(), g_pointList);
    QVector<double> x(g_costHistory.size()), y(g_costHistory.size()); // initialize with entries 0..100
    for (int i = 0; i < x.count(); ++i)
    {
        x[i] = i;
        y[i] = g_costHistory[i];
    }

    QCustomPlot* customPlot = ui.customPlot;
    customPlot->clearGraphs();
    // create graph and assign data to it:
    customPlot->addGraph();
    customPlot->graph(0)->setData(x, y);
    customPlot->rescaleAxes();

    // give the axes some labels:
    customPlot->xAxis->setLabel("iteration");
    customPlot->yAxis->setLabel("cost");
    customPlot->replot();

	updateScene();
}

void MainWindow::on_actionPointSave_triggered() {
    SaveListOfPointVectors(g_pointList, "pointSamples.txt");
    updateScene();
}

void MainWindow::on_actionPointLoad_triggered() {
    LoadListOfPointVectors("pointSamples.txt", g_pointList);
    updateScene();
}

void MainWindow::on_actionPlayTrace_triggered() {
    g_pointHistoryIdx = 0;
    m_Timer.start(ui.animationInterval->value());
}


void MainWindow::on_animationInterval_valueChanged(int v) {
    if (m_Timer.isActive()) {
        m_Timer.stop();
        m_Timer.start(ui.animationInterval->value());
    }
}


void MainWindow::particleAnimationTimeout() {
    if (g_pointHistoryIdx >= g_pointHistory.size()) {
        ui.statusBar->showMessage(QString("Animation Done.."));
        m_Timer.stop();
    } else {
        playScene();
        ui.statusBar->showMessage(QString("Playing %1 frame ..").arg(g_pointHistoryIdx));
        g_pointHistoryIdx++;
    }
}

void MainWindow::on_actionLoadTrace_triggered() {
    /*
	QString fileName = QFileDialog::getOpenFileName(this, tr("Load Text File"), "/data", tr("Texts (*.txt)"));
	if (fileName != NULL) {
        g_pointHistory.load(fileName.toAscii().data());
	}
    */
}

void MainWindow::on_actionSaveTrace_triggered() {
    /*
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Text File"), "/data", tr("Texts (*.txt)"));
	if (fileName != NULL) {
        g_pointHistory.save(fileName.toAscii().data());
	}
    */
}

void MainWindow::on_actionShowPlotWindow_triggered() {
    ui.dockWidget_2->show();
}

void MainWindow::updateScene() {
	unsigned currentImage = ui.listWidget->currentRow();
	gs.clear();

    if (ui.actionShowImage->isChecked() && currentImage < g_gradmagBitmapList.size()) {
        BitmapType::Pointer bitmap2;
        bitmap2 = g_gradmagBitmapList[currentImage];
        gs.addPixmap(getPixmap(bitmap2));
    }
    if (ui.actionShowShapeMask->isChecked() && currentImage < g_maskBitmapList.size()) {
        if (!g_maskBitmapList.empty()) {
            BitmapType::Pointer bitmap = g_maskBitmapList[currentImage];
            gs.addPixmap(getPixmap(bitmap));
        }
    } else if (ui.actionShowShapeDistanceMap->isChecked() && currentImage < g_distanceMapBitmapList.size()) {
        if (!g_distanceMapBitmapList.empty()) {
            BitmapType::Pointer bitmap = g_distanceMapBitmapList[currentImage];
            gs.addPixmap(getPixmap(bitmap));
        }
    }


	if (g_pointList.size() > 0) {
		int nPoints = g_pointList[0].size() / POINT_DIMENSIONS;
		itk::Function::HSVColormapFunction<double, itk::RGBAPixel<unsigned char> >::Pointer colorFunc = itk::Function::HSVColormapFunction<double, itk::RGBAPixel<unsigned char> >::New();
		colorFunc->SetMinimumInputValue(0);
		colorFunc->SetMaximumInputValue(nPoints);

		for (int j = 0; j < nPoints; j++) {
			double x = g_pointList[currentImage][POINT_DIMENSIONS * j];
			double y = g_pointList[currentImage][POINT_DIMENSIONS * j + 1];

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
                itk::RGBAPixel<unsigned char> rgb = colorFunc->operator()(j);
                pointColor = QColor(rgb[0], rgb[1], rgb[2]);
            }

			gs.addEllipse(x, y, 1, 1, QPen(pointColor),
                          QBrush(pointColor, Qt::SolidPattern));
		}
	}

    /*
	currentImage = ui.listWidget->currentRow();

    bool drawPhantoms = false && (g_phantomParticles.size() > 0);
    if (drawPhantoms) {
        std::vector<float>& phantoms = g_phantomParticles[currentImage];
        for (unsigned j = 0; j < phantoms.size(); j += 2) {
            gs.addEllipse(phantoms[j], phantoms[j+1], 1, 1, QPen(Qt::white), QBrush(Qt::white, Qt::SolidPattern));
        }
    }
    */
}



void MainWindow::playScene() {
    const unsigned nSubj = g_pointList.size();
    const unsigned nAnimationFrames = g_pointHistory.size();
	if (nAnimationFrames > 0 && nSubj > 0 && g_pointHistoryIdx < nAnimationFrames) {
        const unsigned nVars = g_pointList[0].size();
		const unsigned nPoints = nVars / POINT_DIMENSIONS;
        const unsigned currentImage = ui.listWidget->currentRow();
        if (currentImage >= nSubj) {
            cout << "Can't play due to current image change" << endl;
            return;
        }
        if (g_pointHistory[g_pointHistoryIdx].size() != nSubj * nVars) {
            cout << "Can't play frames due to the mismatch of number of particles ..." << endl;
            return;
        }
        const unsigned nSubjectOffset = currentImage * nVars;

        gs.clear();

        BitmapType::Pointer bitmap = g_bitmapList[currentImage];
        gs.addPixmap(getPixmap(bitmap));


        // hope to change particle colors to represent kappa values
		itk::Function::JetColormapFunction<double, itk::RGBAPixel<unsigned char> >::Pointer colorFunc = itk::Function::JetColormapFunction<double, itk::RGBAPixel<unsigned char> >::New();
		colorFunc->SetMinimumInputValue(0);
		colorFunc->SetMaximumInputValue(255);
		for (unsigned j = 0; j < nPoints; j++) {
			double x = g_pointHistory[g_pointHistoryIdx][nSubjectOffset + POINT_DIMENSIONS*j];
			double y = g_pointHistory[g_pointHistoryIdx][nSubjectOffset + POINT_DIMENSIONS*j + 1];

            ImageType::IndexType xyIdx;
            xyIdx[0] = x;
            xyIdx[1] = y;
            ImageType::PixelType xyVal = g_imageList[currentImage]->GetPixel(xyIdx);
            itk::RGBAPixel<unsigned char> rgb = colorFunc->operator()(xyVal);
			QColor pointColor(rgb[0], rgb[1], rgb[2]);
			gs.addEllipse(x, y, 1, 1, QPen(pointColor),
                          QBrush(pointColor, Qt::SolidPattern));
		}
	}
}

void MainWindow::resizeEvent(QResizeEvent* event) {
}

void MainWindow::addImage(QString filename) {
	loadImage(filename);
	ui.listWidget->addItem(filename);
	// ui.listWidget->setCurrentRow(ui.listWidget->count() - 1);
}

void MainWindow::on_graphicsView_mousePressed(QMouseEvent* event) {
    int currentImage = ui.listWidget->currentRow();
    if (currentImage < 0) {
        return;
    }

    QPointF xy = ui.graphicsView->mapToScene(event->pos());
    ImageType::Pointer sourceImage = g_distanceMapList[currentImage];
    ImageType::IndexType sourceIdx;
    sourceIdx[0] = xy.x();
    sourceIdx[1] = xy.y();
    ImageType::PixelType pixel = sourceImage->GetPixel(sourceIdx);
    cout << sourceIdx[0] << "," << sourceIdx[1] << "," << pixel << endl;
}
