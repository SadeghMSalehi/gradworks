#include "mainwindow.h"
#include "QGraphicsView"
#include "QGraphicsScene"
#include "QFileDialog"
#include "QImage"
#include "QDebug"
#include "QGraphicsPixmapItem"
#include "vector"

#include "SurfaceEntropyCostFunction.h"
//#include "SurfaceEntropyNLP.h"
#include "imageParticleTypes.h"
#include "NBody.h"


#define RADIUS 60
#define NUMBER_OF_POINTS 500
#define POINT_DIMENSIONS 2

std::vector<BitmapType::Pointer> bitmapList;
std::vector<BitmapType::Pointer> gradMagBitmapList;
std::vector<ImageType::Pointer> imageList;
std::vector<GradientImageType::Pointer> gradientImageList;
std::vector<ImageType::Pointer> gradMagImageList;

unsigned int g_pointHistoryIdx;
unsigned int g_numberOfPoints;
unsigned int g_numberOfIterations;
MatrixType g_pointHistory;
typedef std::vector<OptimizerType::ParametersType> ParametersHistoryType;
std::vector<double> g_costHistory;

MatrixType pointList;

class OptimizerProgress: public itk::Command {
private:
    ParametersHistoryType* m_ParametersHistory;
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


    void SetParticlesHistory(ParametersHistoryType* parametersHistory) {
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
		const OptimizerType* realCaller =
        dynamic_cast<const OptimizerType*>(caller);
		if (realCaller == NULL) {
			return;
		}
        bool traceSteps = false;
        if (traceSteps) {
            cout << realCaller->GetCurrentIteration() << "; Cost = "
            << realCaller->GetValue() << "; Parameters = "
            << realCaller->GetCurrentPosition() << endl;
        }
        if (++m_Counter % 1 == 0) {
            m_ParametersHistory->push_back(realCaller->GetCurrentPosition());
            g_costHistory.push_back(realCaller->GetValue());
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

BitmapType::Pointer loadImage(QString f) {
	BitmapType::Pointer rgbImage;
	try {
		itkcmds::itkImageIO<ImageType> io;
		ImageType::Pointer img = io.ReadImageT(f.toAscii().data());
        imageList.push_back(img);

		ScalarToRGBFilter::Pointer rgbFilter = ScalarToRGBFilter::New();
		rgbFilter->SetInput(img);
		rgbFilter->SetNumberOfThreads(1);
		rgbFilter->Update();
		rgbImage = rgbFilter->GetOutput();
		bitmapList.push_back(rgbImage);

        GradientImageFilter::Pointer gradFilter = GradientImageFilter::New();
        gradFilter->SetInput(img);
        gradFilter->SetSigma(.5);
        gradFilter->Update();
        gradientImageList.push_back(gradFilter->GetOutput());

        VectorMagnitudeImageFilter::Pointer magFilter = VectorMagnitudeImageFilter::New();
        magFilter->SetInput(gradientImageList.back());
        magFilter->Update();
        gradMagImageList.push_back(magFilter->GetOutput());


		ScalarToRGBFilter::Pointer rgbMagFilter = ScalarToRGBFilter::New();
		rgbMagFilter->SetInput(magFilter->GetOutput());
		rgbMagFilter->SetNumberOfThreads(1);
		rgbMagFilter->Update();
        BitmapType::Pointer rgbGradMagImage = rgbMagFilter->GetOutput();
        gradMagBitmapList.push_back(rgbGradMagImage);

	} catch (itk::ExceptionObject& ex) {
		cout << ex.what() << endl;
	}
	return rgbImage;
}

QPixmap getPixmap(BitmapType::Pointer bitmap) {
	BitmapType::SizeType bitmapSz = bitmap->GetBufferedRegion().GetSize();
	QImage qImg = QImage((unsigned char*) bitmap->GetBufferPointer(),
                         bitmapSz[0], bitmapSz[1], QImage::Format_ARGB32);
	return QPixmap::fromImage(qImg);
}

void drawGrid(QGraphicsScene* gs, double x0, double y0, double xSpacing,
              double ySpacing) {
	BitmapType::SizeType bmpSz = bitmapList[0]->GetBufferedRegion().GetSize();
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

    //    QSize gvSize = ui.graphicsView->size();
    //    ui.graphicsView->setSceneRect(0, 0, gvSize.width(), gvSize.height());
    //    qDebug() << "Scene Rect: " << ui.graphicsView->sceneRect();
	ui.graphicsView->setScene(&gs);
	ui.graphicsView->setRenderHints(
                                    QPainter::Antialiasing | QPainter::SmoothPixmapTransform);
	ui.graphicsView->scale(2, 2);
    //	addImage(tr("/base/imageParticles/data/12.T2.slice.nrrd"));
    //	addImage(tr("/base/imageParticles/data/13.T2.slice.nrrd"));

#ifdef __APPLE__
    addImage(tr("/data/2D/circle.jpg"));
    addImage(tr("/data/2D/circle.jpg"));
#else
    addImage(tr("/base/imageParticles/data/circle.jpg"));
    addImage(tr("/base/imageParticles/data/circle.jpg"));
#endif
    QObject::connect(&m_Timer, SIGNAL(timeout()), this, SLOT(on_timer_timeout()));
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
 * assume image sizes are same
 *
 */
void MainWindow::on_actionDeploy_triggered() {
	if (bitmapList.size() == 0) {
		return;
	}

	BitmapType::SizeType bmpSz = bitmapList[0]->GetBufferedRegion().GetSize();
	bool uniformDistribution = false;
	/*
	 * uniform distribution
	 */
	if (uniformDistribution) {
		int x0 = 0, y0 = 0;
		int nX = 10, nY = 20;
		int nPoints = nX * nY;
		int distX = (int) round((bmpSz[0] - x0) / (nX - 1));
		int distY = (int) round((bmpSz[1] - y0) / (nY - 1));
		pointList.zeros(bitmapList.size(), nPoints * 2);
		std::vector<int> xCoords, yCoords;
		xCoords.resize(nX);
		yCoords.resize(nY);
		for (int i = 0; i < nX; i++) {
			xCoords[i] = x0 + i * distX;
		}
		for (int j = 0; j < nY; j++) {
			yCoords[j] = y0 + j * distY;
		}
		for (unsigned int k = 0; k < bitmapList.size(); k++) {
			for (int j = 0; j < nY; j++) {
				for (int i = 0; i < nX; i++) {
					pointList(k, 2 * (j * nX + i)) = xCoords[i];
					pointList(k, 2 * (j * nX + i) + 1) = yCoords[j];
				}
			}
		}
	} else {
        g_numberOfPoints = ui.numberOfParticles->value();
		pointList.randn(bitmapList.size(), g_numberOfPoints * POINT_DIMENSIONS);
		BitmapType::IndexType imageCenter;
		imageCenter[0] = bmpSz[0] / 2;
		imageCenter[1] = bmpSz[1] / 2;
		for (unsigned int i = 0; i < pointList.n_cols; i += POINT_DIMENSIONS) {
			for (unsigned int k = 0; k < bitmapList.size(); k++) {
				pointList(k, i) = 10 * pointList(k, i) + imageCenter[0];
				pointList(k, i + 1) = 10 * pointList(k, i + 1) + imageCenter[1];
				if (pointList(k, i) < 0) {
					pointList(k, i) = imageCenter[0];
				} else if (pointList(k, i) >= bmpSz[1]) {
					pointList(k, i) = imageCenter[1];
				}
				if (pointList(k, i + 1) < 0) {
					pointList(k, i + 1) = imageCenter[0];
				} else if (pointList(k, i + 1) >= bmpSz[1]) {
					pointList(k, i + 1) = imageCenter[1];
				}
			}
		}
	}

	updateScene();
}

void MainWindow::on_listWidget_itemClicked(QListWidgetItem *item) {
}

void MainWindow::on_listWidget_currentRowChanged(int currentRow) {
	updateScene();
}

void MainWindow::on_actionRun_triggered() {
	if (bitmapList.size() == 0) {
		return;
	}
	BitmapType::SizeType imgSz = bitmapList[0]->GetBufferedRegion().GetSize();

    bool useIpOpt = false;
    if (useIpOpt) {
        /*
         SurfaceEntropyNLP* nlp = new SurfaceEntropyNLP();
         nlp->SetInitialPoints(pointList);
         arma::vec imageRegion;
         imageRegion.zeros(2);
         imageRegion[0] = imgSz[0];
         imageRegion[1] = imgSz[1];
         nlp->SetRegion(imageRegion);
         nlp->SetCenter(imageRegion/2);
         nlp->SetRadius(arma::vec("40 40"));
         nlp->SetNumberOfPoints(pointList.n_cols / 2);

         Ipopt::SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
         app->Options()->SetNumericValue("tol", 1e-9);
         app->Options()->SetStringValue("mu_strategy", "adpative");
         app->Options()->SetStringValue("hessian_approximation", "limited-memory");

         Ipopt::ApplicationReturnStatus status;
         status = app->Initialize();
         if (status != Solve_Succeeded) {
         printf("Error during initializtion");
         }

         status = app->OptimizeTNLP(Ipopt::SmartPtr<TNLP>(nlp));
         if (status == Solve_Succeeded) {
         printf("Problem Solved");
         } else if (status == Invalid_Number_Detected) {
         cout << "Invalid Number Detected" << endl;
         } else {
         printf("Problem Failed: %d", status);
         }
         arma::vec points = nlp->GetResultPoints();
         for (int j = 0; j < points.n_elem; j++) {
         pointList.at(0, j) = points.at(j);
         }
         updateScene();        //delete nlp;
         return;
         */
    }



    g_numberOfIterations = ui.numberOfIterations->value();
    g_costHistory.clear();

    cout << "Number of Iterations:" << g_numberOfIterations << endl;
	int currentImage = ui.listWidget->currentRow();


	const static int nDim = POINT_DIMENSIONS;
	typedef SurfaceEntropyCostFunction<nDim> CostFunction;
	CostFunction::PointContainerType initialPoints;
	initialPoints.resize(pointList.n_cols / nDim);
	for (unsigned int i = 0; i < initialPoints.size(); i++) {
		CostFunction::PointType pi;
		for (int j = 0; j < nDim; j++) {
			pi[j] = pointList(currentImage, i * nDim + j);
		}
		initialPoints[i] = pi;
	}

	CostFunction::Pointer costFunc = CostFunction::New();
    costFunc->SetImage(gradMagImageList[currentImage]);
    costFunc->SetUseAdaptiveSampling(ui.actionAdaptiveSampling->isChecked());

	costFunc->Initialize(initialPoints);


	CostFunction::ParametersType initialParams =
    costFunc->GetInitialParameters();
	double cost = costFunc->GetValue(initialParams);
	cout << "Cost before optimization: " << cost << endl;

    ParametersHistoryType positionHistory;
	cout  << "Intial Parameters: " << initialParams << endl;
	OptimizerProgress::Pointer progress = OptimizerProgress::New();
    progress->SetParticlesHistory(&positionHistory);
    
	OptimizerType::Pointer opti = OptimizerType::New();
	opti->AddObserver(itk::IterationEvent(), progress.GetPointer());
	opti->SetCostFunction(costFunc);
	opti->SetInitialPosition(initialParams);
	opti->SetRadius(RADIUS);
    opti->SetMaximumStepLength(1);
    
	OptimizerType::VectorType center;
	center[0] = imgSz[0] / 2;
	center[1] = imgSz[1] / 2;
	opti->SetCenter(center);
    opti->SetNumberOfIterations(g_numberOfIterations);
    
	//opti->SetUseUnitLengthGradient(true);
	//opti->SetMaximumIteration(3);
	//opti->SetMaximumLineIteration(10);

    //opti->SetUseUnitLengthGradient(true);


	opti->StartOptimization();

	CostFunction::ParametersType result = opti->GetCurrentPosition();
	for (unsigned int i = 0; i < result.GetSize(); i++) {
		pointList(currentImage, i) = result[i];
	}

    cout << "Number of Iteration Performed: " << opti->GetCurrentIteration() << endl;
    cout << "Number of traces: " << positionHistory.size() << endl;
    g_pointHistory.zeros(positionHistory.size(), pointList.n_cols);
    for (unsigned int i = 0; i < g_pointHistory.n_rows; i++) {
        for (unsigned int j = 0; j < g_pointHistory.n_cols; j++) {
            g_pointHistory.at(i,j) = positionHistory[i][j];
        }
    }

    // generate some data:
    QVector<double> x(g_costHistory.size()), y(g_costHistory.size()); // initialize with entries 0..100
    for (int i = 0; i < x.count(); ++i)
    {
        x[i] = i;
        y[i] = g_costHistory[i];

        cout << y[i] << endl;
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

void MainWindow::on_actionNBodySimulation_triggered() {
	if (pointList.n_elem == 0 || bitmapList.size() < 1) {
		return;
	}
	g_numberOfIterations = ui.numberOfIterations->value();

	cout << "Running nbody simulation: " << g_numberOfIterations << endl;
	int currentImage = ui.listWidget->currentRow();
	int nPoints = pointList.n_cols / POINT_DIMENSIONS;
	BitmapType::SizeType imgSz = bitmapList[0]->GetBufferedRegion().GetSize();

	NBodySystem nBodySystem;
	nBodySystem.setBounds(0, 0, imgSz[0], imgSz[1]);
	for (int i = 0; i < pointList.n_cols; i += POINT_DIMENSIONS) {
		Body b;
		b.x = pointList.at(currentImage, i);
		b.y = pointList.at(currentImage, i+1);
		b.mass = 1;
		nBodySystem.addBody(b);
	}

	cout << "particles added" << endl;
	for (unsigned int i = 0; i < g_numberOfIterations; i++) {
		nBodySystem.advance(0.01);
	}

	cout << "iteration done" << endl;
	for (int i = 0; i < nPoints; i++) {
		pointList.at(currentImage, POINT_DIMENSIONS*i) = nBodySystem.bodies[i].x;
		pointList.at(currentImage, POINT_DIMENSIONS*i+1) = nBodySystem.bodies[i].y;
	}

    pointList.print();

	updateScene();
}

void MainWindow::on_actionPointSave_triggered() {
    pointList.save("pointSamples.txt");
    updateScene();
}

void MainWindow::on_actionPointLoad_triggered() {
    pointList.load("pointSamples.txt");
    updateScene();
}

void MainWindow::on_actionPlayTrace_triggered() {
    g_pointHistoryIdx = 0;
    m_Timer.start(ui.animationInterval->value());
}

void MainWindow::on_actionLoadTrace_triggered() {
	QString fileName = QFileDialog::getOpenFileName(this, tr("Load Text File"), "/data", tr("Texts (*.txt)"));
	if (fileName != NULL) {
        g_pointHistory.load(fileName.toAscii().data());
	}
}

void MainWindow::on_actionSaveTrace_triggered() {
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Text File"), "/data", tr("Texts (*.txt)"));
	if (fileName != NULL) {
        g_pointHistory.save(fileName.toAscii().data());
	}
}

void MainWindow::on_actionShowPlotWindow_triggered() {
    ui.dockWidget_2->show();
}

void MainWindow::on_timer_timeout() {
    if (g_pointHistoryIdx >= g_pointHistory.n_rows) {
        m_Timer.stop();
    } else {
        playScene();
        g_pointHistoryIdx++;
    }
}

void MainWindow::updateScene() {
	int currentImage = ui.listWidget->currentRow();
	gs.clear();
	BitmapType::Pointer bitmap = gradMagBitmapList[currentImage];
	QGraphicsPixmapItem* pixItem = gs.addPixmap(getPixmap(bitmap));
	pixItem->pos();

	if (pointList.n_cols > 0) {
		int nPoints = pointList.n_cols / POINT_DIMENSIONS;
        cout << "Number of points: " << nPoints << endl;
		itk::Function::HSVColormapFunction<double, itk::RGBAPixel<unsigned char> >::Pointer colorFunc = itk::Function::HSVColormapFunction<double, itk::RGBAPixel<unsigned char> >::New();
		colorFunc->SetMinimumInputValue(0);
		colorFunc->SetMaximumInputValue(nPoints);

		for (int j = 0; j < nPoints; j++) {
			double x = pointList(currentImage, POINT_DIMENSIONS * j);
			double y = pointList(currentImage, POINT_DIMENSIONS * j + 1);
			itk::RGBAPixel<unsigned char> rgb = colorFunc->operator()(j);
			QColor pointColor(rgb[0], rgb[1], rgb[2]);
			gs.addEllipse(x, y, 1, 1, QPen(pointColor),
                          QBrush(pointColor, Qt::SolidPattern));
		}
	}

	currentImage = ui.listWidget->currentRow();
    bool drawBoundary = false;
    if (drawBoundary) {
        BitmapType::SizeType bmpSz = bitmapList[0]->GetBufferedRegion().GetSize();
        gs.addEllipse(::round(bmpSz[0]/2)- RADIUS,::round(bmpSz[1]/2)-RADIUS, RADIUS*2, RADIUS*2, QPen(Qt::white));
    }
}

void MainWindow::playScene() {
	if (g_pointHistory.n_rows > 0 && g_pointHistoryIdx < g_pointHistory.n_rows) {
        gs.clear();
        
        int currentImage = ui.listWidget->currentRow();
        BitmapType::Pointer bitmap = gradMagBitmapList[currentImage];
        QGraphicsPixmapItem* pixItem = gs.addPixmap(getPixmap(bitmap));
        pixItem->pos();

		int nPoints = g_pointHistory.n_cols / POINT_DIMENSIONS;
		itk::Function::HSVColormapFunction<double, itk::RGBAPixel<unsigned char> >::Pointer colorFunc = itk::Function::HSVColormapFunction<double, itk::RGBAPixel<unsigned char> >::New();
		colorFunc->SetMinimumInputValue(0);
		colorFunc->SetMaximumInputValue(nPoints);
		for (int j = 0; j < nPoints; j++) {
			double x = g_pointHistory.at(g_pointHistoryIdx, POINT_DIMENSIONS*j);
			double y = g_pointHistory.at(g_pointHistoryIdx, POINT_DIMENSIONS*j + 1);
			itk::RGBAPixel<unsigned char> rgb = colorFunc->operator()(j);
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
	ui.listWidget->setCurrentRow(ui.listWidget->count() - 1);
}
