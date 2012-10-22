#include "mainwindow.h"
#include "QGraphicsView"
#include "QGraphicsScene"
#include "QFileDialog"
#include "QImage"
#include "QDebug"
#include "QGraphicsPixmapItem"
#include "vector"
#include "itkImageIO.h"
#include "itkScalarToARGBColormapImageFilter.h"
#include "itkMyFRPROptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkSphereBoundedGradientDescentOptimizer.h"
#include "itkARGBColormapFunction.h"
#include "armadillo"
#include "SurfaceEntropyCostFunction.h"

typedef itk::RGBAPixel<unsigned char> RGBAPixel;
typedef itk::RGBPixel<unsigned char> RGBPixel;
typedef itk::Image<double, 3> ImageType;
typedef itk::Image<RGBAPixel, 3> BitmapType;
typedef itk::ScalarToARGBColormapImageFilter<ImageType, BitmapType> ScalarToRGBFilter;
typedef itk::SphereBoundedGradientDescentOptimizer OptimizerType;
typedef arma::mat MatrixType;

std::vector<BitmapType::Pointer> bitmapList;
MatrixType pointList;

class OptimizerProgress: public itk::Command {
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
		cout << realCaller->GetCurrentIteration() << "; Cost = "
				<< realCaller->GetValue() << "; Parameters = "
				<< realCaller->GetCurrentPosition() << endl;
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
		ScalarToRGBFilter::Pointer rgbFilter = ScalarToRGBFilter::New();
		rgbFilter->SetInput(img);
		rgbFilter->SetNumberOfThreads(1);
		rgbFilter->Update();
		rgbImage = rgbFilter->GetOutput();
		bitmapList.push_back(rgbImage);
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
	ui.graphicsView->scale(4, 4);

	addImage(tr("/base/imageParticles/data/12.T2.slice.nrrd"));
	addImage(tr("/base/imageParticles/data/13.T2.slice.nrrd"));
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

	int currentImage = ui.listWidget->currentRow();
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
		for (int k = 0; k < bitmapList.size(); k++) {
			for (int j = 0; j < nY; j++) {
				for (int i = 0; i < nX; i++) {
					pointList(k, 2 * (j * nX + i)) = xCoords[i];
					pointList(k, 2 * (j * nX + i) + 1) = yCoords[j];
				}
			}
		}
	} else {
		pointList.randn(bitmapList.size(), 300 * 2);
		BitmapType::IndexType imageCenter;
		imageCenter[0] = bmpSz[0] / 2;
		imageCenter[1] = bmpSz[1] / 2;
		for (int i = 0; i < pointList.n_cols; i += 2) {
			for (int k = 0; k < bitmapList.size(); k++) {
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

void MainWindow::updateScene() {
	int currentImage = ui.listWidget->currentRow();
	gs.clear();
	BitmapType::Pointer bitmap = bitmapList[currentImage];
	QGraphicsPixmapItem* pixItem = gs.addPixmap(getPixmap(bitmap));
	pixItem->pos();

	if (pointList.n_cols > 0) {
		int nPoints = pointList.n_cols / 2;
		itk::Function::HSVColormapFunction<unsigned char, itk::RGBAPixel<unsigned char> >::Pointer colorFunc = itk::Function::HSVColormapFunction<unsigned char, itk::RGBAPixel<unsigned char> >::New();
		colorFunc->SetMinimumInputValue(0);
		colorFunc->SetMaximumInputValue(nPoints);
		for (int j = 0; j < nPoints; j++) {
			double x = pointList(currentImage, 2 * j);
			double y = pointList(currentImage, 2 * j + 1);
			itk::RGBAPixel<unsigned char> rgb = colorFunc->operator()(j / 2);
			QColor pointColor(rgb[0], rgb[1], rgb[2]);
			gs.addEllipse(x, y, 1, 1, QPen(pointColor),
					QBrush(pointColor, Qt::SolidPattern));
		}
	}

	currentImage = ui.listWidget->currentRow();
	BitmapType::SizeType bmpSz = bitmapList[0]->GetBufferedRegion().GetSize();

	gs.addEllipse(::round(bmpSz[0]/2)- 40,::round(bmpSz[1]/2)-40, 80, 80, QPen(Qt::white));
}

void MainWindow::on_actionRun_triggered() {
	if (bitmapList.size() == 0) {
		return;
	}

	BitmapType::SizeType imgSz = bitmapList[0]->GetBufferedRegion().GetSize();

	int currentImage = ui.listWidget->currentRow();

	const static int nDim = 2;
	typedef SurfaceEntropyCostFunction<nDim> CostFunction;
	CostFunction::PointContainerType initialPoints;
	initialPoints.resize(pointList.n_cols / nDim);
	for (int i = 0; i < initialPoints.size(); i++) {
		CostFunction::PointType pi;
		for (int j = 0; j < nDim; j++) {
			pi[j] = pointList(currentImage, i * nDim + j);
		}
		initialPoints[i] = pi;
	}

	CostFunction::Pointer costFunc = CostFunction::New();
	costFunc->Initialize(initialPoints);

	CostFunction::ParametersType initialParams =
			costFunc->GetInitialParameters();
	double cost = costFunc->GetValue(initialParams);
	cout << "Cost before optimization: " << cost << endl;

	cout  << "Intial Parameters: " << initialParams << endl;
	OptimizerProgress::Pointer progress = OptimizerProgress::New();

	OptimizerType::Pointer opti = OptimizerType::New();
	//opti->AddObserver(itk::IterationEvent(), progress);
	opti->SetCostFunction(costFunc);
	opti->SetInitialPosition(initialParams);
	opti->SetRadius(40);
	OptimizerType::VectorType center;
	center[0] = imgSz[0] / 2;
	center[1] = imgSz[1] / 2;
	opti->SetCenter(center);
	//opti->SetUseUnitLengthGradient(true);
	//opti->SetMaximumIteration(3);
	//opti->SetMaximumLineIteration(10);
	opti->SetNumberOfIterations(1000);
	opti->StartOptimization();

	CostFunction::ParametersType result = opti->GetCurrentPosition();
	for (int i = 0; i < result.GetSize(); i++) {
		pointList(currentImage, i) = result[i];
	}

	cout << result << endl;
	updateScene();

}

void MainWindow::resizeEvent(QResizeEvent* event) {
}

void MainWindow::addImage(QString filename) {
	loadImage(filename);
	ui.listWidget->addItem(filename);
	ui.listWidget->setCurrentRow(ui.listWidget->count() - 1);
}
