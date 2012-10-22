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
#include "armadillo"
#include "SurfaceEntropyCostFunction.h"

typedef itk::RGBAPixel<unsigned char> RGBAPixel;
typedef itk::RGBPixel<unsigned char> RGBPixel;
typedef itk::Image<double,3> ImageType;
typedef itk::Image<RGBAPixel,3> BitmapType;
typedef itk::ScalarToARGBColormapImageFilter<ImageType,BitmapType> ScalarToRGBFilter;
typedef arma::mat MatrixType;

std::vector<BitmapType::Pointer> bitmapList;
MatrixType pointList;

template<typename T>
class ImageParticlesAlgorithm {
public:
    typedef std::vector<typename T::Pointer> ImageList;
    typedef arma::mat MatrixType;

    ImageParticlesAlgorithm(ImageList imgList, MatrixType pointList) {
        m_ImageList = imgList;
        m_PointList = pointList;
    }

    virtual ~ImageParticlesAlgorithm() {
    }

    void Run() {
        
    }
    
private:
    ImageList m_ImageList;
    MatrixType m_PointList;

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
    QImage qImg = QImage((unsigned char*) bitmap->GetBufferPointer(), bitmapSz[0], bitmapSz[1], QImage::Format_ARGB32);
    return QPixmap::fromImage(qImg);
}

void drawGrid(QGraphicsScene* gs, double x0, double y0, double xSpacing, double ySpacing) {
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

MainWindow::MainWindow(QWidget *parent): QMainWindow(parent) {
    ui.setupUi(this);
//    QSize gvSize = ui.graphicsView->size();
//    ui.graphicsView->setSceneRect(0, 0, gvSize.width(), gvSize.height());
//    qDebug() << "Scene Rect: " << ui.graphicsView->sceneRect();
    ui.graphicsView->setScene(&gs);
    ui.graphicsView->setRenderHints(QPainter::Antialiasing | QPainter::SmoothPixmapTransform);
    ui.graphicsView->scale(4, 4);

    addImage(tr("/data/sliceImages/12.T2.slice.nrrd"));
    addImage(tr("/data/sliceImages/13.T2.slice.nrrd")); 
}

MainWindow::~MainWindow() {
}

void MainWindow::on_actionOpen_triggered() {
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), "/data", tr("Images (*.nrrd)"));
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
                pointList(k, 2*(j*nX + i)) = xCoords[i];
                pointList(k, 2*(j*nX + i) + 1) = yCoords[j];
            }
        }
    }

    SurfaceEntropyCostFunction<2>::Pointer costFunc = SurfaceEntropyCostFunction<2>::New();

    
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
        for (int j = 0; j < pointList.n_cols / 2; j++) {
            double x = pointList(0, 2*j);
            double y = pointList(0, 2*j+1);
            gs.addEllipse(x+1, y+1, 1, 1, QPen(Qt::red), QBrush(Qt::red, Qt::SolidPattern));
        }
    }
}

void MainWindow::on_actionRun_triggered() {
    if (bitmapList.size() == 0) {
        return;
    }

    ImageParticlesAlgorithm<BitmapType> alg(bitmapList, pointList);
    alg.Run();
}

void MainWindow::resizeEvent(QResizeEvent* event) {
}

void MainWindow::addImage(QString filename) {
    loadImage(filename);
    ui.listWidget->addItem(filename);
    ui.listWidget->setCurrentRow(ui.listWidget->count() - 1);
}