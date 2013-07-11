//
//  mainwindow.cpp
//  laplacePDE
//
//  Created by Joohwi Lee on 11/13/12.
//
//

#include "mainwindow.h"
#include "armadillo"

#include "qutils.h"
#include <QGraphicsPixmapItem>


#define FLOAT_MAX 1e4

using namespace std;

static arma::Mat<float> _data, _dataT;
static arma::Mat<float> _mapid;
static arma::Mat<float> _imageBuffer;
static QImage _image;
static QFloatImageItem* _mapItem;
static QGraphicsPathItem* _traceItem;
static QPainterPath _tracePath;
static QGraphicsEllipseItem* _plotMarker;

uint __hsv100[100] = {
    0xFFFF0000,0xFFFF0F00,0xFFFF1F00,0xFFFF2E00,0xFFFF3D00,0xFFFF4D00,0xFFFF5C00,0xFFFF6B00,0xFFFF7A00,0xFFFF8A00,0xFFFF9900,0xFFFFA800,0xFFFFB800,0xFFFFC700,0xFFFFD600,0xFFFFE500,0xFFFFF500,0xFFFAFF00,0xFFEBFF00,0xFFDBFF00,0xFFCCFF00,0xFFBDFF00,0xFFADFF00,0xFF9EFF00,0xFF8FFF00,0xFF80FF00,0xFF70FF00,0xFF61FF00,0xFF52FF00,0xFF42FF00,0xFF33FF00,0xFF24FF00,0xFF14FF00,0xFF05FF00,0xFF00FF0A,0xFF00FF19,0xFF00FF29,0xFF00FF38,0xFF00FF47,0xFF00FF57,0xFF00FF66,0xFF00FF75,0xFF00FF85,0xFF00FF94,0xFF00FFA3,0xFF00FFB3,0xFF00FFC2,0xFF00FFD1,0xFF00FFE0,0xFF00FFF0,0xFF00FFFF,0xFF00F0FF,0xFF00E0FF,0xFF00D1FF,0xFF00C2FF,0xFF00B2FF,0xFF00A3FF,0xFF0094FF,0xFF0085FF,0xFF0075FF,0xFF0066FF,0xFF0057FF,0xFF0047FF,0xFF0038FF,0xFF0029FF,0xFF0019FF,0xFF000AFF,0xFF0500FF,0xFF1400FF,0xFF2400FF,0xFF3300FF,0xFF4200FF,0xFF5200FF,0xFF6100FF,0xFF7000FF,0xFF8000FF,0xFF8F00FF,0xFF9E00FF,0xFFAD00FF,0xFFBD00FF,0xFFCC00FF,0xFFDB00FF,0xFFEB00FF,0xFFFA00FF,0xFFFF00F5,0xFFFF00E6,0xFFFF00D6,0xFFFF00C7,0xFFFF00B8,0xFFFF00A8,0xFFFF0099,0xFFFF008A,0xFFFF007A,0xFFFF006B,0xFFFF005C,0xFFFF004D,0xFFFF003D,0xFFFF002E,0xFFFF001F,0xFFFF000F
};


MainWindow::MainWindow(QWidget* parent) {
    ui.setupUi(this);
    ui.graphicsView->setScene(&_scene);

    _mapid.load("./MapID.mat");
    cout << _mapid.n_cols << "," << _mapid.n_rows << endl;


    connect(ui.slider, SIGNAL(valueChanged(int)), SLOT(selectValue(int)));
    connect(ui.clearGraphButton, SIGNAL(pressed()), SLOT(clearGraphs()));
    connect(ui.loadSeries, SIGNAL(pressed()), SLOT(loadSeries()));

    _mapItem = new QFloatImageItem();
    _mapItem->setImage(_imageBuffer.memptr(), _imageBuffer.n_rows, _imageBuffer.n_cols);
    _mapItem->setRange(ui.minValue->value(), ui.maxValue->value());
    _mapItem->setInteraction(this);

    _scene.addItem(_mapItem);
    _traceItem = _scene.addPath(_tracePath);


    _traceItem->setPen(QPen(Qt::green));

    QTransform matrix;
    matrix.scale(5, 5);
    ui.graphicsView->setTransform(matrix);

    ui.customPlot->setInteraction(QCustomPlot::iSelectLegend, true);
    ui.customPlot->setInteraction(QCustomPlot::iSelectPlottables, true);
    connect(ui.customPlot, SIGNAL(selectionChangedByUser()), SLOT(graphSelected()));

    _plotMarker = new QGraphicsEllipseItem(_mapItem);
    _plotMarker->setVisible(false);
}

MainWindow::~MainWindow() {

}

void MainWindow::loadSeries() {
    QString seriesFile = __fileManager.openFile(0, this, "Open Series File...", "*.mat");

    if (seriesFile != "") {
        cout << "Loading start ..." << flush;
        _data.load(seriesFile.toStdString().c_str());
        cout << " done. " << endl;

        _dataT = arma::trans(_data);
        cout << _data.n_cols << "," << _data.n_rows << endl;

        _imageBuffer.set_size(_mapid.n_rows, _mapid.n_cols);
        for (int i = 0; i < _mapid.n_rows; i++) {
            for (int j = 0; j < _mapid.n_cols; j++) {
                if (_mapid(i,j) != _mapid(i,j)) {
                    _imageBuffer(i,j) = NAN;
                } else {
                    _imageBuffer(i,j) = 0;
                }
            }
        }

        ui.slider->setMaximum(_data.n_cols - 1);

        ui.minValue->setValue(_data.min());
        ui.maxValue->setValue(_data.max());

        _mapItem->setImage(_imageBuffer.memptr(), _imageBuffer.n_rows, _imageBuffer.n_cols);
        _mapItem->setRange(ui.minValue->value(), ui.maxValue->value());
        _mapItem->setInteraction(this);

        selectValue(0);
    }
}

void MainWindow::selectValue(int n) {
    ui.yearMonth->setText(QString("%1/%2").arg(n / 12 + ui.startingYear->value()).arg(n % 12 + 1, 2, 10, QChar('0')));

    int k = 0;
    for (int i = 0; i < _imageBuffer.n_rows; i++) {
        for (int j = 0; j < _imageBuffer.n_cols; j++) {

            if (std::isnan(_imageBuffer(i,j))) {
                continue;
            }
            float w = _data(k,n);
            _imageBuffer(i, j) = w;
            k++;
        }
    }

    _mapItem->setRange(ui.minValue->value(), ui.maxValue->value());
    _mapItem->refresh();
}


float MainWindow::crossCorrelation(float* x, float* y, int n) {
    int i,j;
    double mx,my,sx,sy,sxy,denom,r;

    /* Calculate the mean of the two series x[], y[] */
    mx = 0;
    my = 0;
    for (i=0;i<n;i++) {
        mx += x[i];
        my += y[i];
    }
    mx /= n;
    my /= n;

    /* Calculate the denominator */
    sx = 0;
    sy = 0;
    for (i=0;i<n;i++) {
        sx += (x[i] - mx) * (x[i] - mx);
        sy += (y[i] - my) * (y[i] - my);
    }
    denom = sqrt(sx*sy);

    /* Calculate the correlation series */
    sxy = 0;
    for (i=0; i<n; i++) {
        j = i;
        if (j < 0 || j >= n)
            continue;
        else
            sxy += (x[i] - mx) * (y[j] - my);
    }
    r = sxy / denom;
    return r;
}

void MainWindow::selectCorr(int idx) {
    int mapId = _mapid(idx) - 1;

    int k = 0;
    for (int i = 0; i < _imageBuffer.n_rows; i++) {
        for (int j = 0; j < _imageBuffer.n_cols; j++) {

            if (std::isnan(_imageBuffer(i,j))) {
                continue;
            }
            float corr = crossCorrelation(_dataT.colptr(mapId), _dataT.colptr(k), _dataT.n_rows);
            _imageBuffer(i, j) = corr;
            k++;
        }
    }

    _mapItem->setRange(-1, 1);
    _mapItem->refresh();
}

void MainWindow::clearGraphs() {
    ui.customPlot->clearGraphs();
    ui.customPlot->replot();

    _tracePath = QPainterPath();
    _traceItem->setPath(_tracePath);
}

void MainWindow::graphSelected() {
    QList<QCPGraph*> graphs = ui.customPlot->selectedGraphs();
    for (int i = 0; i < graphs.size(); i++) {
        int x = graphs[i]->property("x").value<int>();
        int y = graphs[i]->property("y").value<int>();

        _plotMarker->setPos(QPointF(x - 0.5f, y - 0.5f));
        _plotMarker->setRect(0, 0, 1, 1);
        _plotMarker->setPen(QPen(Qt::magenta));
        _plotMarker->setVisible(true);

        int graphid = graphs[i]->property("graphid").value<int>();
        ui.customPlot->legend->item(graphid)->setSelected(true);
    }
    if (graphs.size() == 0) {
        _plotMarker->setVisible(false);
        for (int i = 0; i < ui.customPlot->legend->itemCount(); i++) {
            ui.customPlot->legend->item(i)->setSelected(false);
        }
    }
}

void MainWindow::mousePressed(QFloatImageItem* self, QGraphicsSceneMouseEvent* event) {
    QPoint pos = event->pos().toPoint();

    int lg = pos.x();
    int lt = pos.y();

    if (lg < 91) {
        lg = -180 + (91 - lg) * 2;
    } else {
        lg = 180 - (lg - 91) * 2;
    }

    if (lt < 45) {
        lt = (45 - lt) * 2;
    } else {
        lt = -(lt - 45) * 2;
    }

    cout << "Long: " << lg << " / Lati: " << lt << endl;

    if (ui.showCorr->isChecked()) {
        int idx = pos.y() * _imageBuffer.n_rows + pos.x();
        selectCorr(idx);

        clearGraphs();
    } else if (ui.showMonth->isChecked()) {

        if (_tracePath.elementCount() == 0) {
            _tracePath.moveTo(pos);
        } else {
            _tracePath.lineTo(pos);
        }
        _traceItem->setPath(_tracePath);

        int idx = pos.y() * _imageBuffer.n_rows + pos.x();

        ui.customPlot->xAxis->setLabel("month");
        ui.customPlot->yAxis->setLabel("Temp");

        QVector<double> xTicks;
        QVector<QString> xLabels;
        for (int i = 0; i < _dataT.n_rows; i += _dataT.n_rows / 10) {
            xTicks << i;
            xLabels << QString("%1/%2").arg(i / 12 + ui.startingYear->value()).arg(i % 12 + 1, 2, 10, QChar('0'));
        }

        ui.customPlot->xAxis->setAutoTicks(false);
        ui.customPlot->xAxis->setAutoTickStep(false);
        ui.customPlot->xAxis->setAutoTickLabels(false);
        ui.customPlot->xAxis->setTickVector(xTicks);
        ui.customPlot->xAxis->setTickVectorLabels(xLabels);


        int mapId = _mapid(idx) - 1;
        float *colData = _dataT.colptr(mapId);

        QVector<double> xdata, ydata;
        for (int i = 0; i < _dataT.n_rows; i++) {
            xdata.append(i);
            ydata.append(colData[i]);
        }

        ui.customPlot->addGraph(ui.customPlot->xAxis, ui.customPlot->yAxis);

        int lastGraph = ui.customPlot->graphCount() - 1;
        ui.customPlot->graph(lastGraph)->setData(xdata, ydata);

        QPen lineColor(__hsv100[(lastGraph) % 100] & 0xaaffffff);
        ui.customPlot->graph(lastGraph)->setPen(lineColor);

        ui.customPlot->graph(lastGraph)->setName(QString("(%1,%2)").arg(lg).arg(lt));
        ui.customPlot->graph(lastGraph)->setSelectable(true);

        ui.customPlot->graph(lastGraph)->setProperty("x", QVariant(pos.x()));
        ui.customPlot->graph(lastGraph)->setProperty("y", QVariant(pos.y()));
        ui.customPlot->graph(lastGraph)->setProperty("mapid", QVariant(mapId));
        ui.customPlot->graph(lastGraph)->setProperty("graphid", QVariant(lastGraph));


        ui.customPlot->rescaleAxes();
        ui.customPlot->replot();
    }
}

void MainWindow::mouseMoved(QFloatImageItem* self, QGraphicsSceneMouseEvent* event) {

}

void MainWindow::mouseReleased(QFloatImageItem* self, QGraphicsSceneMouseEvent* event) {

}

void MainWindow::hoverEntered(QFloatImageItem* self, QGraphicsSceneHoverEvent* event) {

}

void MainWindow::hoverMoved(QFloatImageItem* self, QGraphicsSceneHoverEvent* event) {

}

void MainWindow::hoverLeft(QFloatImageItem* self, QGraphicsSceneHoverEvent* event) {

}

void MainWindow::on_actionShowPlotWindow_triggered() {
    ui.dockWidget->setVisible(true);
}

void MainWindow::on_showPlotLegend_toggled(bool plotLegend) {
    ui.customPlot->legend->setVisible(plotLegend);
    ui.customPlot->replot();
}