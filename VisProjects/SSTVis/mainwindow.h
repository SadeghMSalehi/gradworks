#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui/QMainWindow>
#include <QGraphicsScene>
#include <QGraphicsImageItem.h>
#include "ui_mainwindow.h"

typedef QGraphicsImageItem<float> QFloatImageItem;


class MainWindow: public QMainWindow, public QFloatImageItem::InteractionType {
    Q_OBJECT
public:
    MainWindow(QWidget* parent = NULL);
    ~MainWindow();

    float crossCorrelation(float* col1, float* col2, int nRows);
    
public slots:
    void selectValue(int n);
    void selectCorr(int idx);
    void clearGraphs();
    void loadSeries();

    void graphSelected();
    
    void on_actionShowPlotWindow_triggered();
    void on_showPlotLegend_toggled(bool);
    
protected:
    virtual void mousePressed(QFloatImageItem* self, QGraphicsSceneMouseEvent* event);
    virtual void mouseMoved(QFloatImageItem* self, QGraphicsSceneMouseEvent* event);
    virtual void mouseReleased(QFloatImageItem* self, QGraphicsSceneMouseEvent* event);
    virtual void hoverEntered(QFloatImageItem* self, QGraphicsSceneHoverEvent* event);
    virtual void hoverMoved(QFloatImageItem* self, QGraphicsSceneHoverEvent* event);
    virtual void hoverLeft(QFloatImageItem* self, QGraphicsSceneHoverEvent* event);


private:
    Ui_MainWindow ui;
    QGraphicsScene _scene;

};

#endif
