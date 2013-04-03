//
//  bigViewWindow.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/31/13.
//
//

#ifndef __ParticleGuidedRegistration__bigViewWindow__
#define __ParticleGuidedRegistration__bigViewWindow__

#include <iostream>

#include "ui_bigViewWindow.h"
#include <QGraphicsScene>

class QDoubleSpinBox;

namespace pi {
    template<class T>
    class ImageDisplayCollection;
}

class BigViewWindow: public QMainWindow {
    Q_OBJECT
public:
    BigViewWindow(QWidget* parent = NULL);
    virtual ~BigViewWindow();

signals:
    void fileDropped(QString);
    void multipleFileDropeed(QList<QString>);

public slots:
    void openFile(QString fileName = "");
    void openFiles(QList<QString> files);
    void zoomIn();
    void zoomOut();
    void volumeSelected(int);
    void changeIntensity();
    void changeDirection();
    void openDualViewer();

protected:
    void dragEnterEvent(QDragEnterEvent *event);
    void dropEvent(QDropEvent *event);

private:
    void setupUi();
    void setupShortcuts();
    void connectSignalSlots();
    void initialize();
    void centerToDesktop();
    
private:
    QMovie* _loadingMovie;
    QDoubleSpinBox* _lowIntensitySpinBox;
    QDoubleSpinBox* _highIntensitySpinBox;

    Ui::BigViewWindow ui;
    pi::AIRDisplayCollection* _images;
    pi::SliceDirectionEnum _sliceDirection;
};

#endif /* defined(__ParticleGuidedRegistration__bigViewWindow__) */