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

public slots:
    void openFile(QString);
    void zoomIn();
    void zoomOut();

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
    Ui::BigViewWindow ui;
    pi::AIRDisplayCollection* _images;
    pi::SliceDirectionEnum _sliceDirection;
};
#endif /* defined(__ParticleGuidedRegistration__bigViewWindow__) */
