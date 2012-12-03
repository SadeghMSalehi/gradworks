#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

#include <QtGui/QWidget>
#include "QPaintEvent"

#include "ui_imageviewer.h"
#include "itkMyCore.h"

class ImageViewer : public QWidget
{
    Q_OBJECT

public:
    ImageViewer(QWidget *parent = 0);
    ~ImageViewer();

public slots:
	void toggleDraw(void);

protected:
    void paintEvent(QPaintEvent* event);
    void drawSample(QPaintEvent* event);
    void drawSlice(QPaintEvent* event);
private:
    Ui::ImageViewerClass ui;
    bool _draw;
    itkMyCore* _core;
};

#endif // IMAGEVIEWER_H
