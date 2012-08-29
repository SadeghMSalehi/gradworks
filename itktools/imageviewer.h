#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

#include <QtGui/QWidget>
#include "QPaintEvent"

#include "ui_imageviewer.h"

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

private:
    Ui::ImageViewerClass ui;
    bool _draw;
};

#endif // IMAGEVIEWER_H
