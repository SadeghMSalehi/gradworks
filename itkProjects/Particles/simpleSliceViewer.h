//
//  simpleSliceViewer.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/2/13.
//
//

#ifndef __ParticleGuidedRegistration__simpleSliceViewer__
#define __ParticleGuidedRegistration__simpleSliceViewer__

#include <iostream>
#include <QMainWindow>
#include <QGraphicsScene>
#include <itkImage.h>
#include "itkARGBSliceExtractImageFilter.h"
#include "ui_simpleSliceViewer.h"

class SimpleSliceViewer: public QMainWindow {
    Q_OBJECT

public:
    typedef itk::Image<double,3> ImageType;
    typedef ImageType::Pointer ImagePointer;

    SimpleSliceViewer(QWidget* parent = NULL);
    virtual ~SimpleSliceViewer();

public slots:
    void openImage(QString imageFile);
    void updateDisplay();

protected:
    
    QGraphicsScene _scene;
    Ui::SimpleSliceViewer ui;
    itk::ARGBSliceExtractImageFilter<ImageType>::Pointer _sliceExtractor;
};
#endif /* defined(__ParticleGuidedRegistration__simpleSliceViewer__) */
