//
//  QGraphicsParticleItems.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/28/13.
//
//

#include "QGraphicsParticleItems.h"
#include <QGraphicsScene>
#include "itkARGBColorFunction.h"
#include <itkRGBAPixel.h>
#include "piParticleCore.h"

using namespace pi;

// create particle item
typedef itk::RGBAPixel<unsigned char> RGBA;
typedef itk::Function::HSVColormapFunction<float, RGBA> HSVFunction;

HSVFunction::Pointer _hsvFunc = HSVFunction::New();


QGraphicsParticleItems::QGraphicsParticleItems() {
    _scene = NULL;
    _parentItem = NULL;
    _subject = NULL;
    _refSubject = NULL;
    _isParticleSelected = false;
    _particleSelectedId = -1;
    _listener = NULL;
    _useParticleX = true;
}

QGraphicsParticleItems::~QGraphicsParticleItems() {

}


void QGraphicsParticleItems::useParticleX(bool value) {
    _useParticleX = value;
}

int QGraphicsParticleItems::getSelectedParticleId() {
    if (_isParticleSelected) {
        return _particleSelectedId;
    } else {
        return -1;
    }
}

void QGraphicsParticleItems::hideParticles(bool on) {
    _isHide = on;
    updateParticles();
}

void QGraphicsParticleItems::setListener(QGraphicsEllipseEventItem::Listener *listener) {
    _listener = listener;

    const int n = _particleItems.size();
    for (int j = 0; j < n; j++) {
        _particleItems[j]->setListener(_listener);
    }
}

void QGraphicsParticleItems::setScene(QGraphicsScene* scene) {
    _scene = scene;
}

void QGraphicsParticleItems::setParentItem(QGraphicsItem *parentItem) {
    _parentItem = parentItem;
}

void QGraphicsParticleItems::setReferenceSubject(pi::ParticleSubject *refSubject) {
    _refSubject = refSubject;
}

void QGraphicsParticleItems::useScalars(bool use) {
    _useScalars = use;
}

void QGraphicsParticleItems::setScalars(VNLVector scalars) {
    _scalars = scalars;
}

void QGraphicsParticleItems::createParticles(pi::ParticleSubject *subject) {
    _subject = subject;

    clearParticles();

    const int n = _subject->GetNumberOfPoints();
    _hsvFunc->SetMinimumInputValue(0);
    _hsvFunc->SetMaximumInputValue(n);
//    _hsvFunc->SetMaximumInputValue(5);
    _particleItems.resize(n);

    for (int j = 0; j < n; j++) {
        const double r = 5;
        RGBA color = _hsvFunc->operator()(j);
        _particleItems[j] = new QGraphicsEllipseEventItem(_parentItem);
        _particleItems[j]->setListener(_listener);
        _particleItems[j]->setData(0, QVariant(j));
        _particleItems[j]->setRect(-r/2.0, -r/2.0, r, r);
        _particleItems[j]->setZValue(10);
        _particleItems[j]->setOpacity(0.3);
        _particleItems[j]->setPen(Qt::NoPen);
        _particleItems[j]->setBrush(QBrush(qRgb(color[0], color[1], color[2]), Qt::SolidPattern));
//        _particleItems[j]->setBrush(QBrush(Qt::blue, Qt::SolidPattern));
        if (_isHide) {
            _particleItems[j]->hide();
        }
    }
}

void QGraphicsParticleItems::updateParticles() {
    if (_subject == NULL) {
        return;
    }

    int numberOfParticleItems = _particleItems.size();
    const int numberOfSubjectPoints = _subject->GetNumberOfPoints();
    
    if (numberOfParticleItems != numberOfSubjectPoints) {
        createParticles(_subject);
        numberOfParticleItems = numberOfSubjectPoints;
    }

    bool useColorFromScalar = _useScalars && _scalars.size() == numberOfSubjectPoints;
    if (useColorFromScalar) {
        _hsvFunc->SetMinimumInputValue(_scalars.min_value());
        _hsvFunc->SetMaximumInputValue(_scalars.max_value());
    }

    assert(_particleItems.size() == _subject->GetNumberOfPoints());
    
    QBrush grayBrush(Qt::gray, Qt::SolidPattern);
    for (int j = 0; j < numberOfParticleItems; j++) {
        Particle& p = _subject->operator[](j);
        if (_isHide) {
            _particleItems[j]->hide();
        } else {
            _particleItems[j]->show();
        }
        
        RealIndex x;
        if (_useParticleX) {
            _subject->ComputeIndexX(p, x);
        } else {
            _subject->ComputeIndexZ(p, x);
        }
        _particleItems[j]->setPos(x[0], x[1]);

        if (p.enabled) {
            if (_isParticleSelected && j != _particleSelectedId) {
                _particleItems[j]->setBrush(grayBrush);
                _particleItems[j]->setOpacity(0.3);
            } else {
                RGBA color;
                if (useColorFromScalar) {
                    color = _hsvFunc->operator()(_scalars[j]);
                } else {
                    color = _hsvFunc->operator()(j);
                }
                if (_refSubject != NULL && _refSubject->operator[](j).outlier) {
                    _particleItems[j]->setBrush(QBrush(Qt::white, Qt::SolidPattern));
                    _particleItems[j]->setOpacity(1);
                } else {
                    _particleItems[j]->setBrush(QBrush(qRgb(color[0], color[1], color[2]), Qt::SolidPattern));
                    _particleItems[j]->setOpacity(0.3);
                }
            }
        } else {
            _particleItems[j]->hide();
        }
    }
}

void QGraphicsParticleItems::clearParticles() {
    if (_subject == NULL) {
        return;
    }

    // check if need to create particle item
    const int n = _particleItems.size();
    for (int j = 0; j < n; j++) {
        _scene->removeItem(_particleItems[j]);
    }
}

void QGraphicsParticleItems::selectParticle(int particleId) {
    if (_isParticleSelected) {
        if (particleId == _particleSelectedId) {
            _isParticleSelected = false;
            _particleSelectedId = -1;
        } else {
            _particleSelectedId = particleId;
        }
    } else {
        _particleSelectedId = particleId;
        _isParticleSelected = true;
    }

    updateParticles();
}

QGraphicsEllipseEventItem* QGraphicsParticleItems::getItem(int i) {
    return _particleItems[i];
}
