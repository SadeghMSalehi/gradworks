//
//  airSuperZoom.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/30/13.
//
//

#ifndef __ParticleGuidedRegistration__airSuperZoom__
#define __ParticleGuidedRegistration__airSuperZoom__

#include <iostream>
#include "ui_testDialog.h"
#include <QGraphicsScene>

namespace air {
    class SuperZoom: public QDialog {
        Q_OBJECT
    public:
        SuperZoom(QWidget* p);
        virtual ~SuperZoom();

    public slots:
        void actionTriggered();

    protected:
        void showEvent(QShowEvent* event);
        void closeEvent(QCloseEvent* event);
        
    private:
        Ui::Dialog ui;
        QGraphicsScene _scene;
    };
};
#endif /* defined(__ParticleGuidedRegistration__airSuperZoom__) */
