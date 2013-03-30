//
//  airAlgorithmManager.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/29/13.
//
//

#ifndef __ParticleGuidedRegistration__airAlgorithmManager__
#define __ParticleGuidedRegistration__airAlgorithmManager__

#include <iostream>
#include "ui_pviewAIRWindow.h"
#include "piImageSlice.h"

class AIRWindow;

namespace air {

enum AlgorithmType { KMEANS };
enum AlgorithmResult { SUCCESS, FAIL };

class AlgorithmManager: public QObject {
    Q_OBJECT

private:
    AIRWindow* _window;
    Ui::AIRWindow& ui;

    void setupConnection();
    
public:
    AlgorithmManager(AIRWindow* parent, Ui::AIRWindow& u);
    virtual ~AlgorithmManager();

public slots:
    void executeKmeans2D();

signals:
    void algorithmFinished(AlgorithmType, AlgorithmResult, pi::AIRLabel::Pointer output);

};

}

#endif /* defined(__ParticleGuidedRegistration__airAlgorithmManager__) */