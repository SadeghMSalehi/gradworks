//
//  airAlgorithmManager.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/29/13.
//
//

#include "airAlgorithmManager.h"
#include "pviewAIRWindow.h"

using namespace pi;

namespace air {
    AlgorithmManager::AlgorithmManager(AIRWindow* parent, Ui::AIRWindow& u): _window(parent), ui(u)  {
        setupConnection();
    }

    AlgorithmManager::~AlgorithmManager() {
        
    }

    void AlgorithmManager::setupConnection() {
        connect(ui.kMeansRun, SIGNAL(clicked()), this, SLOT(executeKmeans2D()));
    }

    void AlgorithmManager::executeKmeans2D() {
        if (!_window->IsImage1Loaded()) {
            emit algorithmFinished(KMEANS, FAIL, NULL);
            return;
        }

        AIRImageDisplay gray = _window->imageDisplays[0];
        AIRLabelSlice label = ui.graphicsView->getLabelSlice();

        cout << "executeKmeans()" << endl;
    }
}