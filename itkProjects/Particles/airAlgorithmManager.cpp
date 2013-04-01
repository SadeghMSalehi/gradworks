//
//  airAlgorithmManager.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/29/13.
//
//

#include "airAlgorithmManager.h"
#include "airImageAlgorithm.h"
#include "pviewAIRWindow.h"
#include <itkIsolatedConnectedImageFilter.h>


using namespace pi;

typedef air::ImageAlgorithm<AIRImage,AIRLabel> Algo;
Algo __algo;

namespace air {
    AlgorithmManager::AlgorithmManager(AIRWindow* parent, Ui::AIRWindow& u): _window(parent), ui(u)  {
        setupConnection();
    }

    AlgorithmManager::~AlgorithmManager() {
        
    }

    void AlgorithmManager::setupConnection() {
        connect(ui.kMeansRun, SIGNAL(clicked()), this, SLOT(executeKmeans2D()));
        connect(ui.isolatedRegionGrowingRun, SIGNAL(clicked()), this, SLOT(executeIsoRG()));
    }

    void AlgorithmManager::executeKmeans2D() {
        if (!_window->IsImage1Loaded()) {
            emit algorithmFinished(KMEANS, NORUN, NULL);
            return;
        }

        AIRImageDisplay gray = _window->imageDisplays[0];
        AIRLabelSlice label = ui.graphicsView->getLabelSlice();

        cout << "executeKmeans()" << endl;
    }

    void AlgorithmManager::executeIsoRG() {
        if (!_window->IsImage1Loaded()) {
            emit algorithmFinished(KMEANS, NORUN, NULL);
            return;
        }

        AIRClass fgId = _window->getForegroundLabel();
        AIRClass bgId = _window->getBackgroundLabel();

        AIRImageDisplay gray = _window->imageDisplays[0];
        AIRImage::Pointer srcImg = gray.srcImg;
        AIRLabel::Pointer labelImg = ui.graphicsView->getLabelVolume();
        AIRLabel::Pointer labelOutput = __algo.ExecuteIsolatedConnectedImageFilter(labelImg, srcImg, fgId, bgId);

        __airLabelIO.WriteImage("/tmpfs/temp.nrrd", labelOutput);
    }
}