//
//  piParticleTrainer.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 7/21/13.
//
//

#include <itkLabelStatisticsImageFilter.h>

#include "piParticleTrainer.h"
#include "piParticleCollision.h"

using namespace pi;

ParticleTrainer::ParticleTrainer(ParticleSystem* system) {
    _system = system;
}

ParticleTrainer::~ParticleTrainer() {

}


// warning: must set number of particles before run this
void ParticleTrainer::trainParticles() {
    sampleParticles();
    for (int i = 0; i < _numLabels; i++) {
        spreadParticles(i);
    }
}

void ParticleTrainer::sampleParticles() {
    _numLabels = 0;
    ParticleSystem system = *_system;


    int nPoints = 0;
    for (int j = 0; j < _counts.size(); j++) {
        nPoints += _counts[j];
    }

    _numSubjs = _system->GetNumberOfSubjects();
    for (int i = 0; i < _numSubjs; i++) {
        LabelImage::Pointer labelImage = system[i].GetLabel();

        typedef itk::LabelStatisticsImageFilter<RealImage, LabelImage> LabelStatFilterType;
        LabelStatFilterType::Pointer labelStat = LabelStatFilterType::New();
        labelStat->SetInput(system[i].GetImage(0));
        labelStat->SetLabelInput(labelImage);
        labelStat->Update();
        
        _numLabels = labelStat->GetNumberOfLabels();
        
        LabelPoints indices;
        indices.resize(_numLabels);

        for (int j = 0; j < _numLabels; j++) {
            indices[j].reserve(1000);
        }

        LabelImageIteratorType iter(labelImage, labelImage->GetBufferedRegion());
        iter.GoToBegin();

        LabelImage::PointType indexPoint;
        while (!iter.IsAtEnd()) {
            int label = iter.Get();
            if (label > 0) {
                labelImage->TransformIndexToPhysicalPoint(iter.GetIndex(), indexPoint);
                indices[label].push_back(indexPoint);
            }
            ++iter;
        }

        for (int j = 0; j < _numLabels; j++) {
            std::random_shuffle(indices[j].begin(), indices[j].end());
        }

        system[i].NewParticles(nPoints);

        int l = 0;
        for (int j = 0; j < _numLabels; j++) {
            for (int k = 0; k < _counts[j]; k++) {
                fordim (m) {
                    system[i][l++].x[m] = indices[j][k][m];
                }
            }
        }
    }
}

void ParticleTrainer::spreadParticles(int labelId) {
    for (int i = 0; i < _numSubjs; i++) {
        ParticleCollision collisionHandler;
        ParticleSubject& subj = _system->GetSubjects()[i];
        for (int j = 0; j < _numLabels; j++) {
            LabelImage::Pointer extractedLabel = extractLabel(subj.GetLabel(), labelId);
            collisionHandler.subject = &subj;
            collisionHandler.SetLabelImage(extractedLabel);
            collisionHandler.UpdateImages();
        }
    }
}


LabelImage::Pointer ParticleTrainer::extractLabel(LabelImage::Pointer labelImage, int labelId) {
    
}