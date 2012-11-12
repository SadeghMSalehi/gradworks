//
//  imageParticleTools.cpp
//  imageParticles
//
//  Created by Joohwi Lee on 11/12/12.
//
//
#include "imageParticleTypes.h"

void SaveListOfPointVectors(ListOfPointVectorType& inputPoints, const char* fileName) {

}

void LoadListOfPointVectors(const char* fileName, ListOfPointVectorType& outputPoints) {
    
}

void ConvertParametersToListOfPointVectors(OptimizerParameters& inputParams, int nSubj, int nVars, ListOfPointVectorType& outputPoints) {
    int k = 0;
    outputPoints.clear();
    for (int j = 0; j < nSubj; j++) {
        PointVectorType subjectParticles;
        subjectParticles.reserve(nVars);
        for (int i = 0; i < nVars; i++) {
            subjectParticles.push_back(inputParams[k++]);
        }
        outputPoints.push_back(subjectParticles);
    }
}


void ConvertListOfPointVectorsToParameters(ListOfPointVectorType& inputPoints, OptimizerParameters& outputParams) {
    if (inputPoints.size() == 0) {
        return;
    }
    int nSubj = inputPoints.size();
    int nVars = inputPoints[0].size();
    outputParams.SetSize(nSubj*nVars);
    int k = 0;
    for (int j = 0; j < nSubj; j++) {
        for (int i = 0; i < nVars; i++) {
            outputParams[k++] = inputPoints[j][i];
        }
    }
}
