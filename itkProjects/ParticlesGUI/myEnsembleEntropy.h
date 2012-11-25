//
//  myEnsembleEntropy.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/24/12.
//
//

#ifndef __ParticlesGUI__myEnsembleEntropy__
#define __ParticlesGUI__myEnsembleEntropy__

#include <iostream>

#include "myImageContainer.h"
#include "itkOptimizerCommon.h"
#include "vnlCommon.h"

class myEnsembleEntropy {
public:
    myEnsembleEntropy();
    ~myEnsembleEntropy();

    void SetImageList(ImageContainer::List* imageList);
    void SetInitialPositions(const OptimizerParametersType& params, int nSubject, int nPoints, int nParams);
    void ComputePositionalEnsemble(const OptimizerParametersType& params, double& cost, VNLMatrix& deriv) const;

    void SetPatchSize(int n) { m_PatchSize = n; }
    void SetTransformType(int t) { m_TransformType = t; }
    void SetTransformTypeToRigid() { m_TransformType = 1; }
    void EstimateRigidParameters(VNLMatrix& transformParams, const OptimizerParametersType& params, int target, int source) const;
    void SetVariableCounts(int s, int p, int np);
private:
    //void EstimateRigidParameters(VNLMatrix& transformParams, const OptimizerParametersType& params, int target, int source) const;
    
    ImageContainer::List* m_ImageList;
    int m_PatchSize;
    int m_TransformType;
    int m_nSubjects;
    int m_nPoints;
    int m_nParams;
    
    // assume there's only rigid transformations
    const int m_nTransformParams = 3;

};

#endif /* defined(__ParticlesGUI__myEnsembleEntropy__) */
