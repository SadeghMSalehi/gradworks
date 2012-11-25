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

class myEnsembleEntropy {
public:
    typedef vnl_matrix<double> MatrixType;
    typedef vnl_vector<double> VectorTYpe;

    void SetImageList(const ImageContainer::List* imageList);
    void SetInitialPositions(const OptimizerParametersType& params);
    void Compute(const OptimizerParametersType& params, double& cost, MatrixType& deriv) const;

    void SetPatchSize(int n) { m_PatchSize = n; }
    void SetTransformType(int t) { m_TransformType = t; }
    void SetTransformTypeToRigid() { m_TransformType = 1; }

private:
    ImageContainer::List* m_ImageList;
    int m_PatchSize;
    int m_TransformType;

};

#endif /* defined(__ParticlesGUI__myEnsembleEntropy__) */
