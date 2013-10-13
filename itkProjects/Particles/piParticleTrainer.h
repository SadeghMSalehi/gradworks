//
//  piParticleTrainer.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 7/21/13.
//
//

#ifndef __ParticleGuidedRegistration__piParticleTrainer__
#define __ParticleGuidedRegistration__piParticleTrainer__

#include <iostream>
#include <vector>
#include <itkExtractImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>

#include "piParticleCore.h"

namespace pi {

typedef std::pair<int, double> ParticleWeightPair;
typedef std::vector<ParticleWeightPair> ParticleWeightVector;
typedef std::pair<int, ParticleWeightVector> NeighborWeightPair;
typedef std::vector<NeighborWeightPair> NeighborWeightMatrix;
typedef std::vector<NeighborWeightMatrix> NeighborCollection;

class ParticleTrainer {
public:
    ParticleTrainer(ParticleSystem* system);
    ~ParticleTrainer();

    void setNumberOfParticles(IntVector& counts);
    void trainParticles();
    void saveParticles();
    void loadParticles();
    void initialClosestCorrespondences();
    void establishCorrespondences();
    void sampleParticles();


    // remove outliers from graph-based correspondences
    // outlier check is eulidean-distance based
    void removeOutliers();

    void computeNearestNeighbors(NeighborWeightMatrix& neighbors, int subjId, int labelId);
    void computeDelaunayTriangulation(int subjId, int labelId, int pointSetId);

    void setNumberOfPointSets(int n);
    NeighborWeightMatrix& getPointSets(int n);

private:
    typedef itk::LabelStatisticsImageFilter<RealImage, LabelImage> LabelStatFilterType;
    typedef itk::ExtractImageFilter<LabelImage, LabelImage> LabelExtractFilterType;
    typedef std::vector<LabelStatFilterType::BoundingBoxType> BoundingBoxVector;
    typedef std::vector<BoundingBoxVector> LabelBoundingBoxes;
    typedef std::vector<LabelImage::PointType> PointVector;
    typedef std::vector<PointVector> LabelPoints;

    void spreadParticles(int labelId);
    LabelImage::Pointer extractLabel(LabelImage::Pointer labelImage, int imageId, int labelId);

private:
    ParticleSystem* _system;
    int _numLabels;
    int _numSubjs;
    IntVector _counts;
    LabelBoundingBoxes _labelBounds;
    vnl_matrix<int> _incidenceGraph;
    NeighborCollection _neighborCollection;
};

}
#endif /* defined(__ParticleGuidedRegistration__piParticleTrainer__) */
