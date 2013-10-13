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
#include "piParticleForces.h"
#include "piImageIO.h"

#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkDelaunay2D.h"
#include "vtkIdList.h"
#include "vtkCell.h"

#include "set"

using namespace pi;

struct ParticleWeightComparator {
    bool operator()(const ParticleWeightPair& a, const ParticleWeightPair& b) {
        return a.second < b.second;
    }
};


ParticleTrainer::ParticleTrainer(ParticleSystem* system) {
    _system = system;
}

ParticleTrainer::~ParticleTrainer() {

}

void ParticleTrainer::setNumberOfParticles(IntVector &counts) {
    _counts = counts;
}

// warning: must set number of particles before run this
void ParticleTrainer::trainParticles() {
    sampleParticles();
    for (int i = 0; i < _numLabels; i++) {
        spreadParticles(i);
    }

    saveParticles();
}

void ParticleTrainer::sampleParticles() {
    _numLabels = 0;
    ParticleSystem system = *_system;


    int nPoints = 0;
    for (int j = 0; j < _counts.size(); j++) {
        nPoints += _counts[j];
    }

    _numSubjs = _system->GetNumberOfSubjects();
    _labelBounds.resize(_numSubjs);

    // Random particle sampling for every subject
    for (int i = 0; i < _numSubjs; i++) {
        LabelImage::Pointer labelImage = system[i].GetLabel();

        LabelStatFilterType::Pointer labelStat = LabelStatFilterType::New();
        labelStat->SetInput(system[i].GetImage(0));
        labelStat->SetLabelInput(labelImage);
        labelStat->Update();

        // count # of labels
        _numLabels = labelStat->GetNumberOfLabels();

        LabelPoints labelPoints;
        labelPoints.resize(_numLabels);
        _labelBounds[i].resize(_numLabels);

        for (int j = 1; j < _numLabels; j++) {
            labelPoints[j].reserve(1000);
            _labelBounds[i][j] = labelStat->GetBoundingBox(j);
        }

        LabelImageIteratorType iter(labelImage, labelImage->GetBufferedRegion());
        iter.GoToBegin();

        LabelImage::PointType indexPoint;
        while (!iter.IsAtEnd()) {
            int label = iter.Get();
            if (label > 0) {
                labelImage->TransformIndexToPhysicalPoint(iter.GetIndex(), indexPoint);
                labelPoints[label].push_back(indexPoint);
            }
            ++iter;
        }

        // randomly sample points except for background
        for (int j = 1; j < _numLabels; j++) {
            std::random_shuffle(labelPoints[j].begin(), labelPoints[j].end());
        }

        // spread particles inside each label
        ParticleSubject& subj = _system->GetSubjects()[i];

        // create particles
        subj.NewParticles(nPoints);

        // assign random points to the particle
        int l = 0;
        for (int j = 1; j < _numLabels; j++) {
            // be careful about the label index
            for (int k = 0; k < _counts[j-1]; k++) {
                fordim (m) {
                    subj[l].x[m] = labelPoints[j][k][m];
                }
                subj[l].label = j;
                l++;
            }
        }

//        ParticleMultiCollision collisionHandler;
//        collisionHandler.subject = &subj;
//        collisionHandler.Initialize(subj.GetLabel());


        /*
        for (int j = 1; j < _numLabels; j++) {
            LabelImage::Pointer extractedLabel = extractLabel(subj.GetLabel(), i, j);

            ParticleCollision collisionHandler;
            collisionHandler.subject = &subj;

            char fname[128];
            sprintf(fname, "/tmpfs/mask_%02d_%02d.nrrd", i, j);
            collisionHandler.binaryMaskCache = string(fname);
            sprintf(fname, "/tmpfs/dmap_%02d_%02d.nrrd", i, j);
            collisionHandler.distanceMapCache = string(fname);

            collisionHandler.SetLabelImage(extractedLabel);
            collisionHandler.UpdateImages();

            // run again with preprocessing interval
            // give enough time range
            DataReal t0 = 0;
            DataReal dt = 0.05;
            DataReal t1 = 500;

            EntropyInternalForce internalForce;

            // iterate over
            for (DataReal t = t0; t < t1; t += dt) {
                cout << "processing time: " << t << endl;

                // clear forces
                for (int k = 0; k < nPoints; k++) {
                    Particle& pi = subj[k];
                    if (pi.label != j) {
                        continue;
                    }
                    forfill(pi.f, 0);
                }
                
                // particle force computation
                collisionHandler.ConstrainPoint(subj, j);
                internalForce.ComputeForce(subj, j);
                collisionHandler.ProjectForceAndVelocity(subj, j);

                // update system
                for (int l = 0; l < nPoints; l++) {
                    Particle& p = subj[l];
                    if (p.label != j) {
                        continue;
                    }
                    p.UpdateStatus(dt);

                    IntIndex pIdx;
                    fordim (k) {
                        pIdx[k] = p.x[k] + 0.5;
                    }
                    if (!collisionHandler.IsBufferInside(pIdx)) {
                        cout << "\nStop system: out of region" << endl;
                        return;
                    }
                }
            }
        }
         */
    }
}


void ParticleTrainer::spreadParticles(int labelId) {
}


LabelImage::Pointer ParticleTrainer::extractLabel(LabelImage::Pointer labelImage, int imageId, int labelId) {
    LabelImage::RegionType wholeRegion;
    LabelImage::RegionType extractRegion;
    LabelImage::IndexType idx1, idx2;

    wholeRegion = labelImage->GetBufferedRegion();

    // assume that boundary buffer 20 pixels
    for (int i = 0; i < DIMENSIONS; i++) {
        idx1[i] = std::max(_labelBounds[imageId][labelId][2*i] - 10, wholeRegion.GetIndex(i)) ;
        idx2[i] = std::min(_labelBounds[imageId][labelId][2*i+1] + 10, wholeRegion.GetUpperIndex()[i]);
    }
    extractRegion.SetIndex(idx1);
    extractRegion.SetUpperIndex(idx2);

    LabelExtractFilterType::Pointer extractionFilter = LabelExtractFilterType::New();
    extractionFilter->SetInput(labelImage);
//    extractionFilter->SetExtractionRegion(extractRegion);
    extractionFilter->SetExtractionRegion(wholeRegion);
    extractionFilter->Update();
    LabelImage::Pointer outputImage = extractionFilter->GetOutput();

    LabelImageIteratorType iter(outputImage, outputImage->GetBufferedRegion());
    iter.GoToBegin();
    while (!iter.IsAtEnd()) {
        if (iter.Get() != labelId) {
            iter.Set(0);
        }
        ++iter;
    }

    ImageIO<LabelImage> io;
    char fname[128];
    sprintf(fname, "/tmpfs/label_%02d_%02d.nii.gz", imageId, labelId);
    io.WriteImage(fname, outputImage);
    return outputImage;
}

void ParticleTrainer::saveParticles() {
    const int n = _system->GetNumberOfSubjects();

    for (int i = 0; i < n; i++) {
        ParticleSubject& subj = _system->operator[](i);
        char fileOut[128];
        sprintf(fileOut, "/tmpfs/particles_%02d.txt", i);
        ofstream of(fileOut);
        subj.WriteParticles(of);
        of.close();
    }
}

void ParticleTrainer::loadParticles() {
    int nsubjs = 0;
    for (int i = 0; i < 10; i++) {
        char fileIn[128];
        sprintf(fileIn, "/tmpfs/particles_%02d.txt", i);
        ifstream fin(fileIn);
        if (!fin.good()) {
            break;
        }
        nsubjs++;
    }

    cout << "# of particle subjects: " << nsubjs << endl;
    _system->InitializeSystem(nsubjs, 0);
    for (int i = 0; i < nsubjs; i++) {
        char fileIn[128];
        sprintf(fileIn, "/tmpfs/particles_%02d.txt", i);
        ifstream fin(fileIn);

        _system->operator[](i).ReadParticles(fin, -1);
        _system->operator[](i).WriteParticles(cout);
    }
}


void ParticleTrainer::computeNearestNeighbors(NeighborWeightMatrix& neighbors, int subjId, int labelId) {
    ParticleSubject& subj = _system->operator[](subjId);
    int npoints = subj.GetNumberOfPoints();

    const int kNearest = 3;
    neighbors.clear();

    ParticleWeightVector tempNeighbors;

    for (int i = 0; i < npoints; i++) {

        Particle& pi = subj[i];
        if (pi.label != labelId) {
            continue;
        }
        tempNeighbors.clear();
        for (int j = 0; j < npoints; j++) {
            Particle& pj = subj[j];
            if (pj.label != labelId) {
                continue;
            }
            if (i == j) {
                continue;
            } else {
                double w = pi.Dist2(pj);
                double a = -1.0 / w;
                tempNeighbors.push_back(ParticleWeightPair(j, a));
            }
        }

        std::sort(tempNeighbors.begin(), tempNeighbors.end(), ParticleWeightComparator());

        neighbors.push_back(NeighborWeightPair(i, ParticleWeightVector()));
        for (int j = 0; j < kNearest; j++) {
            neighbors[i].second.push_back(tempNeighbors[j]);
        }
    }
}


void ParticleTrainer::initialClosestCorrespondences() {
    int nSubjs = _system->GetNumberOfSubjects();
    int nPoints = _system->operator[](0).GetNumberOfPoints();

    ParticleSubject& subj0 = _system->operator[](0);
    for (int s = 1; s < nSubjs; s++) {
        // construct laplacian matrix
        ParticleSubject& subjs = _system->operator[](s);
        for (int i = 0; i < nPoints; i++) {
            Particle& p0 = subj0[i];
            p0.correspondence = i;
            double distMin = 1e9;
            int correspondenceId = -1;
            for (int j = 0; j < nPoints; j++) {
                Particle& ps = subjs[j];
                if (ps.correspondence >= 0 || ps.label != p0.label) {
                    continue;
                }
                double dist2 = p0.Dist2(ps);
                if (dist2 < distMin) {
                    distMin = dist2;
                    correspondenceId = j;
                }
            }
            subjs[correspondenceId].correspondence = i;
            subjs[correspondenceId].correspondenceScore = distMin;
//            cout << "Minimum Distance:" << distMin << endl;
        }

        subjs.SortByCorrespondence();
    }


    const bool printCorrespondence = false;
    if (printCorrespondence) {
        for (int i = 0; i < nPoints; i++) {
            for (int s = 0; s < nSubjs; s++) {
                cout << _system->operator[](s)[i].correspondence << " ";
                fordim (k) {
                    cout << _system->operator[](s)[i].x[k] << " ";
                }
                cout << "; ";
            }
            cout << endl;
        }
    }
}

void ParticleTrainer::establishCorrespondences() {

    int nsubjs = _system->GetNumberOfSubjects();
    int npoints = _system->operator[](0).GetNumberOfPoints();
    int nPointsPerLabel = 0;
    
    ParticleSubject& subj0 = _system->operator[](0);
    for (int i = 0; i < npoints; i++) {
        if (subj0[i].label == 1) {
            nPointsPerLabel ++;
        }
    }

    VNLDoubleMatrix laplacian(nPointsPerLabel * 2, nPointsPerLabel * 2);
    laplacian.fill(0);

    ParticleWeightVector neighbors;
    const int kNeighbors = 15;

    // kNearest method
    {
        /*

    // l-index for Laplacian
    int m = 0;
    for (int i = 0; i < npoints; i++) {
        Particle& pi = subj0[i];
        if (pi.label != 1) {
            continue;
        }
//        double sumW = 0;
        neighbors.clear();
        int n = 0;
        for (int j = 0; j < npoints; j++) {
            Particle& pj = subj0[j];
            if (pj.label != 1) {
                continue;
            }
            if (i == j) {
                laplacian[m][n] = 0;
            } else {
                double w = pi.Dist2(pj);
                double a = -1.0 / w;
                neighbors.push_back(ParticleWeightPair(n, a));
//                sumW += a;
            }
            n++;
        }

        std::sort(neighbors.begin(), neighbors.end(), ParticleWeightComparator());
        double sum = 0;
        for (int j = 0; j < kNeighbors; j++) {
            laplacian[m][neighbors[j].first] = neighbors[j].second;
            sum += neighbors[j].second;
        }
        laplacian[m][m] = -sum;

        m++;
    }
         */
    }


    // laplacian for first point sets
    NeighborWeightMatrix& pointSet0 = _neighborCollection[0];
    for (int i = 0; i < pointSet0.size(); i++) {
        int srcId = pointSet0[i].first;
        ParticleWeightVector& neighborPoints = pointSet0[i].second;
        double sum = 0;
        for (int j = 0; j < neighborPoints.size(); j++) {
            ParticleWeightPair& neighborPair = neighborPoints[j];
            int dstId = neighborPair.first;
            laplacian[srcId][dstId] = neighborPair.second;
            sum += neighborPair.second;
        }
        laplacian[srcId][srcId] = -sum;
    }


    for (int s = 1; s < nsubjs; s++) {
        ParticleSubject& subj = _system->operator[](s);

        NeighborWeightMatrix& pointSet = _neighborCollection[s];
        for (int i = 0; i < pointSet.size(); i++) {
            int srcId = pointSet[i].first + npoints;
            ParticleWeightVector& neighborPoints = pointSet[i].second;
            double sum = 0;
            for (int j = 0; j < neighborPoints.size(); j++) {
                ParticleWeightPair& neighborPair = neighborPoints[j];
                int dstId = neighborPair.first + npoints;
                laplacian[srcId][dstId] = neighborPair.second;
                sum += neighborPair.second;
            }
            laplacian[srcId][srcId] = -sum;
        }

        /*
        // construct laplacian matrix
        ParticleSubject& subj = _system->operator[](s);

        int nPointsPerSubj = subj.GetNumberOfPoints();

        // l-index for Laplacian
        int m = nPointsPerLabel * s;
        for (int i = 0; i < nPointsPerSubj; i++) {
            Particle& pi = subj[i];
            if (pi.label != 1) {
                continue;
            }

            // neighbor loop
            neighbors.clear();
            int n = nPointsPerLabel * s;
            for (int j = 0; j < nPointsPerSubj; j++) {
                Particle& pj = subj[j];
                if (pj.label != 1) {
                    continue;
                }
                if (i == j) {
                    laplacian[m][n] = 0;
                } else {
                    double w = pi.Dist2(pj);
                    double a = -1.0 / w;
                    neighbors.push_back(ParticleWeightPair(n, a));
                }
                n++;
            }

            std::sort(neighbors.begin(), neighbors.end(), ParticleWeightComparator());
            double sum = 0;
            for (int j = 0; j < kNeighbors; j++) {
                laplacian[m][neighbors[j].first] = neighbors[j].second;
                sum += neighbors[j].second;
            }
            laplacian[m][m] = -sum;
            m++;
        }
         */

        {
            int m = 0;
            // check to connect two graphs
            bool useDualConnection = true;
            if (useDualConnection) {
                m = 0;
                for (int i = 0; i < npoints; i += 5) {
                    Particle p0i = subj0[i];
                    Particle psi = subj[i];

                    if (p0i.label == 1 && psi.label == 1) {
                        double dsi = p0i.Dist2(psi);
                        laplacian[m][m + nPointsPerLabel] = -1.0 / dsi;
                        laplacian[m + nPointsPerLabel][m] = -1.0 / dsi;
                        laplacian[m][m] += 1/dsi;
                        laplacian[m + nPointsPerLabel][m + nPointsPerLabel] += 1/dsi;
                        m++;
                    }
                }
            }
        }

        cout << "==================" << endl;
        cout << laplacian << endl;
        cout << "==================" << endl;

        vnl_symmetric_eigensystem<double> eig(laplacian);
        cout << "Eigenvalues: " << eig.D << endl;
        VNLDoubleVector v1 = eig.get_eigenvector(1);
        VNLDoubleVector v2 = eig.get_eigenvector(2);
        VNLDoubleVector v3 = eig.get_eigenvector(3);
        VNLDoubleVector v4 = eig.get_eigenvector(4);

        const int kDegree = 8;

        int l = 0;
        for (int i = 0; i < npoints; i++) {
            Particle& pi = subj0[i];
            if (pi.label != 1) {
                continue;
            }

//            pi.correspondenceScore = v1[l]*v1[l]/eig.D[1] + v2[1]*v2[1]/eig.D[2] + v3[1]*v3[1]/eig.D[3];
            double score = 0;
            for (int k = 1; k <= kDegree; k++) {
                score += ((eig.V[l][k] * eig.V[l][k]) / eig.D[k]);
            }
            pi.correspondenceScore = score;
            l++;
        }

        for (int i = 0; i < npoints; i++) {
            Particle& pi = subj[i];
            if (pi.label != 1) {
                continue;
            }

//            pi.correspondenceScore = v1[l]*v1[l]/eig.D[1] + v2[1]*v2[1]/eig.D[2] + v3[1]*v3[1]/eig.D[3];
            double score = 0;
            for (int k = 1; k <= kDegree; k++) {
                score += ((eig.V[l][k] * eig.V[l][k]) / eig.D[k]);
            }
            pi.correspondenceScore = score;
            l++;
        }

        vnl_vector<int> marker(npoints);
        marker.fill(0);

        for (int i = 0; i < npoints; i++) {
        	Particle& pi = subj0[i];
        	double distMin = 1e9;
        	int correspondenceId = -1;
        	for (int j = 0; j < npoints; j++) {
        		Particle& pj = subj[j];
        		if (marker[j] > 0) {
        			continue;
        		}
//        		double dist = (pi.correspondenceScore - pj.correspondenceScore);
//        		double dist = (v1[i] - v1[j+npoints])*(v1[i] - v1[j+npoints])/eig.D[1] + (v2[i] - v2[j+npoints])*(v2[i] - v2[j+npoints])/eig.D[2];



                double dist2 = 0;
                for (int k = 1; k <= kDegree; k++) {
                    double vki = eig.V[i][k];
                    double vkj = eig.V[j+npoints][k];
                    dist2 += (vki - vkj)*(vki - vkj)/eig.D[k];
                }

                double dist = dist2;

        		if (dist < distMin) {
        			correspondenceId = j;
        			distMin = dist;
        		}
        	}
            if (correspondenceId < 0) {
                cout << i << " has no correspondence" << endl;
            } else {
                subj[correspondenceId].correspondence = i;
                marker[correspondenceId] = 1;
            }
        }

        subj.SortByCorrespondence();
    }
}


void ParticleTrainer::removeOutliers() {
    int nsubjs = _system->GetNumberOfSubjects();
    int npoints = _system->operator[](0).GetNumberOfPoints();
    int nPointsPerLabel = 0;


    // variable to mark up outliers
    vnl_vector<int> outliers(npoints);
    outliers.fill(false);

    // check number of points per label
    ParticleSubject& subj0 = _system->operator[](0);
    for (int i = 0; i < npoints; i++) {
        if (subj0[i].label == 1) {
            nPointsPerLabel ++;
        }
        subj0[i].outlier = false;
    }

    // variables to store distance between corresponding points
    VNLVector dist(npoints);
    VNLVector dist2(npoints), distz(npoints);

    // loop over subjects and compare with the template
    for (int i = 1; i < nsubjs; i++) {
        ParticleSubject& subji = _system->operator[](i);
        dist.fill(0);
        for (int j = 0; j < npoints; j++) {
            double d = 0;
            fordim (k) {
                d += (subj0[j].x[k]-subji[j].x[k])*(subj0[j].x[k]-subji[j].x[k]);
            }
            dist[j] = sqrt(d);
        }
        double mean = dist.mean();
        for (int j = 0; j < npoints; j++) {
            dist2[j] = dist[j] * dist[j];
        }
        double x2mean = dist2.mean();
        double stdev = sqrt(x2mean - mean * mean);

        for (int j = 0; j < npoints; j++) {
            distz[j] = (dist[j] - mean) / stdev;
            if (std::abs(distz[j]) > 0.5) {
                // outlier
                outliers[j] = true;
                subj0[j].outlier = true;
                subji[j].enabled = false;
                subj0[j].enabled = false;
            }
        }
    }

}

void ParticleTrainer::setNumberOfPointSets(int n) {
    _neighborCollection.resize(n);
}

NeighborWeightMatrix& ParticleTrainer::getPointSets(int n) {
    return _neighborCollection[n];
}

void ParticleTrainer::computeDelaunayTriangulation(int subjId, int labelId, int pointSetId) {
    NeighborWeightMatrix& matrix = _neighborCollection[pointSetId];

    ParticleSubject& subj = _system->operator[](subjId);
    int npoints = subj.GetNumberOfPoints();

    vtkPoints* points = vtkPoints::New();
    points->SetNumberOfPoints(npoints);
    for (int i = 0; i < npoints; i++) {
        points->SetPoint(i, subj[i].x);
    }

    vtkPolyData* pointsData = vtkPolyData::New();
    pointsData->SetPoints(points);

    vtkDelaunay2D* triFilter = vtkDelaunay2D::New();
    triFilter->SetInput(pointsData);
    triFilter->Update();

    vtkPolyData* triangles = triFilter->GetOutput();
    vtkIdList* cellIds = vtkIdList::New();
    for (int i = 0; i < npoints; i++) {
        Particle& pi = subj[i];
        if (pi.label != labelId) {
            continue;
        }
        triangles->GetPointCells(i, cellIds);

        ParticleWeightVector connects;
        std::set<int> neighborIds;
        for (int j = 0; j < cellIds->GetNumberOfIds(); j++) {
            int id = cellIds->GetId(j);
            vtkCell* cell = triangles->GetCell(id);
            if (cell->GetNumberOfEdges() <= 0) {
                continue;
            }
            for (int k = 0; k < cell->GetNumberOfEdges(); k++) {
                vtkCell* edge = cell->GetEdge(k);
                vtkIdList* pointIdList = edge->GetPointIds();
                if (pointIdList->GetId(0) == i || pointIdList->GetId(1) == i) {
                    continue;
                }
                neighborIds.insert(pointIdList->GetId(0));
                neighborIds.insert(pointIdList->GetId(1));
            }
        }

        for (std::set<int>::iterator ii = neighborIds.begin(); ii != neighborIds.end(); ii++) {
            int nId = *ii;
            Particle& pn = subj[nId];
            double dist = pi.Dist2(pn);
            connects.push_back(ParticleWeightPair(nId, -1.0/dist));
        }

        matrix.push_back(NeighborWeightPair(i, connects));
    }
    triFilter->Delete();
}
