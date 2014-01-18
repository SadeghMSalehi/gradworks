//
//  ParticleWarp.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/31/13.
//
//

#include "piParticleWarp.h"
#include <itkDisplacementFieldTransform.h>
#include <itkBSplineScatteredDataPointSetToImageFilter.h>
#include <itkBSplineTransform.h>
#include <itkWarpImageFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkPoints.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkPolyDataWriter.h>
#include <vtkTriangle.h>

namespace pi {

    const int __SplineOrder = 3;

    typedef itk::BSplineScatteredDataPointSetToImageFilter
    <DisplacementFieldPointSetType, DisplacementFieldType> BSplineFilterType;
    typedef BSplineFilterType::WeightsContainerType WeightsContainerType;
    typedef itk::BSplineTransform<PointReal,__Dim, __SplineOrder> BSplineTransform;
    typedef itk::WarpImageFilter<RealImage, RealImage, DisplacementFieldType> WarpImageFilterType;
    typedef itk::WarpImageFilter<LabelImage, LabelImage, DisplacementFieldType> WarpLabelFilterType;

    void ParticleWarp::setParameters(pi::ConfigFile &config) {
        controlSpacing = config["particles.bspline-transform.control-point-spacing"];
    }

    void ParticleWarp::estimateBsplineWarp(Px::Vector &src, Px::Vector &dst) {
        if (reference.IsNull()) {
            cout << "Cannot estimate grid size without reference image!" << endl;
            return;
        }

        // how itk::PointSet manage the number of points?
        if (m_FieldPoints.IsNull()) {
            m_FieldPoints = DisplacementFieldPointSetType::New();
            m_FieldPoints->Initialize();
        }

        int n = src.size();
        for (int i = 0; i < n; i++) {
            IntPointSetType::PointType srcPoint;
            IntPointSetType::PointType dstPoint;
            fordim(j) {
                srcPoint[j] = src[i][j];
                dstPoint[j] = dst[i][j];
            }
            VectorType vector;
            fordim(j) {
                vector[j] = dst[i][j] - src[i][j];
            }

            // if there is no warp, check below
//            cout << srcPoint << " => " << dstPoint << " : " << vector << endl;


            m_FieldPoints->SetPoint(i, srcPoint);
            m_FieldPoints->SetPointData(i, vector);

        }

        LabelImage::SizeType imageSize = reference->GetBufferedRegion().GetSize();
        LabelImage::SpacingType imageSpacing = reference->GetSpacing();
        LabelImage::PointType imageOrigin = reference->GetOrigin();
        LabelImage::DirectionType imageDirection = reference->GetDirection();


        BSplineFilterType::Pointer bspliner = BSplineFilterType::New();
        BSplineFilterType::ArrayType numControlPoints;

        // control point is given by its number, not in physical unit!
        int numOfLevels = 3;
        for (int i = 0; i < imageSize.GetSizeDimension(); i++) {
            numControlPoints[i] = imageSize[i] / controlSpacing + __SplineOrder;
        }

        try {
            // debug: reparameterized point component is outside
            bspliner->SetOrigin(imageOrigin);
            bspliner->SetSpacing(imageSpacing);
            bspliner->SetSize(imageSize);
            bspliner->SetDirection(imageDirection);
            bspliner->SetGenerateOutputImage(true);
            bspliner->SetNumberOfLevels(numOfLevels);
            bspliner->SetSplineOrder(__SplineOrder);
            bspliner->SetNumberOfControlPoints(numControlPoints);
            bspliner->SetInput(m_FieldPoints);

            /*
            if (m_UseWeights && n == m_Weights.size()) {
                BSplineFilterType::WeightsContainerType::Pointer weights = BSplineFilterType::WeightsContainerType::New();

                weights->Reserve(n);

                for (int i = 0; i < n; i++) {
                    // how to choose src or dst weight?
                    weights->SetElement(i, m_Weights[i]);
                }
                bspliner->SetPointWeights(weights.GetPointer());
            }
             */

            bspliner->Update();
            displacementField = bspliner->GetOutput();
        } catch (itk::ExceptionObject& e) {
            e.Print(std::cout);
        }
    }


    LabelImage::Pointer ParticleWarp::warpLabel(LabelImage::Pointer input) {
        if (displacementField.IsNull()) {
            return input;
        }
        WarpLabelFilterType::Pointer warpFilter = WarpLabelFilterType::New();
        warpFilter->SetInput(input);
        warpFilter->SetDisplacementField(displacementField);
        warpFilter->SetOutputParametersFromImage(input);
        warpFilter->Update();
        return warpFilter->GetOutput();
    }

    RealImage::Pointer ParticleWarp::warpImage(RealImage::Pointer input) {
        if (displacementField.IsNull()) {
            return input;
        }
        WarpImageFilterType::Pointer warpFilter = WarpImageFilterType::New();
        warpFilter->SetInput(input);
        warpFilter->SetDisplacementField(displacementField);
        warpFilter->SetOutputParametersFromImage(input);
        warpFilter->Update();
        return warpFilter->GetOutput();
    }

#pragma mark ParticleMesh Implementations
    void ParticleMesh::constructNeighbors(int regionId, int nPx, PxSubj& subj, double cutoff, PxGlobal::Neighbors& neighbors) {

        // total number of particles (>= nPx)
        const int tnPx = subj.size();

        // idmap from vtkPoint to subject-particle
        // because vtkPoint will have new id to run delaunay filter
        IntVector idmap(tnPx);
        std::fill(idmap.begin(), idmap.end(), -1);



        // setup vtkPoints and store vtkPoint-ParticleId mapping
        vtkPoints* points = vtkPoints::New();
        points->SetNumberOfPoints(nPx);
        for (int i = 0, t = 0; i < tnPx; i++) {
            if (subj.attrs[i].label == regionId) {
                double p[3] = { 0, };
                fordim (k) {
                    p[k] = subj.particles[i][k];
                }
                points->SetPoint(t, p);
                idmap[t] = i;
                t++;
            }
        }

        // run vtk delaunay filter
        vtkPolyData* pointsData = vtkPolyData::New();
        pointsData->SetPoints(points);


        // assume delaunay2d filter will give correct results regardless of z-coordinate
        vtkDelaunay2D* triFilter = vtkDelaunay2D::New();
        triFilter->SetInput(pointsData);
        triFilter->Update();

        // create triangle filter to retrieve edges
        vtkPolyData* triangles = triFilter->GetOutput();

//        vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
//        writer->SetInput(triangles);
//        writer->SetFileName("/tmpfs/mesh.vtk");
//        writer->Write();
//        writer->Delete();

//        triangles->Print(cout);

        // loop over points and identify neighbor points
        typedef std::set<int> IntSet;
        typedef std::vector<IntSet> IntSetVector;

        IntSetVector nbrIds;
        nbrIds.resize(tnPx);

        for (int j = 0; j < triangles->GetNumberOfCells(); j++) {
            vtkCell* cells = triangles->GetCell(j);
            vtkIdList* points = cells->GetPointIds();
            for (int k = 0; k < points->GetNumberOfIds(); k++) {
                int pxk = idmap[ points->GetId(k) ];
                for (int l = k + 1; l < points->GetNumberOfIds(); l++) {
                    int pxl = idmap [ points->GetId(l) ];
                    double dkl = subj.particles[pxk].dist2(subj.particles[pxl]);
                    if (dkl <= cutoff * cutoff) {
                        nbrIds[pxk].insert(pxl);
                        nbrIds[pxl].insert(pxk);
                    }
                }
            }
        }

        // iterate over IntSetVector and copy to neighbors
        // some degree of duplication is required due to avoid 'set' structure
        // is 'set' also efficiently iteratable?
        // question for std implementation
        IntSetVector::const_iterator iter = nbrIds.begin();
        for (int i = 0; iter != nbrIds.end(); iter++, i++) {
            neighbors[i].clear();
            IntSet::const_iterator niter = iter->begin();
            for (; niter != iter->end(); niter++) {
                neighbors[i].push_back(*niter);
            }
        }

        /*
         vtkIdList* cellIds = vtkIdList::New();
         for (int i = 0; i < nPx; i++) {
            triangles->GetPointCells(i, cellIds);

            // store in set to avoid duplicate ids
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

            // loop over found neighbors
            for (std::set<int>::iterator ii = neighborIds.begin(); ii != neighborIds.end(); ii++) {
                int nId = *ii;
                int pxId = idmap[nId];
                neighbors[i].push_back(pxId);
            }
            triFilter->Delete();
        }
         */
        triFilter->Delete();
    }
}