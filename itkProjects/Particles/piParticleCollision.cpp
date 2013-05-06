//
//  piParticleCollision.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/11/13.
//
//

#include "piParticleCollision.h"
#include "itkZeroCrossingImageFilter.h"
#include "itkInvertIntensityImageFilter.h"
#include "piImageIO.h"
#include "piImageProcessing.h"
#include "boost/algorithm/string/predicate.hpp"

#if DIMENSIONS == 3
#define IsEqual(x, y) (x[0]==y[0]&&x[1]==y[1]&&x[2]==y[2])
#else
#define IsEqual(x, y) (x[0]==y[0]&&x[1]==y[1])
#endif

#define x2string1(x) (x[0]+1)<<","<<(x[1]+1)<<","<<(x[2]+1)

namespace pi {
    bool ParticleCollision::ComputeContactPoint(DataReal *p0, DataReal *p1, ContactPoint& cp) {
        // do binary search between x0 and x1;
        IntIndex x0, x1, xm;
        RealIndex y0, y1, ym;

        fordim (k) {
            // round up
            x0[k] = (int) (p0[k] + 0.5);
            x1[k] = (int) (p1[k] + 0.5);
            y0[k] = p0[k];
            y1[k] = p1[k];
        }
        const bool inx0 = IsRegionInside(x0);
        const bool inx1 = IsRegionInside(x1);
        cp.status = NO_CROSSING;
        if ((inx0 && inx1) || (!inx0 && !inx1)) {
            cp.status = NO_CROSSING;
            return false;
        }
        bool isx0In = inx0 && !inx1;
        if (IsCrossing(x0)) {
            fordim(k) {
                cp.cp[k] = y0[k];
            }
            cp.status = STARTING_CONTACT;
            return true;
        } else if (IsCrossing(x1)) {
            fordim(k) {
                cp.cp[k] = y1[k];
            }
            cp.status = ENDING_CONTACT;
            return true;
        }
        cp.status = CROSSING_CONTACT;
        for (int cnt = 0; true; cnt++) {
            fordim (k) {
                ym[k] = (y0[k] + y1[k]) / 2.0;
                xm[k] = ym[k] + 0.5;
                x0[k] = y0[k] + 0.5;
                x1[k] = y1[k] + 0.5;
            }
            if (IsEqual(xm, x0) || IsEqual(xm, x1)) {
                // found crossing
                fordim (k) {
                    cp.cp[k] = y0[k];
                }
                return true;
            }
            const bool isCrossing = IsCrossing(xm);
            if (isCrossing) {
                // found crossing
                fordim (k) {
                    cp.cp[k] = ym[k];
                }
                return true;
            } else {
                // set up for next loop
                const bool isInside = IsRegionInside(xm);
                if (isInside == isx0In) {
                    fordim (k) {
                        y0[k] = ym[k];
                    }
                } else {
                    fordim (k) {
                        y1[k] = ym[k];
                    }
                }
            }
        }
        return true;
    }

    bool ParticleCollision::ComputeNormal(DataReal *contactPoint, DataReal* normalOutput) {
        RealImage::PointType point;
        fordim (k) {
            point[k] = contactPoint[k];
        }
        GradientPixel normal = m_NormalPicker->Evaluate(point);
        fordim (k) {
            normalOutput[k] = -normal[k];
        }

        return true;
    }

    DataReal ParticleCollision::ComputeDistance(DataReal *x1, DataReal *cp) {
        DataReal dist = 0;
        fordim (k) {
            dist = (x1[k] - cp[k]) * (x1[k] - cp[k]);
        }
        return sqrt(dist);
    }

    void ParticleCollision::ComputeClosestBoundary(Particle& pi, DataReal *x1, DataReal *contactPoint) {
        RealImage::PointType point;
        fordim (k) {
            point[k] = x1[k];
        }

        // compute offset as interpolated value
        // wondering if this is not going to cause any problem
        VectorType offset = m_DistOffsetPicker->Evaluate(point);
        fordim (k) {
            contactPoint[k] = x1[k] + offset[k] * m_ImageSpacing[k];
        }

#ifndef BATCH
        if (::abs(offset[0]) > 10 || ::abs(offset[1]) > 10 || ::abs(offset[2]) > 10)  {
            cout << "too large projection offset = [" << offset[0] << "," << offset[1] << "," << offset[2] << "] at [" << x1[0] << "," << x1[1] << "," << x1[2] << "]" << "; particle subj = " << pi.subj << " id = " << pi.idx << endl;
        }
        fordim(k) {
            point[k] = contactPoint[k];
        }
        if (!m_DistOffsetPicker->IsInsideBuffer(point)) {
            cout << "fail to find closest boundary" << endl;
            // rollback
            fordim (k) {
                contactPoint[k] = x1[k];
            }
        };
#endif
    }

    void ParticleCollision::SetLabelImage(LabelImage::Pointer labelImage) {
        m_LabelImage = labelImage;
    }

    void ParticleCollision::SetBinaryMask(LabelImage::Pointer labelImage) {
        m_BinaryMask = labelImage;
    }

    void ParticleCollision::UpdateImages() {
        // check if binary mask should be updated
        if (m_BinaryMask.IsNotNull() && m_LabelImage.IsNotNull()
                && (m_BinaryMask->GetMTime() > m_LabelImage->GetMTime())) {
            return;
        }

        ImageProcessing proc;
        ImageIO<LabelImage> io;

        if (!LoadBinaryMask(binaryMaskCache)) {
            cout << "binary mask creation error!!" << endl;
            return;
        };
        if (!LoadDistanceMap(distanceMapCache)) {
            cout << "distance map creation error!!" << endl;
            return;
        }

        typedef itk::InvertIntensityImageFilter<LabelImage,LabelImage> InvertFilterType;
        InvertFilterType::Pointer invertFilter = InvertFilterType::New();
        invertFilter->SetInput(m_BinaryMask);
        invertFilter->SetMaximum(255);
        invertFilter->Update();
        LabelImage::Pointer invertedMask = invertFilter->GetOutput();

        typedef itk::ZeroCrossingImageFilter<LabelImage,LabelImage> FilterType;
        FilterType::Pointer filter = FilterType::New();
        filter->SetInput(invertedMask);
        filter->SetForegroundValue(1);
        filter->SetBackgroundValue(0);
        filter->Update();
        m_ZeroCrossing = filter->GetOutput();

        m_Gradient = proc.ComputeGaussianGradient(m_BinaryMask);

        m_CrossingPicker = NNLabelInterpolatorType::New();
        m_CrossingPicker->SetInputImage(m_ZeroCrossing);

        m_RegionPicker = NNLabelInterpolatorType::New();
        m_RegionPicker->SetInputImage(m_BinaryMask);

        m_NormalPicker = GradientInterpolatorType::New();
        m_NormalPicker->SetInputImage(m_Gradient);

        m_DistOffsetPicker = NNVectorImageInterpolatorType::New();
        m_DistOffsetPicker->SetInputImage(m_DistanceMap);

    }


    bool ParticleCollision::LoadBinaryMask(std::string maskCache) {
        if (m_BinaryMask.IsNotNull()) {
            return true;
        }
        ImageIO<LabelImage> io;
        if (io.FileExists(maskCache.c_str())) {
            // use cache regardless of using mask smoothing
            m_BinaryMask = io.ReadImage(maskCache.c_str());
            return true;
        }
        // to create binary mask, the label image is mandatory
        if (m_LabelImage.IsNull()) {
            return false;
        }
        ImageProcessing proc;
        m_BinaryMask = proc.ThresholdToBinary(m_LabelImage);
        if (applyMaskSmoothing) {
            m_BinaryMask = proc.SmoothLabelMap(m_BinaryMask);
        }
        
        // store cache
        if (maskCache != "") {
            io.WriteImage(maskCache.c_str(), m_BinaryMask);
        }
        return m_BinaryMask.IsNotNull();
    }

    bool ParticleCollision::LoadDistanceMap(string filename) {
        ImageIO<VectorImage> io;
        if (m_DistanceMap.IsNotNull()) {
            // if already exist, then do nothing
            return false;
        }
        if (io.FileExists(filename.c_str())) {
            m_DistanceMap = io.ReadImage(filename.c_str());
        } else {
            ImageProcessing proc;
            m_DistanceMap = proc.ComputeDistanceMap(m_BinaryMask);
            if (filename != "") {
                io.WriteImage(filename.c_str(), m_DistanceMap, !boost::algorithm::ends_with(filename, ".nrrd"));
            }
        }
        return true;
    }

    void ParticleCollision::ConstrainPoint(ParticleSubject& subj) {
        m_ImageSpacing = subj.GetLabel()->GetSpacing();

        const int nPoints = subj.GetNumberOfPoints();
        for (int i = 0; i < nPoints; i++) {
            Particle &p = subj.m_Particles[i];
            VNLVector normal(__Dim, 0);

            LabelImage::IndexType idx;
            subj.ComputeIndexX(p, idx);

            const bool isValidRegion = IsRegionInside(idx);
            const bool isContacting = IsCrossing(idx);

            if (isValidRegion && !isContacting) {
                p.collisionEvent = false;
                continue;
            }

            DataReal contactPoint[__Dim];

            // initialize contact point as current point
            forset (p.x, contactPoint);
            if (!isValidRegion) {
                // important!!
                // this function will move current out-of-region point into the just closest boundary
                // Point compuation is done in physical point
                ComputeClosestBoundary(p, p.x, contactPoint);
            }
            forset (contactPoint, p.x);
            p.collisionEvent = true;
        }
    }

    void ParticleCollision::ProjectForceAndVelocity(pi::ParticleSubject &subj) {
        const int nPoints = subj.GetNumberOfPoints();
        VNLVector normal(__Dim);
        for (int i = 0; i < nPoints; i++) {
            Particle& p = subj[i];
            normal.fill(0);
            if (p.collisionEvent) {
                // project velocity and forces
                ComputeNormal(p.x, normal.data_block());
                normal.normalize();
                DataReal nv = dimdot(p.v, normal);
                DataReal nf = dimdot(p.f, normal);
                fordim (k) {
                    p.v[k] = (p.v[k] - nv * normal[k]);
                    p.f[k] -= nf * normal[k];
                }
            }
        }
    }

    void ParticleCollision::HandleCollision(ParticleSubject& subj) {
        const int nPoints = subj.GetNumberOfPoints();
        for (int i = 0; i < nPoints; i++) {
            Particle &p = subj.m_Particles[i];
            VNLVector normal(__Dim, 0);


            // compute the index of the particle
            LabelImage::IndexType idx;
            subj.ComputeIndexX(p, idx);

            const bool isValidRegion = IsRegionInside(idx);
            const bool isContacting = IsCrossing(idx);
            
            if (isValidRegion && !isContacting) {
                continue;
            }

            DataReal contactPoint[__Dim];
            
            // initialize contact point as current point
            forset (p.x, contactPoint);
            
            if (!isValidRegion) {
                // important!!
                // this function will move current out-of-region point into the closest boundary
                ComputeClosestBoundary(p, p.x, contactPoint);
            }

            // project velocity and forces
            ComputeNormal(contactPoint, normal.data_block());
            normal.normalize();
            DataReal nv = dimdot(p.v, normal);
            DataReal nf = dimdot(p.f, normal);
            fordim (k) {
                p.v[k] = (p.v[k] - nv * normal[k]);
                p.f[k] -= nf * normal[k];
            }
        }
    }
    
    void ParticleCollision::HandleCollision(ParticleSubjectArray& subjs) {
        const int nShapes = subjs.size();

        for (int n = 0; n < nShapes; n++) {
            HandleCollision(subjs[n]);
        }
    }
}