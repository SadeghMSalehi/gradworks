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
#include "itkImageIO.h"
#include "piImageProcessing.h"

#ifdef DIMENSION3
#define IsEqual(x, y) (x[0]==y[0]&&x[1]==y[1]&&x[2]==y[2])
#else
#define IsEqual(x, y) (x[0]==y[0]&&x[1]==y[1])
#endif

#define x2string1(x) (x[0]+1)<<","<<(x[1]+1)<<","<<(x[2]+1)

namespace pi {
    bool ParticleCollision::ComputeContactPoint(double *p0, double *p1, ContactPoint& cp) {
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

    bool ParticleCollision::ComputeNormal(double *cp, double* ret) {
        RealIndex idx;
        fordim (k) {
            idx[k] = cp[k];
        }
        VectorType normal = m_NormalPicker->EvaluateAtContinuousIndex(idx);
        fordim (k) {
            ret[k] = -normal[k];
        }

        return true;
    }

    double ParticleCollision::ComputeDistance(double *x1, double *cp) {
        double dist = 0;
        fordim (k) {
            dist = (x1[k] - cp[k]) * (x1[k] - cp[k]);
        }
        return sqrt(dist);
    }

    void ParticleCollision::ComputeClosestBoundary(double *x1, double * cp) {
        IntIndex idx;
        fordim (k) {
            idx[k] = x1[k];
        }
        VectorType offset = m_DistOffsetPicker->EvaluateAtIndex(idx);
        fordim (k) {
            cp[k] = x1[k] + offset[k];
        }
    }

    void ParticleCollision::SetBinaryMask(LabelImage::Pointer binary) {
        m_BinaryMask = binary;
    }

    void ParticleCollision::UpdateImages() {
        ImageProcessing proc;
        if (m_ApplySmoothing) {
            m_BinaryMask = proc.SmoothLabelMap(m_BinaryMask);
        }
        
        typedef itk::InvertIntensityImageFilter<LabelImage,LabelImage> InvertFilterType;
        InvertFilterType::Pointer invertFilter = InvertFilterType::New();
        invertFilter->SetInput(m_BinaryMask);
        invertFilter->SetMaximum(255);
        invertFilter->Update();
        m_InvertedBinaryMask = invertFilter->GetOutput();

        typedef itk::ZeroCrossingImageFilter<LabelImage,LabelImage> FilterType;
        FilterType::Pointer filter = FilterType::New();
        filter->SetInput(m_InvertedBinaryMask);
        filter->SetForegroundValue(1);
        filter->SetBackgroundValue(0);
        filter->Update();
        m_ZeroCrossing = filter->GetOutput();

        m_Gradient = proc.ComputeNormal(m_BinaryMask);

        m_CrossingPicker = NNLabelInterpolatorType::New();
        m_CrossingPicker->SetInputImage(m_ZeroCrossing);

        m_RegionPicker = NNLabelInterpolatorType::New();
        m_RegionPicker->SetInputImage(m_BinaryMask);

        m_NormalPicker = LinearVectorImageInterpolatorType::New();
        m_NormalPicker->SetInputImage(m_Gradient);

        m_DistanceMap = proc.DistanceMap(m_BinaryMask);
        m_DistOffsetPicker = NNVectorImageInterpolatorType::New();
        m_DistOffsetPicker->SetInputImage(m_DistanceMap);
    }

    void ParticleCollision::Write(std::string b, std::string c) {
        itkcmds::itkImageIO<LabelImage> io;
        io.WriteImageT("/tmpfs/binary.nrrd", m_BinaryMask);
//        io.WriteImageT(b.c_str(), m_InvertedBinaryMask);
        io.WriteImageT("/tmpfs/edge.nrrd", m_ZeroCrossing);
    }

}