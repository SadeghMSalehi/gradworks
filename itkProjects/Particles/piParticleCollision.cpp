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

    bool ParticleCollision::ComputeNormal(DataReal *cp, DataReal* ret) {
        RealIndex idx;
        fordim (k) {
            idx[k] = cp[k];
        }
        GradientPixel normal = m_NormalPicker->EvaluateAtContinuousIndex(idx);
        fordim (k) {
            ret[k] = -normal[k];
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

    void ParticleCollision::ComputeClosestBoundary(DataReal *x1, DataReal * cp) {
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

    void ParticleCollision::UseBinaryMaskSmoothing() {
        m_ApplySmoothing = true;
    }

    void ParticleCollision::UseBinaryMaskSmoothingCache(const char *cacheName) {
        m_BinaryMaskSmoothingCacheName = cacheName;
    }

    void ParticleCollision::UseDistanceMapCache(const char *cacheName) {
        m_DistanceMapCacheName = cacheName;
    }

    void ParticleCollision::UpdateImages() {
        ImageProcessing proc;
        itkcmds::itkImageIO<LabelImage> io;
        if (m_ApplySmoothing) {
            bool generateNewMask = true;
            if (m_BinaryMaskSmoothingCacheName != "") {
                if (io.FileExists(m_BinaryMaskSmoothingCacheName.c_str())) {
                    m_BinaryMask = io.ReadImageT(m_BinaryMaskSmoothingCacheName.c_str());
                    generateNewMask = false;
                }
            }
            if (generateNewMask) {
                m_BinaryMask = proc.SmoothLabelMap(m_BinaryMask);
                if (m_BinaryMaskSmoothingCacheName != "") {
                    io.WriteImageT(m_BinaryMaskSmoothingCacheName.c_str(), m_BinaryMask);
                }
            }
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

        m_Gradient = proc.ComputeGaussianGradient(m_BinaryMask);

        m_CrossingPicker = NNLabelInterpolatorType::New();
        m_CrossingPicker->SetInputImage(m_ZeroCrossing);

        m_RegionPicker = NNLabelInterpolatorType::New();
        m_RegionPicker->SetInputImage(m_BinaryMask);

        m_NormalPicker = GradientInterpolatorType::New();
        m_NormalPicker->SetInputImage(m_Gradient);

        if (m_DistanceMapCacheName != "") {
            LoadDistanceMap(m_DistanceMapCacheName.c_str());
        }

        if (m_DistanceMap.IsNull()) {
            m_DistanceMap = proc.DistanceMap(m_BinaryMask);
            m_DistOffsetPicker = NNVectorImageInterpolatorType::New();
            m_DistOffsetPicker->SetInputImage(m_DistanceMap);
            if (m_DistanceMapCacheName != "") {
                SaveDistanceMap(m_DistanceMapCacheName.c_str());
            }
        }
    }


    bool ParticleCollision::LoadBinaryMask(std::string file) {
        itkcmds::itkImageIO<LabelImage> io;
        if (io.FileExists(file.c_str())) {
            m_BinaryMask = io.ReadImageT(file.c_str());
            return true;
        }
        return false;
    }

    bool ParticleCollision::LoadDistanceMap(const char* filename) {
        itkcmds::itkImageIO<VectorImage> io;
        if (m_DistanceMap.IsNull() && io.FileExists(filename)) {
            m_DistanceMap = io.ReadImageT(filename);
            if (m_DistanceMap.IsNotNull()) {
                m_DistOffsetPicker = NNVectorImageInterpolatorType::New();
                m_DistOffsetPicker->SetInputImage(m_DistanceMap);
                return true;
            }
        }
        return false;
    }

    void ParticleCollision::SaveGradientMagnitude(const char* filename) {
        itkcmds::itkImageIO<DoubleImage> io;
        if (m_Gradient.IsNotNull()) {
            ImageProcessing proc;
            DoubleImage::Pointer mag = proc.ComputeMagnitudeMap(m_Gradient);
            io.WriteImageT(filename, mag);
        }
    }


    void ParticleCollision::SaveDistanceMap(const char* filename) {
        itkcmds::itkImageIO<VectorImage> io;
        if (m_DistanceMap.IsNotNull()) {
            io.WriteImageT(filename, m_DistanceMap);
        }
    }
    
    
    void ParticleCollision::Write(std::string b, std::string c) {
        itkcmds::itkImageIO<LabelImage> io;
        io.WriteImageT("/tmpfs/binary.nrrd", m_BinaryMask);
//        io.WriteImageT(b.c_str(), m_InvertedBinaryMask);
        io.WriteImageT("/tmpfs/edge.nrrd", m_ZeroCrossing);
    }

    void ParticleCollision::HandleCollision(ParticleSubject& subj) {
        const int nPoints = subj.GetNumberOfPoints();
        for (int i = 0; i < nPoints; i++) {
            Particle &p = subj.m_Particles[i];
            VNLVector normal(__Dim, 0);
            LabelImage::IndexType idx;
            fordim (k) {
                idx[k] = p.x[k];
            }
            
            const bool isValidRegion = IsRegionInside(idx);
            const bool isContacting = IsCrossing(idx);
            
            if (isValidRegion && !isContacting) {
                continue;
            }
            
            if (!isValidRegion) {
                // important!!
                // this function will move current out-of-region point into the closest boundary
                ComputeClosestBoundary(p.x, p.x);
            }
            
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
    
    void ParticleCollision::HandleCollision(ParticleSubjectArray& subjs) {
        const int nShapes = subjs.size();

        for (int n = 0; n < nShapes; n++) {
            HandleCollision(subjs[n]);
        }
    }
}