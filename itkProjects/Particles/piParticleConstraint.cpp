//
//  ParticleConstraint.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/23/12.
//
//

#include "vnlCommon.h"
#include "piImageDef.h"
#include "piParticleConstraint.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"


namespace pi {
    typedef itk::BinaryThresholdImageFilter<LabelImage, LabelImage> BinaryThresholdFilterType;
    typedef itk::VectorMagnitudeImageFilter<VectorImage, LabelImage> GradientMagnitudeFilterType;

    template <class TIn, class TOut>
    class InvertLabel {
    public:
        InvertLabel() {}
        ~InvertLabel() {}
        bool operator!=(const InvertLabel&) const { return false; }
        bool operator==(const InvertLabel& other) const { return !(*this != other); }
        inline TOut operator()(const TIn& A) const { return !(A>0); }
    };

    template <class T1, class T2>
    class BinaryThreshold {
        public :
        BinaryThreshold() {}
        ~BinaryThreshold() {}
        bool operator!=(const BinaryThreshold&) const { return false; }
        bool operator==(const BinaryThreshold& other) const { return !(*this != other); }
        inline T2 operator()(const T1& A) const { return (A>0)?1:0; }
    };

    void ParticleConstraint::Clear() {
        m_DistanceMaps.clear();
        m_DistanceMapInterpolators.clear();
        m_InsideDistanceVectorMaps.clear();
        m_OutsideDistanceVectorMaps.clear();
        m_GradientInterpolators.clear();
        m_GradientMaps.clear();
    }

    void ParticleConstraint::SetImageList(LabelVector& imageList) {
        int nSubj = imageList.size();
        Clear();
        for (int i = 0; i < nSubj; i++) {
            LabelImage::Pointer labelMap = imageList[i];

            // create binary image for a mask for a correct distance map
            BinaryThresholdFilterType::Pointer binThreshFilter = BinaryThresholdFilterType::New();
            binThreshFilter->SetInput(labelMap);
            binThreshFilter->SetInsideValue(1);
            binThreshFilter->SetOutsideValue(0);
            binThreshFilter->SetLowerThreshold(1);
            binThreshFilter->SetUpperThreshold(255);
            LabelImage::Pointer binaryMap = binThreshFilter->GetOutput();

//            itkcmds::itkImageIO<LabelImage> io;
//            char buf[128];
//            sprintf(buf, "/tmpfs/binary_%02d.nrrd", i);
//            io.WriteImageT(buf, binaryMap);

            // construct signed distance filter
            SignedDistanceMapFilterType::Pointer distmapFilter = SignedDistanceMapFilterType::New();
            distmapFilter->SetInput(binaryMap);
            distmapFilter->Update();
            m_DistanceMaps.push_back(distmapFilter->GetOutput());
            m_OutsideDistanceVectorMaps.push_back(distmapFilter->GetVectorDistanceMap());

            //m_OutsideDistanceVectorMaps[0].Print(cout);
            // to compute inside offset correctly, invert the label map
            typedef itk::UnaryFunctorImageFilter<LabelImage, LabelImage, InvertLabel<LabelImage::PixelType, LabelImage::PixelType> > InvertImageFilterType;
            InvertImageFilterType::Pointer invertFilter = InvertImageFilterType::New();
            invertFilter->SetInput(binaryMap);
            invertFilter->Update();

            // compute inside offset using inverted image
            DistanceMapFilterType::Pointer insideDistmapFilter = DistanceMapFilterType::New();
            insideDistmapFilter->SetInput(invertFilter->GetOutput());
            insideDistmapFilter->Update();
            m_InsideDistanceVectorMaps.push_back(insideDistmapFilter->GetVectorDistanceMap());

            LinearImageInterpolatorType::Pointer interpol = LinearImageInterpolatorType::New();
            interpol->SetInputImage(m_DistanceMaps.back());
            m_DistanceMapInterpolators.push_back(interpol);

            // compute gradient for outside boundary

            LabelImageGradientFilterType::Pointer gradient = LabelImageGradientFilterType::New();
            gradient->SetInput(binaryMap);
            gradient->SetSigma(.5);
            try {
                gradient->Update();
            } catch (itk::ExceptionObject& e) {
                cout << e.GetDescription() << endl;
                throw;
            }
            m_GradientMaps.push_back(gradient->GetOutput());


            GradientInterpolatorType::Pointer gradientInterpolator = GradientInterpolatorType::New();
            gradientInterpolator->SetInputImage(gradient->GetOutput());
            m_GradientInterpolators.push_back(gradientInterpolator);

            GradientMagnitudeFilterType::Pointer magnitudeFilter = GradientMagnitudeFilterType::New();
            magnitudeFilter->SetInput(gradient->GetOutput());
            magnitudeFilter->Update();
        }
    }

    double ParticleConstraint::GetDistance(int subjId, LinearImageInterpolatorType::ContinuousIndexType &idx) const {
        return m_DistanceMapInterpolators[subjId]->EvaluateAtContinuousIndex(idx);
    }

    bool ParticleConstraint::IsInsideRegion(int subjId, LinearImageInterpolatorType::ContinuousIndexType& idx) const {
        if (subjId < m_DistanceMaps.size()) {
            return m_DistanceMapInterpolators[subjId]->IsInsideBuffer(subjId);
        }
        return false;
    }

    bool ParticleConstraint::IsInsideRegion(int subjId, LinearImageInterpolatorType::IndexType& idx) const {
        if (subjId < m_DistanceMaps.size()) {
            return m_DistanceMapInterpolators[subjId]->IsInsideBuffer(subjId);
        }
        return false;
    }

    ParticleConstraint::DistanceVectorImageType::PixelType ParticleConstraint::GetInsideOffset(int subjId, DoubleImage::IndexType& idx) const {
        return m_InsideDistanceVectorMaps[subjId]->GetPixel(idx);
    }

    ParticleConstraint::DistanceVectorImageType::PixelType ParticleConstraint::GetOutsideOffset(int subjId, DoubleImage::IndexType& idx) const {
        return m_OutsideDistanceVectorMaps[subjId]->GetPixel(idx);
    }

    ParticleConstraint::GradientPixelType ParticleConstraint::GetGradient(int subjId, GradientInterpolatorType::ContinuousIndexType& idx) const {
        //    return m_GradientMaps[subjId]->GetPixel(idx);
        return m_GradientInterpolators[subjId]->EvaluateAtContinuousIndex(idx);
    }

    void ParticleConstraint::ApplyConstraint(ParticleSubjectArray& shapes) {        
        // boundary constraint
        const int nSubj = shapes.size();
        const int nPoints = shapes[0].m_Particles.size();

        // iterate over all subjects
        for (int n = 0; n < nSubj; n++) {
            // iterate over all particles
            for (int i = 0; i < nPoints; i++) {
                Particle& pi = shapes[n][i];
                LinearImageInterpolatorType::ContinuousIndexType nidx;
                LinearImageInterpolatorType::IndexType idx;
                
                fordim(k) {
                    nidx[k] = idx[k] = pi.x[k];
                }

                // case when a particle escaped from the region
                if (!IsInsideRegion(n, nidx)) {
                    // actually, this shouldn't happen
                    fordim(k) {
                        pi.v[k] = 0;
                    }
                    continue;
                }

                // otherwise, compute distance from the boundary
                double dist = GetDistance(n, nidx);
                if (dist >= 0) {
                    OffsetType offset = GetOutsideOffset(n, idx);
                    fordim(k) {
                        pi.x[k] += offset[k];
                        // update current position to boundary
                        idx[k] = nidx[k] = pi.x[k];
                    }
                }

                double m_COR = 0.1;
                VectorType g = GetGradient(n, nidx);
                VNLVectorRef normal(__Dim, g.GetDataPointer());
                VNLVectorRef iVel(__Dim, pi.v);
                VNLVectorRef iForce(__Dim, pi.f);

                double normalMagnitude = normal.two_norm();
                if (normalMagnitude > 0.4) {
                    normal.normalize();

                    // velocity should be zero toward normal direction
                    double normalSpeed = dot_product(iVel, normal);
                    if (normalSpeed < 0) {
                        VNLVector newVelocity = iVel - 2 * normalSpeed * normal;
                        newVelocity *= m_COR;
                        // how to know current timestep?
                        //                        newVelocity /= 0.1;
                        iForce.fill(0);
                        iVel.copy_in(newVelocity.data_block());
                    }

                    // remove normal term of the current force
                    double normalForce = dot_product(normal, iForce);
                    if (normalForce < 0) {
                        VNLVector newForce = iForce + normalForce * normal;
                        iForce.copy_in(newForce.data_block());
                    }
                }
            }
        }
    }
}