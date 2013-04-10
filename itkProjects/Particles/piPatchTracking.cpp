//
//  piPatchTracking.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/9/13.
//
//

#include "piPatchTracking.h"
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkExtractImageFilter.h>
#include <itkResampleImageFilter.h>

using namespace itk;
using namespace std;

namespace pi {

    template <class T>
    class OptimizerProgress: public itk::Command {
    public:
        /** Standard class typedefs. */
        typedef OptimizerProgress Self;
        typedef itk::Command Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        itkTypeMacro(OptimizerProgress, itk::Command);
        itkNewMacro(OptimizerProgress);

        virtual void Execute(itk::Object *caller, const itk::EventObject & event) {
            this->Execute((const itk::Object*) caller, event);
        }

        virtual void Execute(const itk::Object *caller,
                             const itk::EventObject & event) {
            const T* realCaller = dynamic_cast<const T*>(caller);
            if (realCaller == NULL) {
                return;
            }
            cout << realCaller->GetCurrentIteration() << ": " << realCaller->GetCurrentPosition() << endl;
        }

    protected:
        OptimizerProgress() {}

        virtual ~OptimizerProgress() {}

    private:
        OptimizerProgress(const Self &);        //purposely not implemented
        void operator=(const Self &); //purposely not implemented
    };

    PatchTracking::PatchTracking() {
        
    }

    void PatchTracking::setImage(int i, RealImage::Pointer image) {
        _images[i] = image;
    }

    void PatchTracking::setInitialRegion(int i, RealImage::RegionType region) {
        _initialRegion[i] = region;
        itk::ContinuousIndex<double,2> idx;
        for (int j = 0; j < 2; j++) {
            idx[j] = (region.GetIndex(j) + region.GetSize(j) / 2.0);
        }
//        _images[i]->TransformContinuousIndexToPhysicalPoint(idx, _centerOfRegion[i]);

        extractPatch(i);
    }

    RealImage::Pointer PatchTracking::getPatch(int i) {
        return _patches[i];
    }

    RealImage::IndexType& PatchTracking::getFinalIndex(int i) {
        return _finalIndex[i];
    }

    void PatchTracking::extractPatch(int i) {
        typedef itk::ExtractImageFilter<RealImage, RealImage> ExtractFilter;
        ExtractFilter::Pointer filter = ExtractFilter::New();
        filter->SetInput(_images[i]);
        filter->SetExtractionRegion(_initialRegion[i]);
        filter->Update();
        _patches[i] = filter->GetOutput();
        _patches[i]->DisconnectPipeline();
    }

    // set up optimization routine
    void PatchTracking::beginTracking() {
        ParametersType initialParams;
        initialParams.SetSize(2);
        initialParams.Fill(0);

        _transform = TransformType::New();

        _costFunc = CostFunctionType::New();
        _costFunc->SetFixedImage(_patches[0]);
        _costFunc->SetFixedImageRegion(_initialRegion[0]);
        _costFunc->SetMovingImage(_images[1]);
        _costFunc->UseAllPixelsOn();

        _costFunc->SetTransform(dynamic_cast<CostFunctionType::TransformType*>(_transform.GetPointer()));
        _costFunc->SetInterpolator(LinearImageInterpolatorType::New());
        _costFunc->Initialize();


        OptimizerType::ScalesType scales;
        scales.SetSize(2);
        scales.Fill(1);

        OptimizerProgress<OptimizerType>::Pointer progress = OptimizerProgress<OptimizerType>::New();

        RealImage::SpacingType spacing = _images[0]->GetSpacing();
        _optimizer = OptimizerType::New();
        _optimizer->SetCostFunction(_costFunc);
        _optimizer->SetInitialPosition(initialParams);
        _optimizer->SetScales(scales);
        _optimizer->SetMinimumStepLength(spacing[0]*0.1);
        _optimizer->SetMaximumStepLength(spacing[0]);
        _optimizer->SetNumberOfIterations(100);
        _optimizer->SetRelaxationFactor(0.5);
        _optimizer->AddObserver(itk::IterationEvent(), progress);
        _optimizer->StartOptimization();

        cout << "Iterations: " << _optimizer->GetCurrentIteration() << endl;
        cout << "Stop Reason: " << _optimizer->GetStopConditionDescription() << endl;

        const ParametersType& params = _optimizer->GetCurrentPosition();
        _transform->SetParameters(params);
        
        typedef itk::ResampleImageFilter<RealImage,RealImage> ResampleFilter;
        ResampleFilter::Pointer resampler = ResampleFilter::New();
        resampler->SetInput(_images[1]);
        resampler->SetTransform(dynamic_cast<ResampleFilter::TransformType*>(_transform.GetPointer()));
        resampler->SetReferenceImage(_patches[0]);
        resampler->UseReferenceImageOn();
        resampler->GraftOutput(_patches[1]);
        resampler->Update();

        TransformType::InputPointType point[2];
        for (int j = 0; j < 2; j++) {
            point[0][j] = _initialRegion[0].GetIndex(j);
            point[1][j] = _initialRegion[0].GetIndex(j) + _initialRegion[0].GetSize(j);
        }
        
        for (int i = 0; i < 2; i++) {
            _finalPoints[i] = _transform->TransformPoint(point[i]);
            _images[0]->TransformPhysicalPointToIndex(_finalPoints[i], _finalIndex[i]);
        }
    }

    void PatchTracking::stepTracking() {

    }
}