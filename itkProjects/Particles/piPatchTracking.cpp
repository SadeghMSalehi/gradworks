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

    void PatchTracking::setPatchRegion(RealImage::RegionType region) {
        _patchRegion[0] = region;
        _patchRegion[1] = region;
        extractPatch(0);
        extractPatch(1);
        transformPatchRegion();
    }

    void PatchTracking::translatePatchRegion(RealImage::RegionType region) {
        _patchRegion[0] = region;
        transformPatchRegion();
        extractPatch(0, false);
        resamplePatch();
    }

    RealImage::Pointer PatchTracking::getPatch(int i) {
        return _patches[i];
    }

    QPolygonF& PatchTracking::getPatchPolygon(int i) {
        return _patchPolygon[i];
    }

    void PatchTracking::resamplePatch() {
        typedef itk::ResampleImageFilter<RealImage,RealImage> ResampleFilter;
        ResampleFilter::Pointer resampler = ResampleFilter::New();
        resampler->SetInput(_images[1]);
        if (_transform.IsNotNull()) {
            resampler->SetTransform(dynamic_cast<ResampleFilter::TransformType*>(_transform.GetPointer()));
        }
        resampler->SetReferenceImage(_patches[0]);
        resampler->UseReferenceImageOn();
//        resampler->GraftOutput(_patches[1]);
        resampler->Update();
        _patches[1] = resampler->GetOutput();
    }

    void PatchTracking::extractPatch(int i, bool newPatch) {
        typedef itk::ExtractImageFilter<RealImage, RealImage> ExtractFilter;
        ExtractFilter::Pointer filter = ExtractFilter::New();
        filter->SetInput(_images[i]);
        filter->SetExtractionRegion(_patchRegion[i]);
        try {
            filter->Update();
        } catch (itk::ExceptionObject& ex) {
            ex.Print(cout);
        }
        _patches[i] = filter->GetOutput();
        _patches[i]->DisconnectPipeline();
    }

    // set up optimization routine
    void PatchTracking::beginTracking() {
        _transform = TransformType::New();

        ParametersType initialParams;
        initialParams.SetSize(_transform->GetNumberOfParameters());
        initialParams.Fill(0);

        _costFunc = CostFunctionType::New();
        _costFunc->SetFixedImage(_patches[0]);
        _costFunc->SetFixedImageRegion(_patchRegion[0]);
        _costFunc->SetMovingImage(_images[1]);
        _costFunc->UseAllPixelsOn();

        _costFunc->SetTransform(dynamic_cast<CostFunctionType::TransformType*>(_transform.GetPointer()));
        _costFunc->SetInterpolator(LinearImageInterpolatorType::New());
        _costFunc->Initialize();


        OptimizerType::ScalesType scales;
        scales.SetSize(initialParams.Size());
        scales.Fill(1);

        OptimizerProgress<OptimizerType>::Pointer progress = OptimizerProgress<OptimizerType>::New();

        _optimizer = OptimizerType::New();
        _optimizer->SetCostFunction(_costFunc);
        _optimizer->SetInitialPosition(initialParams);
        _optimizer->SetScales(scales);
//        _optimizer->AddObserver(itk::IterationEvent(), progress);

        setupOptimizer(_optimizer);
        
        _optimizer->StartOptimization();

        cout << "Iterations: " << _optimizer->GetCurrentIteration() << endl;
        cout << "Stop Reason: " << _optimizer->GetStopConditionDescription() << endl;

        const ParametersType& params = _optimizer->GetCurrentPosition();
        _transform->SetParameters(params);

        resamplePatch();
        transformPatchRegion();
    }

    void PatchTracking::transformPatchRegion() {
        // index to physical point
        RealImage::RegionType& r = _patchRegion[0];
        RealImage::IndexType i[4];
        i[3] = i[1] = i[0] = r.GetIndex();
        i[2][1] = (i[1][1] += r.GetSize(1));
        i[2][0] = (i[3][0] += r.GetSize(0));

        QPolygonF transformedRegion;
        RealImage::PointType p[4];
        for (int j = 0; j < 4; j++) {
            if (_transform.IsNotNull()) {
                _images[0]->TransformIndexToPhysicalPoint(i[j], p[j]);
                p[j] = _transform->TransformPoint(p[j]);
                _images[0]->TransformPhysicalPointToIndex(p[j], i[j]);
            }
            transformedRegion << QPointF(i[j][0], i[j][1]);
        }
        transformedRegion << QPointF(i[0][0], i[0][1]);
        _patchPolygon[1] = transformedRegion;
    }

    void PatchTracking::stepTracking() {

    }

    void PatchTracking::setupOptimizer(RegularStepGradientDescentOptimizer* opti) {
        RealImage::SpacingType spacing = _images[0]->GetSpacing();

        opti->SetMinimumStepLength(spacing[0]*0.1);
        opti->SetMaximumStepLength(spacing[0]);
        opti->SetNumberOfIterations(100);
        opti->SetRelaxationFactor(0.5);
    }

    void PatchTracking::setupOptimizer(FRPROptimizer* opti) {
        RealImage::SpacingType spacing = _images[0]->GetSpacing();

        opti->SetStepLength(spacing[0]*0.1);
    }
}