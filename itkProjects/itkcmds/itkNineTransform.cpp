#include "itkImageIO.h"
#include "itkCompositeTransform.h"
#include "itkTranslationTransform.h"
#include "itkEuler3DTransform.h"
#include "itkScaleVersor3DTransform.h"
#include "itkTransformFileWriter.h"
#include "itkMyMetric.h"
#include "itkMyFRPROptimizer.h"
#include "itkMath.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkOptimizationReporter.h"
#include "itkContinuousIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkResampleImageFilter.h"
#include "iostream"
#include "vector"

using namespace std;

#define angle2rad(a) (a*itk::Math::pi/180.0)

class NineReg {
public:
    typedef itk::Image<unsigned int,3> LabelType;
    typedef itk::Image<float,3> ImageType;
    typedef itk::CompositeTransform<double,3> CompositeTransform;
    typedef itk::TranslationTransform<double> InitialTransform;
    typedef itk::ScaleVersor3DTransform<double> NineTransform;
    typedef itk::TransformFileWriter TransformWriter;
    typedef itk::MyMetric<ImageType, ImageType> Metric;
    typedef itk::MyFRPROptimizer Optimizer;
    typedef itk::LinearInterpolateImageFunction<ImageType> Interpolator;
    typedef itk::NearestNeighborInterpolateImageFunction<ImageType> InterpolatorNN;

    typedef OptimizationReporter<Optimizer> OptiReporter;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;


    ImageType::Pointer _src;
    ImageType::Pointer _dst;
    LabelType::Pointer _dstLabel;

    ImageType::PointType _srcCenter;
    ImageType::PointType _dstCenter;

    InitialTransform::Pointer _initialTransform;
    NineTransform::Pointer _finalTransform;
    
    InitialTransform::ParametersType _initialResult;
    NineTransform::ParametersType _finalResult;
    Metric::FixedImageIndexContainer _labelIndexes;
    InitialTransform::ParametersType _centerOfRotation;
    char* _transformOut;
    char* _resampledOut;

    void LoadFiles(int argc, char* argv[]) {
        itkcmds::itkImageIO<ImageType> io;
        itkcmds::itkImageIO<LabelType> io2;
        _src = io.ReadImageT(argv[1]);
        _dst = io.ReadImageT(argv[2]);
        _dstLabel = io2.ReadImageT(argv[3]);
        _transformOut = argv[4];
        _resampledOut = argv[5];

        cout << "Transform Output: " << _transformOut << endl;
        cout << "Resampled Output: " << _resampledOut << endl;

        ImageType::SizeType szDst = _dst->GetBufferedRegion().GetSize();
        itk::ContinuousIndex<double,3> szIdx;
        for (int i = 0; i < 3; i++) {
            szIdx[i] = szDst[i] / 2.0;
        }
        _dst->TransformContinuousIndexToPhysicalPoint(szIdx, _dstCenter);

        itk::ImageRegionConstIteratorWithIndex<LabelType> labelIter(_dstLabel, _dstLabel->GetBufferedRegion());
        for (labelIter.GoToBegin(); !labelIter.IsAtEnd(); ++labelIter) {
            LabelType::PixelType label = labelIter.Get();
            if (label > 0) {
                _labelIndexes.push_back(labelIter.GetIndex());
            }
        }

        _centerOfRotation.SetSize(ImageType::ImageDimension);
        for (int i = 0; i < 3; i++) {
            _centerOfRotation[i] = _dstCenter[i];
        }
    }

    void RunRigid() {
        _initialTransform = InitialTransform::New();

        _initialTransform->SetFixedParameters(_centerOfRotation);

        OptiReporter::Pointer optiReporter = OptiReporter::New();
        Metric::Pointer metric = Metric::New();
        metric->SetFixedImage(_dst);
        metric->SetMovingImage(_src);
        metric->SetInterpolator(Interpolator::New());
        metric->SetTransform(_initialTransform);
        //metric->SetFixedImageRegion(_dst->GetBufferedRegion());
        metric->SetFixedImageIndexes(_labelIndexes);
//        metric->SetUseAllPixels(true);
        metric->Initialize();

        Optimizer::Pointer opti = Optimizer::New();
        opti->SetCostFunction(metric);
        Optimizer::ScalesType scales;
        scales.SetSize(InitialTransform::ParametersDimension);
        scales.Fill(1);
        //scales[0] = scales[1] = scales[2] = .1;
        
        opti->SetUseUnitLengthGradient(true);
        opti->SetScales(scales);
        opti->SetInitialPosition(_initialTransform->GetParameters());
        opti->AddObserver(itk::StartEvent(), optiReporter);
        opti->AddObserver(itk::IterationEvent(), optiReporter);
        opti->StartOptimization();
        cout << "Current Cost: " << opti->GetCurrentCost() << endl;
        _initialResult = opti->GetCurrentPosition();
        _initialTransform->SetParameters(opti->GetCurrentPosition());
    }

    void RunNine() {
        _finalTransform = NineTransform::New();
//        NineTransform::OutputVectorType initialTranslation;
//        initialTranslation[0] = _initialResult[0];
//        initialTranslation[1] = _initialResult[1];
//        initialTranslation[2] = _initialResult[2];

        // _finalTransform->SetTranslation(initialTranslation);
        _finalTransform->SetFixedParameters(_centerOfRotation);
        
        OptiReporter::Pointer optiReporter = OptiReporter::New();
        Metric::Pointer metric = Metric::New();
        metric->SetFixedImage(_dst);
        metric->SetMovingImage(_src);
        metric->SetInterpolator(InterpolatorNN::New());
        metric->SetTransform(_finalTransform);
        metric->SetFixedImageIndexes(_labelIndexes);
//        metric->SetUseAllPixels(true);
        metric->Initialize();

        Optimizer::Pointer opti = Optimizer::New();
        opti->SetCostFunction(metric);
        Optimizer::ScalesType scales;
        scales.SetSize(NineTransform::ParametersDimension);
        scales.Fill(1);
        scales[0] = scales[1] = scales[2] = 10;
        scales[6] = scales[7] = scales[8] = 100;

        opti->SetMaximumIteration(500);
        opti->SetUseUnitLengthGradient(true);
        opti->SetStepLength(0.5);
        opti->SetScales(scales);
        opti->SetInitialPosition(_finalTransform->GetParameters());
        opti->AddObserver(itk::StartEvent(), optiReporter);
        opti->AddObserver(itk::IterationEvent(), optiReporter);
        opti->StartOptimization();
        cout << "Current Cost: " << opti->GetCurrentCost() << endl;
        _finalResult = opti->GetCurrentPosition();
        _finalTransform->SetParameters(opti->GetCurrentPosition());
    }

    void SaveTransform() {
        TransformWriter::Pointer writer = TransformWriter::New();
        //writer->AddTransform(_initialTransform);
        writer->AddTransform(_finalTransform);
        writer->SetFileName(_transformOut);
        writer->Update();


    }

    void ResampleResult(const char* name) {
        InterpolatorNN::Pointer resampleInterpolator = InterpolatorNN::New();
        ResamplerType::Pointer resample = ResamplerType::New();
        resample->SetTransform(_finalTransform);
        resample->SetInput(_src);
        resample->UseReferenceImageOn();
        resample->SetReferenceImage(_dst);
        resample->SetInterpolator(resampleInterpolator);
        resample->Update();
        ImageType::Pointer resampledImage = resample->GetOutput();
        itkcmds::itkImageIO<ImageType> io;
        io.WriteImageT(name, resampledImage);
    }


    int main(int argc, char* argv[]) {
        LoadFiles(argc, argv);
        // RunRigid();
        RunNine();
        SaveTransform();
        ResampleResult(_resampledOut);
        return 0;
    }
};


int main(int argc, char* argv[]) {
    NineReg nrg;
    nrg.main(argc, argv);
}