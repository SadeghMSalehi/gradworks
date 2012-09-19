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
#include "string"

using namespace std;

#define angle2rad(a) (a*itk::Math::pi/180.0)

template<class TransformType>
class RegistrationEngine {
public:
    typedef itk::Image<unsigned int,3> LabelType;
    typedef itk::Image<float,3> ImageType;
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

    typename TransformType::Pointer _transform;
    typename TransformType::ParametersType _transformResult;
    typename TransformType::ParametersType _centerOfRotation;

    Metric::FixedImageIndexContainer _labelIndexes;

    char* _transformOut;
    char* _resampledOut;
    string _method;

    void LoadFiles(int argc, char* argv[]) {
        itkcmds::itkImageIO<ImageType> io;
        itkcmds::itkImageIO<LabelType> io2;
        _src = io.ReadImageT(argv[2]);
        _dst = io.ReadImageT(argv[3]);
        _dstLabel = io2.ReadImageT(argv[4]);
        _transformOut = argv[5];
        _resampledOut = argv[6];

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

    void RunRegistration() {
        _transform = TransformType::New();
        float params[12] = { 0.994639, -0.00207486, -0.00192788, -0.000832746, 1.00641, -0.00146783, 0.000495955, -0.000724458, 0.999367, 0.989471, 10.2177, 4.1862 };

        typename TransformType::ParametersType initialParams;
        initialParams.SetSize(TransformType::ParametersDimension);
        for (int i = 0; i < TransformType::ParametersDimension; i++) {
            initialParams[i] = params[i];
        }
        _transform->SetParameters(initialParams);
        _transform->SetFixedParameters(_centerOfRotation);

        OptiReporter::Pointer optiReporter = OptiReporter::New();
        Metric::Pointer metric = Metric::New();
        metric->SetFixedImage(_dst);
        bool useIndexes = true;
        if (useIndexes) {
            metric->SetFixedImageIndexes(_labelIndexes);
        } else {
            metric->SetFixedImageRegion(_dst->GetBufferedRegion());
            metric->SetUseAllPixels(true);
        }

        metric->SetMovingImage(_src);
        metric->SetInterpolator(InterpolatorNN::New());
        metric->SetTransform(_transform);
        metric->Initialize();

        Optimizer::Pointer opti = Optimizer::New();
        opti->SetCostFunction(metric);
        Optimizer::ScalesType scales;
        scales.SetSize(TransformType::ParametersDimension);
        scales.Fill(1);
        if (_method == "affine") {
            cout << "apply affine scaling ..." << endl;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    if (i == j) {
                        scales[3*i+j] = 30;
                    } else {
                        scales[3*i+j] = 160;
                    }
                }
            }
            scales[9] = scales[10] = scales[11] = 0.1;
        } else if (_method == "nine") {
            scales[0] = scales[1] = scales[2] = 10;
            scales[6] = scales[7] = scales[8] = 100;
        }

        opti->SetMaximumIteration(1000);
        opti->SetMaximumLineIteration(10);
        opti->SetUseUnitLengthGradient(true);
        opti->SetStepLength(0.25);
        opti->SetScales(scales);
        opti->SetToFletchReeves();
        opti->SetInitialPosition(_transform->GetParameters());
        opti->AddObserver(itk::StartEvent(), optiReporter);
        opti->AddObserver(itk::IterationEvent(), optiReporter);
        opti->StartOptimization();
        cout << "Current Cost: " << opti->GetCurrentCost() << endl;
        _transformResult = opti->GetCurrentPosition();
        _transform->SetParameters(opti->GetCurrentPosition());
    }

    void SaveTransform() {
        TransformWriter::Pointer writer = TransformWriter::New();
        //writer->AddTransform(_initialTransform);
        writer->AddTransform(_transform);
        writer->SetFileName(_transformOut);
        writer->Update();


    }

    void ResampleResult(const char* name) {
        InterpolatorNN::Pointer resampleInterpolator = InterpolatorNN::New();
        ResamplerType::Pointer resample = ResamplerType::New();
        resample->SetTransform(_transform);
        resample->SetInput(_src);
        resample->UseReferenceImageOn();
        resample->SetReferenceImage(_dst);
        resample->SetInterpolator(resampleInterpolator);
        resample->Update();
        ImageType::Pointer resampledImage = resample->GetOutput();
        itkcmds::itkImageIO<ImageType> io;
        io.WriteImageT(name, resampledImage);
    }


    int main(int argc, char* argv[], string method) {
        _method = method;
        LoadFiles(argc, argv);
        RunRegistration();
        SaveTransform();
        ResampleResult(_resampledOut);
        return 0;
    }
};


int main(int argc, char* argv[]) {
    string method(argv[1]);
    if (method == "affine") {
        RegistrationEngine<itk::AffineTransform<double,3> > reg;
        reg.main(argc, argv, method);
    } else if (method == "nine") {
        RegistrationEngine<itk::ScaleVersor3DTransform<double> > reg;
        reg.main(argc, argv, method);
    }
}