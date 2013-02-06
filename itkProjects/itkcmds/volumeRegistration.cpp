#include "itkImageIO.h"
#include "itkCompositeTransform.h"
#include "itkTranslationTransform.h"
#include "itkEuler3DTransform.h"
#include "itkScaleVersor3DTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkTransformFileWriter.h"
#include "itkMyMetric.h"
#include "itkMath.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkOptimizationReporter.h"
#include "itkContinuousIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkResampleImageFilter.h"
#include "itkMyScaleVersor3DTransformOptimizer.h"
#include "itkMyFRPROptimizer.h"

#include "iostream"
#include "vector"
#include "string"

using namespace std;

#define angle2rad(a) (a*itk::Math::pi/180.0)

#define USE_GD_OPTIMIZER


template<class TransformType>
class RegistrationEngine {
public:
    typedef itk::Image<unsigned int,3> LabelType;
    typedef itk::Image<float,3> ImageType;
    typedef itk::TransformFileWriter TransformWriter;
    typedef itk::MyMetric<ImageType, ImageType> Metric;
    // typedef itk::MyFRPROptimizer Optimizer;
    typedef itk::MyScaleVersor3DTransformOptimizer Optimizer;
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
    }

    void RunRegistration() {
        _transform = TransformType::New();
        OptiReporter::Pointer optiReporter = OptiReporter::New();
        Metric::Pointer metric = Metric::New();
        metric->SetFixedImage(_dst);
        bool useIndexes = false;
        if (useIndexes) {
            _centerOfRotation.SetSize(ImageType::ImageDimension);
            for (int i = 0; i < ImageType::ImageDimension; i++) {
                _centerOfRotation[i] = i;
            }
            itk::ImageRegionConstIteratorWithIndex<LabelType> labelIter(_dstLabel, _dstLabel->GetBufferedRegion());
            int nPixels = 0;
            for (labelIter.GoToBegin(); !labelIter.IsAtEnd(); ++labelIter) {
                LabelType::PixelType label = labelIter.Get();
                if (label > 0) {
                    _labelIndexes.push_back(labelIter.GetIndex());
                    for (int i = 0; i < ImageType::ImageDimension; i++) {
                        _centerOfRotation[i] += labelIter.GetIndex()[i];
                    }
                    nPixels ++;
                }
            }
            for (int i = 0; i < ImageType::ImageDimension; i++) {
                _centerOfRotation[i] /= nPixels;
            }
            metric->SetFixedImageIndexes(_labelIndexes);
            _transform->SetFixedParameters(_centerOfRotation);
        } else {
            metric->SetFixedImageRegion(_dst->GetBufferedRegion());
            metric->SetUseAllPixels(true);
            _centerOfRotation.SetSize(ImageType::ImageDimension);
            for (int i = 0; i < 3; i++) {
                _centerOfRotation[i] = _dstCenter[i];
            }
            _transform->SetFixedParameters(_centerOfRotation);
        }

        cout << "Fixed Parameters: " << _centerOfRotation << endl;

        metric->SetMovingImage(_src);
        metric->SetInterpolator(Interpolator::New());
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
                        scales[3*i+j] = 160;
                    } else {
                        scales[3*i+j] = 30;
                    }
                }
            }
            scales[9] = scales[10] = scales[11] = 0.1;
        } else if (_method == "scale") {
            scales[0] = scales[1] = scales[2] = 30;
            scales[3] = scales[4] = scales[5] = .5;
            scales[6] = scales[7] = scales[8] = 100;
        } else if (_method == "similar") {
            scales[0] = scales[1] = scales[2] = 10;
            scales[3] = scales[4] = scales[5] = 0.5;
            scales[6] = 100;
        }
        opti->SetScales(scales);

        const int maxIters = 100;
#ifdef USE_CG_OPTIMIZER
        opti->SetMaximumIteration(maxIters);
        opti->SetMaximumLineIteration(10);
        opti->SetUseUnitLengthGradient(true);
        opti->SetStepLength(1);
        opti->SetToFletchReeves();
#endif

#ifdef USE_GD_OPTIMIZER
        opti->SetNumberOfIterations(maxIters);
        opti->SetMinimumStepLength(1e-4);
        opti->SetMaximumStepLength(3);
        opti->SetRelaxationFactor(.5);
        opti->SetGradientMagnitudeTolerance(1e-4);
#endif

        opti->SetInitialPosition(_transform->GetParameters());
        opti->AddObserver(itk::StartEvent(), optiReporter);
        opti->AddObserver(itk::IterationEvent(), optiReporter);
        opti->StartOptimization();
        cout << "Current Cost: " << opti->GetValue() << endl;
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
    } else if (method == "scale") {
        RegistrationEngine<itk::ScaleVersor3DTransform<double> > reg;
        reg.main(argc, argv, method);
    } else if (method == "similar") {
        RegistrationEngine<itk::Similarity3DTransform<double> > reg;
        reg.main(argc, argv, method);
    }
}