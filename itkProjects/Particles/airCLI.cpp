//
//  airCLI.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/1/13.
//
//

#include "airCLI.h"
#include "airImageAlgorithm.h"

#include <itkExtractImageFilter.h>
#include <itkGenerateImageSource.h>
#include <itkPasteImageFilter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkCenteredRigid2DTransform.h>
#include <itkTranslationTransform.h>
#include <itkCorrelationImageToImageMetricv4.h>
#include <itkBSplineTransform.h>
#include <itkMeanSquaresImageToImageMetricv4.h>
#include <itkMattesMutualInformationImageToImageMetricv4.h>
#include "itkEntropyImageToImageMetricv4.h"
#include <itkGradientDescentOptimizerv4.h>
#include <itkQuasiNewtonOptimizerv4.h>
#include <itkConjugateGradientLineSearchOptimizerv4.h>
#include <itkRegistrationParameterScalesFromIndexShift.h>
#include <itkVersorRigid3DTransform.h>
#include <itkSimilarity3DTransform.h>
#include <itkEuler3DTransform.h>
#include <itkAffineTransform.h>
#include <itkResampleImageFilter.h>
#include <itkCommand.h>

//#include "piParticleCore.h"
//#include "piParticleTrace.h"

#include <QString>

using namespace pi;
using namespace std;
using namespace itk;


namespace air {
    ImageIO<Image> __imageIO;
    ImageIO<Label> __labelIO;
    ImageIO<ImageSlice> __imageSliceIO;
    ImageIO<LabelSlice> __labelSliceIO;

    typedef itk::ImageRegionIteratorWithIndex<Image> ImageIterator;
    typedef itk::ImageRegionIteratorWithIndex<ImageSlice> ImageSliceIterator;
    typedef itk::ImageRegionIteratorWithIndex<Label> LabelIterator;
    typedef itk::ImageRegionIteratorWithIndex<LabelSlice> LabelSliceIterator;
    typedef itk::LinearInterpolateImageFunction<Image> ImageInterpolator;
    typedef itk::NearestNeighborInterpolateImageFunction<Image> NNImageInterpolator;

    typedef air::ImageAlgorithm<Image,Label> ImageAlgo;

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

        virtual void Execute(Object *caller, const EventObject & event) {
            T* realCaller = dynamic_cast<T*>(caller);
            if (realCaller == NULL) {
                return;
            }
            cout << realCaller->GetCurrentIteration() << ": " << realCaller->GetCurrentPosition() << endl;
        }

        virtual void Execute(const Object *caller, const EventObject & event) {

        }


    protected:
        OptimizerProgress() {}

        virtual ~OptimizerProgress() {}

    private:
        OptimizerProgress(const Self &);        //purposely not implemented
        void operator=(const Self &); //purposely not implemented
    };

    
    void CommandLineTools::ExtractSlice(Image::Pointer img, pi::SliceDirectionEnum dir, pi::IntVector range, std::string outputPattern) {

        for (int i = range[0]; i <= range[1]; i++) {
            Image::RegionType extractRegion = img->GetBufferedRegion();
            if (i < 0 || i >= extractRegion.GetSize(dir)) {
                continue;
            }
            extractRegion.SetIndex(dir,i);
            extractRegion.SetSize(dir,0);
            typedef itk::ExtractImageFilter<Image, ImageSlice> Filter;
            Filter::Pointer filter = Filter::New();
            filter->SetInput(img);
            filter->SetExtractionRegion(extractRegion);
            filter->SetDirectionCollapseToGuess();
            filter->Update();
            QString filePattern = QString::fromStdString(outputPattern);
            QString fileOutput = filePattern.arg(i,3,10,QChar('0'));
            __imageSliceIO.WriteImage(fileOutput.toStdString(), filter->GetOutput());
        }

        return;
    }

    void CommandLineTools::PasteSlice(StringVector args, string output) {
        int nSlices = args.size() - 1;
        ImageInfo imageInfo;
        ImageSlice::Pointer firstImage = __imageSliceIO.ReadCastedImage(args[1], imageInfo);
        ImageSlice::RegionType firsrtRegion = firstImage->GetBufferedRegion();
        Image::Pointer pasteImage = __imageIO.NewImageT(firsrtRegion.GetSize(0), firsrtRegion.GetSize(1), nSlices);
        pasteImage->FillBuffer(0);

        ImageIterator pasteIter(pasteImage, pasteImage->GetBufferedRegion());
        pasteIter.GoToBegin();
        for (int i = 0; i < nSlices; i++) {
            cout << "Processing #" << i << ": " << args[i+1] << endl;

            ImageSlice::Pointer sliceImage = __imageSliceIO.ReadCastedImage(args[i+1]);
            if (firsrtRegion != sliceImage->GetBufferedRegion()) {
                cout << "Error Processing #" << i << ": " << args[i+1] << endl;
                return;
            }
            ImageSliceIterator sliceIter(sliceImage, sliceImage->GetBufferedRegion());
            sliceIter.GoToBegin();
            while (!sliceIter.IsAtEnd()) {
                pasteIter.Set(sliceIter.Get());
                ++pasteIter;
                ++sliceIter;
            }
        }

        __imageIO.WriteCastedImage(output, pasteImage, imageInfo.componenttype);
        return;
    }


    void CommandLineTools::PasteLabel(StringVector args, string output) {
        int nSlices = args.size() - 1;
        LabelSlice::Pointer firstImage = __labelSliceIO.ReadCastedImage(args[1]);
        LabelSlice::RegionType sliceRegion = firstImage->GetBufferedRegion();
        Label::Pointer pasteImage = __labelIO.NewImageT(sliceRegion.GetSize(0), sliceRegion.GetSize(1), nSlices);

        LabelIterator pasteIter(pasteImage, pasteImage->GetBufferedRegion());
        pasteIter.GoToBegin();
        for (int i = 0; i < nSlices; i++) {
            cout << "Processing #" << i << ": " << args[i+1] << endl;

            LabelSlice::Pointer sliceImage = __labelSliceIO.ReadCastedImage(args[i+1]);
            LabelSliceIterator sliceIter(sliceImage, sliceImage->GetBufferedRegion());
            sliceIter.GoToBegin();
            while (!sliceIter.IsAtEnd()) {
                pasteIter.Set(sliceIter.Get());
                ++pasteIter;
                ++sliceIter;
            }
        }
        __labelIO.WriteImage(output, pasteImage);
        return;
    }

    void CommandLineTools::RigidRegistration(pi::Options *opts, pi::StringVector args) {
        Image::Pointer fixedImage  = __imageIO.ReadCastedImage(args[0]);
        Image::Pointer movingImage = __imageIO.ReadCastedImage(args[1]);
        
//        typedef itk::TranslationTransform<double,3> TransformType;
//        typedef itk::VersorRigid3DTransform<double> TransformType;
//        typedef itk::Similarity3DTransform<double> TransformType;
        //        typedef itk::BSplineTransform<double,2,4> TransformType;
//        typedef itk::AffineTransform<double,3> TransformType;
        typedef itk::Euler3DTransform<> TransformType;
        TransformType::Pointer transform;

        typedef itk::GradientDescentOptimizerv4 OptimizerType;
        OptimizerType::Pointer optimizer;

        typedef itk::MattesMutualInformationImageToImageMetricv4<Image, Image> CostFunctionType;
        CostFunctionType::Pointer costFunc;

        typedef OptimizerType::ParametersType ParametersType;

        typedef itk::RegistrationParameterScalesFromIndexShift<CostFunctionType> ScaleEstimatorType;
        ScaleEstimatorType::Pointer estimator;

        transform = TransformType::New();

        ParametersType initialParams;
        initialParams.SetSize(transform->GetNumberOfParameters());
        initialParams.Fill(0);

        ParametersType fixedParams;
        fixedParams.SetSize(fixedImage->GetImageDimension());

        itk::ContinuousIndex<double,3> centerIdx;
        Image::PointType centerPoint;

        for (int i = 0; i < fixedParams.GetSize(); i++) {
            centerIdx[i] = fixedImage->GetBufferedRegion().GetIndex(i) + fixedImage->GetBufferedRegion().GetSize(i) / 2.0;
        }
        fixedImage->TransformContinuousIndexToPhysicalPoint(centerIdx, centerPoint);
        for (int i = 0; i < fixedParams.size(); i++) {
            fixedParams[i] = centerPoint[i];
        }
        transform->SetFixedParameters(fixedParams);

        costFunc = CostFunctionType::New();
        costFunc->SetFixedImage(fixedImage);
        costFunc->SetFixedInterpolator(ImageInterpolator::New());
        costFunc->SetMovingImage(movingImage);
        costFunc->SetMovingInterpolator(ImageInterpolator::New());
        costFunc->SetMovingTransform(transform);
        costFunc->SetParameters(initialParams);
        if (dynamic_cast<MattesMutualInformationImageToImageMetricv4<Image, Image>*>(costFunc.GetPointer()) != NULL) {
            costFunc->SetNumberOfHistogramBins(32);
        }
        costFunc->Initialize();

        OptimizerProgress<OptimizerType>::Pointer progress = OptimizerProgress<OptimizerType>::New();

        OptimizerType::ScalesType scales;
        scales.SetSize(transform->GetNumberOfParameters());
        scales.Fill(1);
        
        optimizer = OptimizerType::New();
        optimizer->SetScales(scales);
//        optimizer->SetScalesEstimator(estimator);
        optimizer->SetMetric(costFunc);
        optimizer->AddObserver(IterationEvent(), progress);

        Image::SpacingType spacing = fixedImage->GetSpacing();
        optimizer->SetMaximumStepSizeInPhysicalUnits(spacing[0]*3);

        try {
            ::itk::Object::GlobalWarningDisplayOn();
            optimizer->SetDebug(true);
            costFunc->SetDebug(true);
            optimizer->Print(cout);
            optimizer->StartOptimization();
        } catch (ExceptionObject& ex) {
            ex.Print(cout);
        }

        cout << "Iterations: " << optimizer->GetCurrentIteration() << endl;
        cout << "Stop Reason: " << optimizer->GetStopConditionDescription() << endl;

        const ParametersType& params = optimizer->GetCurrentPosition();
        transform->SetParameters(params);
        
        typedef itk::ResampleImageFilter<Image,Image> ResampleFilter;
        ResampleFilter::Pointer resampler = ResampleFilter::New();
        resampler->SetInput(movingImage);
        if (transform.IsNotNull()) {
            cout << "Setting transform:" << endl;
            transform->Print(cout);
            resampler->SetTransform(dynamic_cast<ResampleFilter::TransformType*>(transform.GetPointer()));
        }
        resampler->SetReferenceImage(fixedImage);
        resampler->UseReferenceImageOn();
        resampler->Update();
        
        Image::Pointer resampled = resampler->GetOutput();
        __imageIO.WriteImage(args[2], resampled);
        __imageIO.WriteSingleTransform(args[3].c_str(), transform);
    }


    int CommandLineTools::Run(Options* parser, StringVector args) {
        if (parser->GetBool("--isoRG")) {
            int fg = parser->GetStringAsInt("-f", 2);
            int bg = parser->GetStringAsInt("-b", 3);
            if (args.size() < 3) {
                cout << "Isolated Connected Filter: -f [fgId=2] -b [bgId=3] [input-image] [input-label] [output-label]" << endl;
                return 0;
            }

            ImageAlgo algo;
            Image::Pointer gray = __imageIO.ReadCastedImage(args[0]);
            Label::Pointer label = __labelIO.ReadCastedImage(args[1]);
            Label::Pointer labelOut = algo.ExecuteIsolatedConnectedImageFilter(label, gray, fg, bg);

            __labelIO.WriteImage(args[2], labelOut);
        } else if (parser->GetBool("--extractSlice")) {
            if (args.size() < 2) {
                cout << "ExtractSlice --dir IJ|JK|KI --range 10,100 [input-image] [output-image-pattern: slice-%1.nii.gz]" << endl;
                return 0;
            }

            Image::Pointer inputImg = __imageIO.ReadCastedImage(args[0]);
            if (inputImg.IsNull()) {
                cout << "can't read " << args[0] << endl;
                return 0;
            }
            string dirString = parser->GetString("--dir");

            SliceDirectionEnum dir = Unknown;
            if (dirString == "IJ") {
                dir = IJ;
            } else if (dirString == "JK") {
                dir = JK;
            } else if (dirString == "KI") {
                dir = KI;
            }

            IntVector rangeValues = parser->GetStringAsIntVector("--range");
            ExtractSlice(inputImg, dir, rangeValues, args[1]);
        } else if (parser->GetBool("--pasteSlice")) {
            if (args.size() < 2) {
                cout << "PasteSlice [output-volume] [input1] [input2] ... [inputN]" << endl;
                return 0;
            }
            PasteSlice(args, args[0]);
        } else if (parser->GetBool("--pasteLabel")) {
            if (args.size() < 2) {
                cout << "PasteLabel [output-volume] [input1] [input2] ... [inputN]" << endl;
                return 0;
            }
            PasteLabel(args, args[0]);
        } else if (parser->GetBool("--rigidRegister")) {
            RigidRegistration(parser, args);
        }
        return 0;
    }


}