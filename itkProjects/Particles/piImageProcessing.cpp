//
//  myImageProcessing.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/11/13.
//
//

#include "piImageProcessing.h"
#include "itkAntiAliasBinaryImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryCrossStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkEllipseSpatialObject.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkMesh.h"
#include "itkVTKPolyDataWriter.h"
#include "itkUnaryFunctorImageFilter.h"
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkMesh.h"
#include "itkBinaryMask3DMeshSource.h"
#include <itkStatisticsImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include "itkMultiResolutionPyramidImageFilter.h"
#include <itkManifoldParzenWindowsPointSetFunction.h>
#include <itkCannyEdgeDetectionImageFilter.h>
#include <itkInvertIntensityImageFilter.h>
#include <itkCenteredRigid2DTransform.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkExtractImageFilter.h>

#include "piImageIO.h"
#include "piImageHistogram.h"
#include "piPowellOpti.h"
#include "piOptionParser.h"
#include "piConfigFile.h"

#if DIMENSIONS == 3
#include "itkBinaryMask3DMeshSource.h"
#endif

namespace pi {
    typedef itk::EllipseSpatialObject<__Dim> EllipseType;
    typedef itk::SpatialObjectToImageFilter<EllipseType, LabelImage> SpatialObjectToImageFilterType;

    ImageIO<RealImage> __realIO;
    ImageIO<LabelImage> __labelIO;

#if DIMENSIONS == 3
    typedef itk::Mesh<PointReal> MeshType;
    typedef itk::BinaryMask3DMeshSource<LabelImage, MeshType> MeshSourceType;
#endif

    typedef itk::ZeroCrossingImageFilter<LabelImage, LabelImage> ZeroCrossingFilterType;
    typedef itk::CastImageFilter<LabelImage ,RealImage > CastToRealFilterType;
    typedef itk::AntiAliasBinaryImageFilter<RealImage, RealImage> AntiAliasFilterType;
    typedef itk::RescaleIntensityImageFilter<RealImage, LabelImage > RescaleFilter;
    typedef itk::ConnectedComponentImageFilter<LabelImage, LabelImage> ConnectedComponentFilterType;
    typedef itk::RelabelComponentImageFilter<LabelImage, LabelImage> RelabelFilterType;
    typedef itk::BinaryThresholdImageFilter<LabelImage, LabelImage> BinaryThreshFilterType;
    typedef itk::BinaryBallStructuringElement<LabelImage::PixelType,__Dim> StructuringElementType;
    typedef itk::BinaryCrossStructuringElement<LabelImage::PixelType,__Dim> AsymStructuringElementType;
    typedef itk::BinaryDilateImageFilter<LabelImage, LabelImage, StructuringElementType> DilateFilterType;
    typedef itk::BinaryErodeImageFilter<LabelImage, LabelImage, StructuringElementType> ErodeFilterType;
    typedef itk::SubtractImageFilter<LabelImage, LabelImage, LabelImage> SubtractFilterType;
    typedef itk::RecursiveGaussianImageFilter<RealImage> GaussianFilterType;
    typedef itk::SignedDanielssonDistanceMapImageFilter<LabelImage, RealImage> SignedDistanceMapFilterType;
    typedef itk::DanielssonDistanceMapImageFilter<LabelImage, RealImage> DistanceMapFilterType;
    typedef DistanceMapFilterType::VectorImageType DistanceVectorImageType;
    typedef DistanceVectorImageType::PixelType DistancePixelType;
    typedef itk::VectorMagnitudeImageFilter<VectorImage, RealImage> GradientMagnitudeFilterType;
    typedef itk::VectorMagnitudeImageFilter<VectorImage2, RealImage2> Gradient2MagnitudeFilterType;
    typedef itk::StatisticsImageFilter<RealImage> RealImageStatisticsFilterType;
    typedef itk::ResampleImageFilter<RealImage, RealImage> ResampleImageFilterType;



    static void die() {
        exit(EXIT_FAILURE);
    }

    static void end() {
        exit(EXIT_SUCCESS);
    }

    class OffsetToVector {
    public:
        bool operator!=(const OffsetToVector& o) {
            return false;
        }
        bool operator==(const OffsetToVector& o) {
            return true;
        }
        inline VectorType operator()(SignedDistanceMapFilterType::VectorImageType::PixelType o) {
            VectorType v;
            fordim (k) {
                v[k] = o[k];
            }
            return v;
        }
    };


    RealImageVector ImageProcessing::ComputeImagePyramid(RealImage::Pointer img, int level) {
        RealImageVector outputs;

        typedef itk::MultiResolutionPyramidImageFilter<RealImage, RealImage> FilterType;
        FilterType::Pointer filter = FilterType::New();
        filter->SetInput(img);
        filter->SetNumberOfLevels(level);
        filter->UseNearestNeighborInterpolatorOn();
        filter->Update();
        for (int i = level - 1; i >= 0; i--) {
            outputs.push_back(filter->GetOutput(i));
        }
        return outputs;
    }

    LabelImageVector ImageProcessing::ComputeImagePyramid(LabelImage::Pointer label, int level) {
        LabelImageVector outputs;

        typedef itk::MultiResolutionPyramidImageFilter<LabelImage, LabelImage> FilterType;
        FilterType::Pointer filter = FilterType::New();
        filter->SetInput(label);
        filter->SetNumberOfLevels(level);
        filter->UseNearestNeighborInterpolatorOff();
        filter->Update();
        for (int i = level - 1; i >= 0; i--) {
            outputs.push_back(filter->GetOutput(i));
        }
        return outputs;
    }




    
    LabelImage::Pointer ImageProcessing::SmoothLabelMap(LabelImage::Pointer img) {
#if 1
        return img;
#else
        cout << "Computing label map smoothing ..." << flush;

        LabelImage::Pointer procImage;

        const int LABEL_VAL = 255;
        CastToRealFilterType::Pointer toReal = CastToRealFilterType::New();
        toReal->SetInput(img);
        toReal->Update();

        AntiAliasFilterType::Pointer antiAliasFilter = AntiAliasFilterType::New();
        antiAliasFilter->SetInput( toReal->GetOutput() );
        antiAliasFilter->SetMaximumRMSError(1);
        antiAliasFilter->SetNumberOfIterations(10);
        antiAliasFilter->SetNumberOfLayers(2);
        antiAliasFilter->Update();

//        itkcmds::itkImageIO<RealImage> iox;
//        iox.WriteImageT("/tmpfs/aadouble.nrrd", antiAliasFilter->GetOutput());

        RescaleFilter::Pointer rescale = RescaleFilter::New();
        rescale->SetInput( antiAliasFilter->GetOutput() );
        rescale->SetOutputMinimum(   0 );
        rescale->SetOutputMaximum( 255 );
        rescale->Update();
        procImage = rescale->GetOutput();

        BinaryThreshFilterType::Pointer binTreshFilter = BinaryThreshFilterType::New();
        binTreshFilter->SetInput(rescale->GetOutput());
        binTreshFilter->SetUpperThreshold(255);
        binTreshFilter->SetLowerThreshold(127);
        binTreshFilter->SetOutsideValue (0);
        binTreshFilter->SetInsideValue (1);
        binTreshFilter->Update();
        procImage = binTreshFilter->GetOutput();

        ImageIO<LabelImage> io;
        io.WriteImage("./proc.nii.gz", procImage);

        ConnectedComponentFilterType::Pointer concompFilter = ConnectedComponentFilterType::New();
        RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
        BinaryThreshFilterType::Pointer binTreshFilter2 = BinaryThreshFilterType::New();
        concompFilter->SetInput(procImage);
        concompFilter->Update();
        procImage = concompFilter->GetOutput();

        relabelFilter->SetInput(concompFilter->GetOutput());
        relabelFilter->Update();
        procImage = relabelFilter->GetOutput();

        binTreshFilter2->SetInput(relabelFilter->GetOutput());
        binTreshFilter2->SetUpperThreshold(1);
        binTreshFilter2->SetLowerThreshold(1);// only largest object (labels are sorted largest to smallest)
        binTreshFilter2->SetOutsideValue (0);
        binTreshFilter2->SetInsideValue (LABEL_VAL);
        binTreshFilter2->Update();
        procImage   = binTreshFilter2->GetOutput();

        StructuringElementType structuringElement;
        structuringElement.SetRadius(1);  // 3x3x3 structuring element
        structuringElement.CreateStructuringElement();

        DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
        ErodeFilterType::Pointer erodeFilter = ErodeFilterType::New();

        dilateFilter->SetDilateValue (LABEL_VAL);
        dilateFilter->SetKernel(structuringElement);
        dilateFilter->SetInput(procImage);
        dilateFilter->Update();

        erodeFilter->SetErodeValue(LABEL_VAL);
        erodeFilter->SetKernel(structuringElement);
        erodeFilter->SetInput(dilateFilter->GetOutput());
        erodeFilter->Update();

        ErodeFilterType::Pointer erodeFilter2 = ErodeFilterType::New();
        erodeFilter2->SetErodeValue(LABEL_VAL);
        erodeFilter2->SetKernel(structuringElement);
        erodeFilter2->SetInput(erodeFilter->GetOutput());
        erodeFilter2->Update();

        procImage = erodeFilter2->GetOutput();

        cout << "done" << endl;
        return procImage;
#endif
    }

    LabelImage::Pointer ImageProcessing::ErodedBorder(pi::LabelImage::Pointer img) {
        LabelImage::Pointer procImage;

        StructuringElementType structuringElement;
        structuringElement.SetRadius(1);  // 3x3x3 structuring element
        structuringElement.CreateStructuringElement();

        ErodeFilterType::Pointer erodeFilter = ErodeFilterType::New();
        erodeFilter->SetErodeValue(255);
        erodeFilter->SetKernel(structuringElement);
        erodeFilter->SetInput(img);
        erodeFilter->Update();
        procImage = erodeFilter->GetOutput();

        SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
        subFilter->SetInput1(img);
        subFilter->SetInput2(procImage);
        subFilter->Update();

        return subFilter->GetOutput();
    }

    GradientImage::Pointer ImageProcessing::ComputeGaussianGradient(LabelImage::Pointer img, double sigma) {
        CastToRealFilterType::Pointer toReal = CastToRealFilterType::New();
        toReal->SetInput(img);
        toReal->Update();

        GaussianGradientFilterType::Pointer gradientFilter = GaussianGradientFilterType::New();
        gradientFilter->SetInput(toReal->GetOutput());
        if (sigma > -1) {
            gradientFilter->SetSigma(sigma);
        }
        gradientFilter->Update();
        return gradientFilter->GetOutput();
    }
    
    GradientImage::Pointer ImageProcessing::ComputeGradient(LabelImage::Pointer img) {
        CastToRealFilterType::Pointer toReal = CastToRealFilterType::New();
        toReal->SetInput(img);
        toReal->Update();

        GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
        gradientFilter->SetInput(toReal->GetOutput());
        gradientFilter->Update();
        return gradientFilter->GetOutput();
    }

    GradientImage::Pointer ImageProcessing::ComputeGaussianGradient(RealImage::Pointer img, double sigma) {
        GaussianGradientFilterType::Pointer gradientFilter = GaussianGradientFilterType::New();
        gradientFilter->SetInput(img);
        if (sigma > -1) {
            gradientFilter->SetSigma(sigma);
        }
        gradientFilter->Update();
        return gradientFilter->GetOutput();
    }

    GradientImage::Pointer ImageProcessing::ComputeGradient(RealImage::Pointer img) {
        GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
        gradientFilter->SetInput(img);
        gradientFilter->Update();
        return gradientFilter->GetOutput();
    }

    RealImage::Pointer ImageProcessing::ComputeGaussianGradientMagnitude(RealImage::Pointer img, double sigma) {
        typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<RealImage> FilterType;
        FilterType::Pointer filter = FilterType::New();
        if (sigma > 0) {
            filter->SetSigma(sigma);
        }
        filter->SetInput(img);
        filter->Update();
        return filter->GetOutput();
    }

    LabelImage::Pointer ImageProcessing::Ellipse(int* outputSize, double *center, double *radius) {
        SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();
        RealImage::SizeType size;
        fordim (k) {
            size[k] = outputSize[k];
        }
        imageFilter->SetSize(size);

        EllipseType::Pointer ellipse = EllipseType::New();
        ellipse->SetDefaultInsideValue(255);
        ellipse->SetDefaultOutsideValue(0);

        EllipseType::ArrayType axes;
        fordim (k) {
            axes[k] = radius[k];
        }
        ellipse->SetRadius(axes);

        EllipseType::TransformType::Pointer transform = EllipseType::TransformType::New();
        transform->SetIdentity();
        TransformType::OutputVectorType translation;
        fordim (k) {
            translation[k] = center[k];
        }
        transform->Translate(translation, false);

        ellipse->SetObjectToParentTransform(transform);
        imageFilter->SetInput(ellipse);
        imageFilter->SetUseObjectValue(true);
        imageFilter->SetOutsideValue(0);
        imageFilter->Update();
        return imageFilter->GetOutput();
    }


    VectorImage::Pointer ImageProcessing::ComputeDistanceMap(LabelImage::Pointer img) {
        cout << "Computing distance map ..." << flush;
        LabelImage::Pointer binaryMap = ThresholdToBinary(img);

        // construct signed distance filter
        SignedDistanceMapFilterType::Pointer distmapFilter = SignedDistanceMapFilterType::New();
        distmapFilter->SetInput(binaryMap);
        distmapFilter->Update();
        SignedDistanceMapFilterType::OutputImagePointer distmap = distmapFilter->GetOutput();

        typedef
        itk::UnaryFunctorImageFilter<SignedDistanceMapFilterType::VectorImageType, VectorImage, OffsetToVector> OffsetToVectorCastFilterType;

        OffsetToVectorCastFilterType::Pointer caster = OffsetToVectorCastFilterType::New();
        caster->SetInput(distmapFilter->GetVectorDistanceMap());
        caster->Update();
        cout << "done" << endl;
        return caster->GetOutput();
    }

    RealImage::Pointer ImageProcessing::ComputeMagnitudeMap(VectorImage::Pointer img) {
        GradientMagnitudeFilterType::Pointer magFilter = GradientMagnitudeFilterType::New();
        magFilter->SetInput(img);
        magFilter->Update();
        return magFilter->GetOutput();
    }

    RealImage2::Pointer ImageProcessing::ComputeMagnitude2Map(VectorImage2::Pointer img) {
        Gradient2MagnitudeFilterType::Pointer magFilter = Gradient2MagnitudeFilterType::New();
        magFilter->SetInput(img);
        magFilter->Update();
        return magFilter->GetOutput();
    }

    RealImage::Pointer ImageProcessing::ComputeMagnitudeMap(GradientImage::Pointer img) {
        typedef itk::VectorMagnitudeImageFilter<GradientImage, RealImage> MagnitudeFilterType;
        MagnitudeFilterType::Pointer magFilter = MagnitudeFilterType::New();
        magFilter->SetInput(img);
        magFilter->Update();
        return magFilter->GetOutput();
    }

    RealImage::Pointer ImageProcessing::ComputeGaussianSmoothing(RealImage::Pointer img, double sigma) {
        typedef itk::SmoothingRecursiveGaussianImageFilter<RealImage,RealImage> GaussianFilterType;
        GaussianFilterType::Pointer filter = GaussianFilterType::New();
        filter->SetInput(img);
        filter->SetSigma(sigma);
        filter->Update();
        RealImage::Pointer output = filter->GetOutput();
        output->DisconnectPipeline();
        return output;
    }

    vtkPolyData* ImageProcessing::ConvertToMesh(LabelImage::Pointer image) {
#if DIMENSIONS == 3
        typedef itk::Mesh<double> MeshType;
        typedef itk::BinaryMask3DMeshSource<LabelImage, MeshType> MeshSourceType;

        MeshSourceType::Pointer meshSource = MeshSourceType::New();
        meshSource->SetInput(image);
        meshSource->SetObjectValue(255);
        meshSource->Update();
        MeshType::Pointer mesh = meshSource->GetOutput();
#endif
        return NULL;
    }


    LabelImage::Pointer ImageProcessing::ZeroCrossing(LabelImage::Pointer src) {
        typedef itk::ZeroCrossingImageFilter<LabelImage, LabelImage> ZeroFilter;
        ZeroFilter::Pointer zero = ZeroFilter::New();
        zero->SetInput(src);
        zero->Update();
        return zero->GetOutput();
    }

    LabelImage::Pointer ImageProcessing::ComputeZeroCrossing(LabelImage::Pointer src) {
        // zero-crossing filter
        typedef itk::InvertIntensityImageFilter<LabelImage,LabelImage> InvertFilterType;
        InvertFilterType::Pointer invertFilter = InvertFilterType::New();
        invertFilter->SetInput(src);
        invertFilter->SetMaximum(255);
        invertFilter->Update();
        LabelImage::Pointer invertedMask = invertFilter->GetOutput();

        typedef itk::ZeroCrossingImageFilter<LabelImage,LabelImage> FilterType;
        FilterType::Pointer filter = FilterType::New();
        filter->SetInput(invertedMask);
        filter->SetForegroundValue(1);
        filter->SetBackgroundValue(0);
        filter->Update();
        return filter->GetOutput();
    }
    
    LabelImage::Pointer ImageProcessing::NormalizeToIntegralType(RealImage::Pointer src, LabelPixel min, LabelPixel max, LabelImage::Pointer mask) {
        ImageIO<LabelImage> io;
        LabelImage::Pointer output = io.CastImageFromS<RealImage>(src);
        output->FillBuffer(0);

        LabelImageIteratorType itermask(mask, mask->GetBufferedRegion());
        LabelImageIteratorType iterout(output, output->GetBufferedRegion());
        RealImageIteratorType iterin(src, src->GetBufferedRegion());

        iterout.GoToBegin();
        iterin.GoToBegin();

        RealImageStatisticsFilterType::Pointer filter = RealImageStatisticsFilterType::New();
        filter->SetInput(src);
        filter->Update();
        DataReal inMax = filter->GetMaximum();
        DataReal inMin = filter->GetMinimum();

        while (!iterout.IsAtEnd()) {
            if (mask.IsNull() || itermask.Get() > 0) {
                DataReal in = iterin.Get();
                LabelPixel out = (max - min) * (in - inMin) / (inMax - inMin) + min;
                iterout.Set(out);
            }
            ++iterout;
            ++iterin;
            ++itermask;
        }

        return output;
    }

    RealImage::Pointer ImageProcessing::NormalizeIntensity(RealImage::Pointer image, LabelImage::Pointer label, double percentile) {
        ImageIO<RealImage> io;
        RealImage::Pointer output = io.CopyImage(image);

        ImageHistogram<RealImage> histo;
        histo.fitPercentile = percentile;
        histo.SetImage(image);

        DataReal* inBuf = image->GetBufferPointer();
        DataReal* outBuf = output->GetBufferPointer();
        const int nPixels = image->GetPixelContainer()->Size();

        for (int i = 0; i < nPixels; i++) {
            outBuf[i] = (histo.rangeMax - histo.rangeMin) * inBuf[i] / (histo.dataMax - histo.dataMin) + histo.rangeMin;
        }

        ImageHistogram<RealImage> histo2;
        histo2.SetImage(output);
        for (int i = 0; i < nPixels; i++) {
            outBuf[i] = histo2.NormalizePixel(outBuf[i]);
        }
        return output;
    }

    LabelImage::Pointer ImageProcessing::ThresholdToBinary(LabelImage::Pointer img) {
        BinaryThreshFilterType::Pointer binThreshFilter = BinaryThreshFilterType::New();
        binThreshFilter->SetInput(img);
        binThreshFilter->SetInsideValue(1);
        binThreshFilter->SetOutsideValue(0);
        binThreshFilter->SetLowerThreshold(1);
        binThreshFilter->SetUpperThreshold(255);
        binThreshFilter->Update();
        return binThreshFilter->GetOutput();
    }


    ImageHistogramFilterType::HistogramPointer ImageProcessing::ComputeHistogram(RealImage::Pointer real, int nbin, DataReal rmin, DataReal rmax) {
        return ImageHistogramFilterType::HistogramPointer();
    }

    string ImageProcessing::ComputeHistogramToString(RealImage::Pointer real, int nbin, DataReal rmin, DataReal rmax) {
        return "";
    }

    LabelImage::Pointer ImageProcessing::ExtractLabelFilter(LabelImage::Pointer inputImage, int extractLabel) {
        BinaryThreshFilterType::Pointer threshFilter = BinaryThreshFilterType::New();
        threshFilter->SetInput(inputImage);
        threshFilter->SetLowerThreshold(extractLabel);
        threshFilter->SetUpperThreshold(extractLabel);
        threshFilter->SetOutsideValue(0);
        threshFilter->SetInsideValue(1);
        threshFilter->Update();
        return threshFilter->GetOutput();
    }

    RealImage::Pointer ImageProcessing::ProcessFeatureDensityImage(RealImage::Pointer realImage, LabelImage::Pointer maskImage, float kernelSigma, float regularizationSigma) {
        typedef itk::PointSet<float, DIMENSIONS> PointSetType;
        typedef PointSetType::PointType PointType;
        typedef PointSetType::PointsContainerPointer PointsContainerPointer;


        // preprocessing edge generation
        typedef itk::CannyEdgeDetectionImageFilter<RealImage, RealImage> FilterType;
        FilterType::Pointer filter = FilterType::New();
        RealImage::SpacingType spacing = realImage->GetSpacing();
        FilterType::ArrayType var;
        for (int i = 0; i < spacing.Size(); i++) {
            var[i] = spacing[i];
        }
        filter->SetInput(realImage);
        filter->SetVariance(var);
        filter->Update();
        RealImage::Pointer edgeImage = filter->GetOutput();


        // point set generation
        PointSetType::Pointer pointSet = PointSetType::New();
        PointsContainerPointer points = pointSet->GetPoints();

        itk::ImageRegionIteratorWithIndex<RealImage> regionIter(edgeImage, edgeImage->GetBufferedRegion());
        regionIter.GoToBegin();
        RealImage::PointType point;

        int i = 0;
        for (; !regionIter.IsAtEnd(); ++regionIter) {
            RealImage::IndexType idx = regionIter.GetIndex();
            if (regionIter.Get() > 0 && maskImage->GetPixel(idx) > 0) {
                edgeImage->TransformIndexToPhysicalPoint(idx, point);
                points->InsertElement(i++, point);
            }
        }


        // point set PDF evaluation
        typedef itk::ManifoldParzenWindowsPointSetFunction<PointSetType> ParzenFuncType;
        ParzenFuncType::Pointer parzenFunc = ParzenFuncType::New();

        parzenFunc->SetKernelSigma(kernelSigma);
        parzenFunc->SetRegularizationSigma(regularizationSigma);
        parzenFunc->SetInputPointSet(pointSet);

        ImageIO<RealImage> realIO;
        RealImage::Pointer outputImage = realIO.NewImage(realImage);
        itk::ImageRegionIteratorWithIndex<RealImage> outputIter(outputImage, outputImage->GetBufferedRegion());
        outputIter.GoToBegin();
        for (; !outputIter.IsAtEnd(); ++outputIter) {
            RealImage::IndexType idx = outputIter.GetIndex();
            RealImage::PointType idxPoint;
            outputImage->TransformIndexToPhysicalPoint(idx, idxPoint);
            outputImage->SetPixel(idx, parzenFunc->Evaluate(idxPoint));
        }

        return outputImage;
    }

    
    void AtlasSimilarityScore::Add(LabelPixel a, LabelPixel b) {
        if (labelMap.size() <= a || labelMap.size() <= b) {
            labelMap.resize(std::max(a,b)+1);
        }
        labelMap[a].L = a;
        labelMap[b].L = b;
        if (a == b) {
            labelMap[a].AB++;
        } else {
            labelMap[a].A++;
            labelMap[b].B++;
        }
    }

    void AtlasSimilarityScore::Compute(LabelImage::Pointer a, LabelImage::Pointer b) {
        int na = a->GetPixelContainer()->Size();
        int nb = b->GetPixelContainer()->Size();

        if (na != nb) {
            cout << "the size of image is different!" << endl;
        }

        const LabelPixel* pa = a->GetBufferPointer();
        const LabelPixel* pb = b->GetBufferPointer();
        for (int i = 0; i < na; i++) {
            Add(*pa, *pb);
            pa++;
            pb++;
        }
    }

    ostream& operator<<(ostream& os, const LabelSimilarityScore& score) {
        if (score.L > 0) {
            os << score.L << " " << score.A << " " << score.AB << " " << score.B << " " << score.total() << " " << score.dice();
        }
        return os;
    }

    ostream& operator<<(ostream& os, const AtlasSimilarityScore& score) {
        for (int i = 0; i < score.labelMap.size(); i++) {
            if (score.labelMap[i].L > 0) {
                os << score.labelMap[i] << endl;
            }
        }
        return os;
    }


#pragma mark -

    /// this code section contains a set of functions used in command line execution
    /// these functions accept arguments parsed from the command line.
    static ImageIO<RealImage> __io;

    void ImageProcessing::main(pi::Options &opts, StringVector &args) {
        doGaussian(opts, args);
        doBlur2(opts, args);
        doEllipse(opts, args);
        doGradMag(opts, args);
        doTransform2(opts, args);
        computeBoundingBox(opts, args);
        doCrop(opts, args);
        doSlice(opts, args);
        deformImage(opts, args);

        // registration test
        doAffineReg(opts, args);
        doGradHist(opts, args);

        testGradHistReg(opts, args);
        return;
    }

    void ImageProcessing::doGaussian(pi::Options &opts, StringVector &args) {
        double sigma = opts.GetReal("--doGaussian", -1);
        if (sigma <= 0) {
            return;
        }

        if (args.size() < 2) {
            cout << "--doGaussian {sigma} input-image output-image" << endl;
            exit(EXIT_FAILURE);
        }

        RealImage::Pointer image = __io.ReadCastedImage(args[0]);
        if (image.IsNull()) {
            exit(EXIT_FAILURE);
        }

        GaussianFilterType::Pointer gaussFilter = GaussianFilterType::New();
        gaussFilter->SetInput(image);
        gaussFilter->SetSigma(sigma);
        gaussFilter->Update();
        __io.WriteImage(args[1], gaussFilter->GetOutput());
        exit(EXIT_SUCCESS);

        return;
    }


    void ImageProcessing::doBlur2(pi::Options &opts, StringVector &args) {
        if (!opts.GetBool("--doBlur2")) {
            return;
        }

        if (args.size() < 3) {
            cout << "--doBlur2 [input-real-image] [output-real-image] [radius]" << endl;
            exit(EXIT_FAILURE);
        }
        float blurRadius = atof(args[2].c_str());

        ImageIO<RealImage2> io;
        RealImage2::Pointer realImage = io.ReadImage(args[0]);

        typedef itk::SmoothingRecursiveGaussianImageFilter<RealImage2, RealImage2> GaussianFilter;
        GaussianFilter::Pointer filter = GaussianFilter::New();
        filter->SetInput(realImage);
        filter->SetSigma(blurRadius);
        filter->Update();
        RealImage2::Pointer outputImage = filter->GetOutput();
        
        io.WriteImage(args[1], outputImage);
        exit(EXIT_SUCCESS);
    }


    void ImageProcessing::doEllipse(pi::Options &opts, StringVector &args) {
        if (!opts.GetBool("--ellipse")) {
            return;
        }

        if (args.size() < __Dim * (3) + 1) {
            cout << "--ellipse output-image [image-size] [ellipse-center] [ellipse-radius] " << endl;
            exit(EXIT_FAILURE);
        }

        int size[__Dim];
        double center[__Dim], radius[__Dim];
        fordim (k) {
            size[k] = atoi(args[1+k].c_str());
            center[k] = atof(args[1+__Dim*1 + k].c_str());
            radius[k] = atof(args[1+__Dim*2 + k].c_str());
        }

        LabelImage::Pointer outputImage = Ellipse(size, center, radius);
        ImageIO<LabelImage> io;
        io.WriteImage(args[0], outputImage);

        exit(EXIT_SUCCESS);
    }

    void ImageProcessing::doGradMag(pi::Options &opts, StringVector &args) {
        ImageProcessing proc;

        if (!opts.GetBool("--gradmag")) {
            return;
        }

        if (args.size() == 0 && opts.GetConfigFile() != "") {
            ConfigFile config(opts.GetConfigFile());

            double sigma = config["gaussian-sigma"];
            std::vector<RealImage::Pointer> outputs;
            for (int i = 0; i < config.imageCount(); i++) {
                RealImage::Pointer output = proc.ComputeGaussianGradientMagnitude(config.image(i), sigma);
                outputs.push_back(output);
            }
            config.writeImages<RealImage>(outputs, "gradient-magnitude");

            exit(EXIT_SUCCESS);
        }

        if (args.size() < 2) {
            cout << "--gradmag input-image output-image [sigma]" << endl;
            exit(EXIT_FAILURE);
        }

        RealImage::Pointer input = __realIO.ReadCastedImage(args[0]);
        double sigma = 1;
        if (args.size() == 3) {
            sigma = atof(args[2].c_str());
        }


        RealImage::Pointer output = proc.ComputeGaussianGradientMagnitude(input, sigma);
        __realIO.WriteImage(args[1], output);

        exit(EXIT_SUCCESS);
    }


    void ImageProcessing::doTransform2(pi::Options &opts, StringVector &args) {
        if (!opts.GetBool("--transform2")) {
            return;
        }

        if (args.size() < 3) {
            cout << "--transform2 input-image output-image r tx ty" << endl;
            exit(EXIT_FAILURE);
        }

        RealImage::Pointer input = __realIO.ReadCastedImage(args[0]);
        double r = 0, tx = 0, ty = 0;
        if (args.size() >= 3) {
            r = atof(args[2].c_str());
        }
        if (args.size() >= 4) {
            tx = atof(args[3].c_str());
        }
        if (args.size() >= 5) {
            ty = atof(args[4].c_str());
        }


        itk::CenteredRigid2DTransform<>::Pointer transform = itk::CenteredRigid2DTransform<>::New();
        // set rotation degree
        transform->SetAngleInDegrees(r);

        // set rotation center
        RealImage::SizeType sz = input->GetBufferedRegion().GetSize();
        itk::CenteredRigid2DTransform<>::InputPointType center;
        RealIndex centerIdx;
        centerIdx[0] = sz[0] / 2;
        centerIdx[1] = sz[1] / 2;
        input->TransformContinuousIndexToPhysicalPoint(centerIdx, center);
        transform->SetCenter(center);
        cout << "Transform Center: " << center << endl;

        // set translation after rotation
        TransformType::OutputVectorType translation;
        translation[0] = tx; translation[1] = ty;
        transform->SetTranslation(translation);

        ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
        resampler->SetInput(input);
        resampler->SetTransform(transform->GetInverseTransform());
        resampler->SetUseReferenceImage(true);
        resampler->SetReferenceImage(input);
        resampler->SetInterpolator(LinearImageInterpolatorType::New());
        resampler->SetDefaultPixelValue(0);
        resampler->Update();
        RealImage::Pointer output = resampler->GetOutput();
        output->DisconnectPipeline();

        __realIO.WriteImage(args[1], output);
        exit(EXIT_SUCCESS);
    }

    void ImageProcessing::computeBoundingBox(pi::Options &opts, StringVector &args) {
        if (!opts.GetBool("--boundingbox")) {
            return;
        }
        if (args.size() < 2) {
            cout << "--boundingbox [label] [crop-coordinate-output] [input-label] ..." << endl;
            die();
        }

        int labelToCrop = atoi(args[0].c_str());
        cout << "Label: " << labelToCrop << endl;

        string outputCoord = args[1];

        typedef itk::LabelStatisticsImageFilter<LabelImage, LabelImage> LabelStatFilter;

        // compute maximum bounding box
        LabelImage::IndexType lowerIndex, upperIndex;
        for (int i = 2; i < args.size(); i++) {
            LabelImage::Pointer label = __labelIO.ReadImage(args[i]);
            LabelStatFilter::Pointer statFilter = LabelStatFilter::New();
            statFilter->SetLabelInput(label);
            statFilter->SetInput(label);
            statFilter->Update();
            cout << label << endl;
            cout << statFilter->GetRegion(labelToCrop) << endl;
            LabelStatFilter::BoundingBoxType box = statFilter->GetBoundingBox(labelToCrop);
            if (i == 2) {
                for (int k = 0; k < box.size() / 2; k++) {
                    lowerIndex[k] = box[2*k];
                    upperIndex[k] = box[2*k+1];
                }
            } else {
                for (int k = 0; k < box.size() / 2; k++) {
                    lowerIndex[k] = std::min(lowerIndex[k], box[2*k]);
                    upperIndex[k] = std::max(upperIndex[k], box[2*k+1]);
                }
            }
        }

        ofstream of(outputCoord);
        of << "[ ";
        for (int j = 0; j < lowerIndex.GetIndexDimension(); j++) {
            cout << lowerIndex[j] << endl;
            of << lowerIndex[j] << ", ";
        }
        for (int j = 0; j < upperIndex.GetIndexDimension(); j++) {
            of << upperIndex[j] << ", ";
        }
        of << "0 ]" << endl;
        of.close();
        end();
    }


    void ImageProcessing::doCrop(Options& opts, StringVector& args) {
        if (opts.GetString("--crop") == "") {
            return;
        }
        if (args.size() < 2 || args.size() % 2 != 0) {
            cout << "--crop crop.txt [--padding n] [input] [output] [input] [output] ..." << endl;
            die();
        }

        RealImage::IndexType regionIndexes[2];

        OptionParser json;
        json.read(opts.GetString("--crop"));

        if (json.size() == 5) {
            json.values(regionIndexes, 2, 2);
        } else if (json.size() == 7) {
            json.values(regionIndexes, 2, 3);
        }

        int paddingSize = opts.GetStringAsInt("--padding", 3);

        RealImage::RegionType region;
        region.SetIndex(regionIndexes[0]);
        region.SetUpperIndex(regionIndexes[1]);
        region.PadByRadius(paddingSize);

        cout << "Crop Region: " << region << endl;

        int nFiles = args.size() / 2;
        for (int j = 0; j < nFiles; j++) {
            ImageInfo imageInfo;
            RealImage::Pointer image = __realIO.ReadCastedImage(args[2*j], imageInfo);
            itk::ExtractImageFilter<RealImage,RealImage>::Pointer filter =  itk::ExtractImageFilter<RealImage,RealImage>::New();
            filter->SetInput(image);
            filter->SetExtractionRegion(region);
            filter->Update();
            RealImage::Pointer outputImage = filter->GetOutput();
            __realIO.WriteCastedImage(args[2*j+1], outputImage, imageInfo.componenttype);
        }
        end();
    }


    void ImageProcessing::doSlice(Options& opts, StringVector& args) {
        if (!opts.GetBool("--slice")) {
            return;
        }
        if (args.size() < 4) {
            cout << "--slice dim index imagefile outputfile" << endl;
            die();
        }

        int dim = atoi(args[0].c_str());
        int slice = atoi(args[1].c_str());

        ImageIO<RealImage3> io;
        ImageInfo info;
        RealImage3::Pointer image = io.ReadCastedImage(args[2], info);
        RealImage2Vector sliceImages = __realImageTools.sliceVolume(image, dim);

        if (slice < sliceImages.size()) {
            ImageIO<RealImage2> wio;
            wio.WriteCastedImage(args[3], sliceImages[slice], info.componenttype);
        } else {
            cout << "slice index is out of range" << endl;
        }
        end();
    }


    void ImageProcessing::deformImage(Options& opts, StringVector& args) {
        if (!opts.GetBool("--deform")) {
            return;
        }
        if (args.size() < 3) {
            cout << "--deform input-image input-deformfield output-image" << endl;
            die();
        }

        ImageIO<DisplacementFieldType> io;
        DisplacementFieldType::Pointer deformField = io.ReadImage(args[1]);
        FieldTransformType::Pointer transform = FieldTransformType::New();
        transform->SetDisplacementField(deformField);

        RealImage::Pointer inputImage = __realIO.ReadImage(args[0]);

        ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
        resampler->SetInput(inputImage);
        resampler->SetTransform(transform);
        resampler->SetUseReferenceImage(true);
        resampler->SetReferenceImage(inputImage);
        resampler->Update();
        __realIO.WriteImage(args[2], resampler->GetOutput());

        end();
    }
    



    class GradientHistogram {
    public:
        GradientHistogram(RealImage::Pointer image) {
            _image = image;
        }

        void compute() {
            _gradImage = _proc.ComputeGaussianGradient(_image);
            _gradMag = _proc.ComputeMagnitudeMap(_gradImage);

            const int nSize = _gradMag->GetPixelContainer()->Size();
            RealImage::PixelType* magBuff = _gradMag->GetBufferPointer();
            GradientImage::PixelType* gradBuff = _gradImage->GetBufferPointer();

            _hist.resize(180);

            const int N = _hist.size();
            _pdf.resize(N);
            _cdf.resize(N);

            std::fill(_hist.begin(), _hist.end(), 0);
            std::fill(_pdf.begin(), _pdf.end(), 0);
            std::fill(_cdf.begin(), _cdf.end(), 0);

            int count = 0;
            for (int i = 0; i < nSize; i++) {
                const double n = magBuff[i];
                if (magBuff[i] > 10) {
                    double y = gradBuff[i][1] / n;
                    double x = gradBuff[i][0] / n;
                    double r = atan2(y, x);
                    if (r < 0) {
                        r += (2*M_PI);
                    }
                    double a = r * 180 / M_PI;
                    int bin = (a + 0.5) / 360 * (_hist.size() - 1);
                    _hist[bin] ++;
                    count ++;
                }
            }

            // threshold for # of gradient pixels (feature pixels)
            if (count > 50) {
                static double mask[] = { 0.0003, 0.1065, 0.7866, 0.1065, 0.0003 };
                double sum = 0;
                for (int i = 0; i < N; i++) {
                    double conv = 0;
                    for (int j = 0; j < 5; j++) {
                        conv += (_hist[(i+j-2+N)%N] * mask[j-2]);
                    }
                    _pdf[i] = conv;
                    sum += conv;
                }
                // normalization to p.d.f and compute c.d.f
                for (int i = 0; i < N; i++) {
                    _pdf[i] = _pdf[i] / sum;
                    if (i > 0) {
                        _cdf[i] = _cdf[i-1] + _pdf[i];
                    } else {
                        _cdf[i] = _pdf[i];
                    }
                }
                cout << "sum = " << _cdf[N-1] << endl;
            }
        }

        double ssd(GradientHistogram& o) {
            if (o._cdf.size() != _cdf.size()) {
                cout << "# of bins mismatch" << endl;
                return -1;
            }

            const int N = _cdf.size();
            double ssd = 0;
            for (int i = 0; i < N; i++) {
                ssd += (_cdf[i] - o._cdf[i])*(_cdf[i] - o._cdf[i]);
            }
            return ssd;
        }

        double emd(GradientHistogram& o) {
            if (o._pdf.size() != _pdf.size()) {
                cout << "# of bins mismatch" << endl;
                return -1;
            }

            const int N = _pdf.size();
            std::vector<double> emd;
            emd.resize(N);
            std::fill(emd.begin(), emd.end(), 0);

            for (int i = 1; i < N; i++) {
                emd[i] = (_pdf[i] + emd[i-1]) - o._pdf[i];
            }

            double emdsum = 0;
            for (int i = 0; i < N; i++) {
                emdsum += std::abs(emd[i]);
            }
            return emdsum;
        }

        void print() {
            for (int i = 0; i < _hist.size(); i++) {
                cout << _pdf[i] << "\t";
            }
            cout << endl;
        }

    private:
        std::vector<int> _hist;
        std::vector<double> _pdf, _cdf;
        RealImage::Pointer _image;
        GradientImage::Pointer _gradImage;
        RealImage::Pointer _gradMag;
        ImageProcessing _proc;
    };

    class GradHistReg {
    public:
        GradHistReg(RealImage::Pointer f, RealImage::Pointer m): _gf(f) {
            _f = f;
            _m = m;
            _gf.compute();

            RealImage::SizeType sz = f->GetBufferedRegion().GetSize();
            _transform = TransformType::New();
            RealIndex centerIdx;
            centerIdx[0] = sz[0] / 2;
            centerIdx[1] = sz[1] / 2;
            f->TransformContinuousIndexToPhysicalPoint(centerIdx, _center);
            _transform->SetCenter(_center);
            cout << "Transform Center: " << _center << endl;
        }

        double operator()(int n, double* x) {
            RealImage::Pointer w = warpImage(n, x);
            GradientHistogram gw(w);
            gw.compute();
            double metric = 0.7 * _gf.ssd(gw) + 0.3 * msd(_f, w);
            cout << x[0] << ", " << x[1] << ", " << x[2] << ", " << metric << endl;
            return metric;
        }

        double msd(RealImage::Pointer f, RealImage::Pointer w) {
            int N = f->GetPixelContainer()->Size();
            RealImage::PixelType* fptr = f->GetBufferPointer();
            RealImage::PixelType* wptr = w->GetBufferPointer();
            double ssd = 0;
            for (int i = 0; i < N; i++) {
                ssd += ((fptr[i] - wptr[i])*(fptr[i] - wptr[i]));
            }
            return sqrt(ssd) / N;
        }

        RealImage::Pointer warpImage(int n, double *x) {
            _transform->SetAngleInDegrees(x[0]);

            TransformType::OutputVectorType tx;
            tx[0] = x[1]; tx[1] = x[2];
            _transform->SetTranslation(tx);

            ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
            resampler->SetInput(_m);
            resampler->SetTransform(_transform->GetInverseTransform());
            resampler->SetUseReferenceImage(true);
            resampler->SetReferenceImage(_f);
            resampler->SetInterpolator(LinearImageInterpolatorType::New());
            resampler->SetDefaultPixelValue(0);
            resampler->Update();
            RealImage::Pointer output = resampler->GetOutput();
            output->DisconnectPipeline();
            return output;
        }
    private:
        RealImage::Pointer _f;
        RealImage::Pointer _m;
        GradientHistogram _gf;

        typedef itk::CenteredRigid2DTransform<double> TransformType;
        TransformType::Pointer _transform;
        TransformType::InputPointType _center;

    };

    class AffineReg {
    public:
        AffineReg(RealImage::Pointer f, RealImage::Pointer m) {
            fixedImage = f;
            movingImage = m;
        }
        double operator()(int n, double *x) {
            double ssd = 0;
            cout << "x = " << flush;
            for (int i = 0; i < n; i++) {
                cout << x[i] << ", ";
            }
            RealImage::Pointer warpedImage = warpImage(n, x);
            RealImage::PixelType* inBuf = fixedImage->GetBufferPointer();
            RealImage::PixelType* outBuf = warpedImage->GetBufferPointer();
            const int nPixels = fixedImage->GetPixelContainer()->Size();
            for (int i = 0; i < nPixels; i++) {
                ssd += ((inBuf[i] - outBuf[i])*(inBuf[i] - outBuf[i]));
            }
            cout << "ssd = " << ssd << endl;
            return ssd;
        }
        RealImage::Pointer warpImage(int n, double *x) {
            typedef itk::AffineTransform<double, __Dim> Transform;
            Transform::Pointer affine = Transform::New();
            Transform::ParametersType params;
            params.SetSize(__Dim * __Dim + __Dim);
            params.Fill(0);
            if (__Dim == 2) {
                double theta = x[0]/18.0*M_1_PI;
                params[0] = cos(theta); params[1] = sin(theta);
                params[2] = -sin(theta); params[3] = cos(theta);
                params[4] = x[1]; params[5] = x[2];
            } else {
                cout << "not supported dimension" << endl;
            }
            affine->SetParameters(params);

            ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
            resampler->SetInput(movingImage);
            resampler->SetTransform(affine->GetInverseTransform());
            resampler->SetUseReferenceImage(true);
            resampler->SetReferenceImage(fixedImage);
            resampler->SetInterpolator(LinearImageInterpolatorType::New());
            resampler->SetDefaultPixelValue(0);
            resampler->Update();
            RealImage::Pointer output = resampler->GetOutput();
            output->DisconnectPipeline();
            return output;
        }

    private:
        RealImage::Pointer fixedImage;
        RealImage::Pointer movingImage;
    };

    void ImageProcessing::doGradHist(pi::Options &opts, StringVector &args) {
        if (!opts.GetBool("--gradhist")) {
            return;
        }

        if (args.size() < 1) {
            cout << "--gradhist input-image" << endl;
            exit(EXIT_FAILURE);
        }

        RealImage::Pointer image = __realIO.ReadCastedImage(args[0]);

        GradientHistogram gradHist(image);
        gradHist.compute();
        gradHist.print();

        exit(EXIT_SUCCESS);
    }

    void ImageProcessing::doAffineReg(pi::Options &opts, StringVector &args) {
        if (!opts.GetBool("--affineReg")) {
            return;
        }

        if (args.size() < 3) {
            cout << "--affineReg fixed-image moving-image warped-image-output" << endl;
            exit(EXIT_FAILURE);
        }

        RealImage::Pointer fixedImage = __realIO.ReadCastedImage(args[0]);
        RealImage::Pointer movingImage = __realIO.ReadCastedImage(args[1]);

        PowellOpti<AffineReg> opti;
        PowellParams initial;
        initial.resize(3);
        initial[0] = initial[1] = initial[2] = 0;

        AffineReg ssd(fixedImage, movingImage);
        opti.minimizeNEWUOA(ssd, initial, 10);

        for (int i = 0; i < initial.size(); i++) {
            cout << initial[i] << endl;
        }

        RealImage::Pointer warpedImage = ssd.warpImage(initial.size(), &initial[0]);
        __realIO.WriteImage(args[2], warpedImage);

        exit(EXIT_SUCCESS);
    }

    void ImageProcessing::testGradHistReg(pi::Options &opts, StringVector &args) {
        if (!opts.GetBool("--testgradreg")) {
            return;
        }

        if (args.size() < 1) {
            cout << "--testgradreg fixed-image moving-image output-image" << endl;
            exit(EXIT_FAILURE);
        }

        RealImage::Pointer f = __realIO.ReadCastedImage(args[0]);
        RealImage::Pointer m = __realIO.ReadCastedImage(args[1]);

        PowellParams x0(3);
        x0[0] = x0[1] = x0[2] = 0;

        PowellOpti<GradHistReg> opti;
        GradHistReg reg(f, m);
        opti.minimizeNEWUOA(reg, x0, 50);

        for (int i = 0; i < x0.size(); i++) {
            cout << x0[i] << endl;
        }

        RealImage::Pointer w = reg.warpImage(x0.size(), &x0[0]);
        __realIO.WriteImage(args[2], w);
        
        exit(EXIT_SUCCESS);
    }
    

}