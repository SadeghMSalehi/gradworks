//
//  piPatchCompare.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/13/13.
//
//

#include "piPatchCompare.h"
#include "piImageIO.h"

#include <vector>
#include <algorithm>

#include <itkImageRegistrationMethod.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkResampleImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkTranslationTransform.h>
#include <itkDemonsRegistrationFilter.h>
#include <vnl/vnl_linear_system.h>

#include "piPowellOpti.h"
#include "piImageProcessing.h"

namespace pi {
    struct PatchSimilarity {
        int originAtlas;
        int originLabel;
        double metric;

        PatchSimilarity(int o, int l, double s) {
            originAtlas = o;
            originLabel = l;
            metric = s;
        }

        PatchSimilarity(const PatchSimilarity& o) {
            originAtlas = o.originAtlas;
            originLabel = o.originLabel;
            metric = o.metric;
        }

        PatchSimilarity& operator=(const PatchSimilarity& o) {
            originAtlas = o.originAtlas;
            originLabel = o.originLabel;
            metric = o.metric;
            return (*this);
        }
    };


    bool patch_compare_descending(const PatchSimilarity& a, const PatchSimilarity& b) {
        return (a.metric < b.metric);
    }


    // main testing routines
    bool PatchCompare::main(pi::Options &parser, StringVector args) {
        if (parser.GetBool("--makePatch")) {
            ImageIO<RealImage> io;
            RealImage::Pointer image = io.ReadCastedImage(args[0]);
            PatchImage::Pointer patchImage = this->buildPatchImage(image);
            ImageIO<PatchImage> io2;
            io2.WriteImage(args[1], patchImage);
            return true;
        } else if (parser.GetBool("--makeGradientPatch")) {
            if (args.size() < 1 + DIMENSIONS) {
                cout << "--makeGradientPatch input-image gx-patch gy-patch gz-patch ..." << endl;
                return true;
            }
            ImageIO<RealImage> io;
            RealImage::Pointer image = io.ReadCastedImage(args[0]);
            PatchImageVector gradientPatchImages;
            buildGradientPatchImage(image, 5, gradientPatchImages);
            for (int i = 0; i < gradientPatchImages.size(); i++) {
                io.WriteImageS<PatchImage>(args[i+1], gradientPatchImages[i]);
            }
            return true;
        } else if (parser.GetBool("--opticalFlow")) {
            if (args.size() < 3) {
                cout << "--opticalFlow fixed-image moving-image output-image" << endl;
                return true;
            }


            ImageIO<RealImage> io;
            RealImage::Pointer fixed = io.ReadCastedImage(args[0]);
            RealImage::Pointer moving = io.ReadCastedImage(args[1]);

            DisplacementFieldType::Pointer flow = computeOpticalFlow(fixed, moving, 0.1);
            io.WriteImageS<DisplacementFieldType>(args[2], flow);
            return true;
        } else if (parser.GetBool("--labelTransferWithPatch")) {
            transferLabelsWithPatch(args, parser.GetString("-o"), parser.GetInt("--searchRadius", 3), parser.GetInt("--kNearest", 3));
            return true;
        } else if (parser.GetBool("--demons")) {
            if (args.size() < 3) {
                cout << "--demons fixed-image moving-image output-image" << endl;
                return true;
            }

            ImageIO<RealImage> io;
            RealImage::Pointer fixed = io.ReadCastedImage(args[0]);
            RealImage::Pointer moving = io.ReadCastedImage(args[1]);

            PatchImage::Pointer nullImage = PatchImage::Pointer(NULL);

            cout << "Perform demons registration ... " << flush;
            performDemonsRegistration(fixed, moving, args[2]);
            cout << "done" << endl;
            
            return true;
        } else if (parser.GetBool("--patchtest") && args.size() > 0) {
            if (args[0] == "interpolation") {
                testInterpolation(args);
            }
            return true;
        }
        return false;
    }

#pragma mark Test Codes
    void PatchCompare::testInterpolation(StringVector &args) {
        ImageIO<RealImage> io;

        if (args.size() < 4) {
            cout << "--patchtest interpolation [image] [x] [y] " << endl;
        }
        RealImage::Pointer image = io.ReadCastedImage(args[1]);
        PatchImage::Pointer patchImage = buildPatchImage(image);
        typedef itk::VectorLinearInterpolateImageFunction<PatchImage> Interpolator;

        double x = atof(args[2].c_str());
        double y = atof(args[3].c_str());

        Interpolator::Pointer intp = Interpolator::New();
        intp->SetInputImage(patchImage);

        ImageIO<RealImage3> io3;
        RealImage3::Pointer realImage = io3.NewImageT(5, 5, 10);

        itk::ImageRegionIteratorWithIndex<RealImage3> iter3(realImage, realImage->GetBufferedRegion());

        iter3.GoToBegin();

        RealIndex idx;
        for (double t= 0; t < 10; t++) {
            x += 0.1;
            y += 0.1;

            idx[0] = x;
            idx[1] = y;

            cout << x << endl;

            PatchImage::PixelType px = intp->EvaluateAtContinuousIndex(idx);
            for (int i = 0; i < 25; i++) {
                iter3.Set(px[i]);
                ++iter3;
            }
        }

        io3.WriteImage("test.nrrd", realImage);
    }

#pragma mark -
    void PatchCompare::setAtlasImages(PatchImageVector atlasImages) {
        _atlasImages = atlasImages;
    }

    void PatchCompare::setAtlasLabels(LabelImageVector atlasLabels) {
        _atlasLabels = atlasLabels;
    }

    void PatchCompare::setTargetImage(RealImage::Pointer target) {
        _targetImage = target;
    }

    LabelImage::Pointer PatchCompare::getTargetLabel() {
        return _targetLabel;
    }

    void PatchCompare::setTargetROI(LabelImage::Pointer target) {
        _targetROI = target;
    }

    void PatchCompare::setParticleSystem(pi::ParticleSystem *system) {
        _system = system;
    }

    void PatchCompare::setTargetRadius(int r) {
        _targetRadius = r;
    }

    // assume isotropic patch size
    void PatchCompare::setPatchRadius(int r) {
        _patchRadius = r;
    }

    // run estimation comparing patches
    void PatchCompare::estimateLabel(int searchRadius, int kNearest) {
        RealImage::RegionType targetRegion = _targetImage->GetBufferedRegion();
        PatchImage::Pointer targetPatch = buildPatchImage((_targetImage));

        ImageIO<LabelImage> labelIO;
        _targetLabel = labelIO.CastImageFromS<RealImage>(_targetImage);
        _targetLabel->FillBuffer(0);
        
        cout << "searchRadius = " << searchRadius << ", kNearest = " << kNearest << endl;
        
        itk::ImageRegionConstIteratorWithIndex<PatchImage> iter(targetPatch, targetRegion);
        // for each particle
        // for each voxel v1 inside the target radius
        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            RealImage::IndexType idx = iter.GetIndex();

            bool insideROI = _targetROI->GetPixel(idx) > 0;
            if (!insideROI) {
                continue;
            }

            // extract patch P1 at v1
            LocalPatch targetPatch = iter.Get();

            // for each voxel v2 inside the search radius of v1

            // define search area
            RealImage::RegionType searchRegion;
            RealImage::IndexType searchIndex = idx;
            RealImage::SizeType searchSize;
            searchSize.Fill(searchRadius);
            RealImage::OffsetType searchOffset;
            searchOffset.Fill(searchRadius/2);
            searchIndex -= searchOffset;
            searchRegion.SetIndex(idx);
            searchRegion.SetSize(searchSize);

            // for each atlas
            std::vector<PatchSimilarity> similarityRanking;
            for (unsigned int i = 0; i < _atlasImages.size(); i++) {
                itk::ImageRegionConstIteratorWithIndex<PatchImage> searchIter(_atlasImages[i], searchRegion);
                for (searchIter.GoToBegin(); !searchIter.IsAtEnd(); ++searchIter) {
                    // extract patch P2 at v2
                    RealImage::IndexType originIdx = searchIter.GetIndex();
                    LocalPatch atlasPatch = searchIter.Get();

                    // compute patch distance between P1 and P2
                    double ssd = 0;
                    for (unsigned int j = 0; j < targetPatch.Size(); j++) {
                        double diff = atlasPatch[j] - targetPatch[j];
                        ssd += (diff * diff);
                    }

                    int originLabel = _atlasLabels[i]->GetPixel(originIdx);
                    similarityRanking.push_back(PatchSimilarity(i, originLabel, ssd));
                }
                // end for
            }
            // end for

            // choose top k nearest patches
            std::sort(similarityRanking.begin(), similarityRanking.end(), patch_compare_descending);

            // estimate label for v1
            std::vector<int> counter(100);
            std::fill(counter.begin(), counter.end(), 0);
            for (unsigned int k = 0; k < kNearest; k++) {
                if (similarityRanking[k].originLabel > counter.size()) {
                    cout << "Too big label index (" << similarityRanking[k].originLabel << ") ignored" << endl;
                    continue;
                }
                counter[similarityRanking[k].originLabel] ++;
            }
            int estimatedLabel = std::distance(counter.begin(), std::max_element(counter.begin(), counter.end()));
            _targetLabel->SetPixel(idx, estimatedLabel);
        }
        // end for
    }

    PatchImage::Pointer PatchCompare::buildPatchImage(RealImage::Pointer image) {
        const int patchWidth = PATCH_SIZE;
        ImageIO<PatchImage> io;

        RealImage::RegionType region = image->GetBufferedRegion();
        RealImage::SizeType imageSize = region.GetSize();

        PatchImage::Pointer outputImage = io.NewImageS<RealImage>(image);
        LocalPatch zeroVector;
        zeroVector.Fill(0);
        outputImage->FillBuffer(zeroVector);

        // set region to contain full patch
        RealImage::IndexType regionIdx;
        regionIdx[0] = patchWidth/2 + 1;
        regionIdx[1] = patchWidth/2 + 1;

        RealImage::SizeType regionSize;
        regionSize[0] = imageSize[0] - patchWidth;
        regionSize[1] = imageSize[1] - patchWidth;

        region.SetIndex(regionIdx);
        region.SetSize(regionSize);

        // iterate over the region
        RealImageIteratorType iter(image, region);
        RealImage::RegionType patchRegion;
        RealImage::OffsetType patchOffset;
        patchOffset.Fill(patchWidth/2);
        RealImage::SizeType patchSize;
        patchSize.Fill(patchWidth);

        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            LocalPatch patch;
            // extract patch
            RealImage::IndexType idx = iter.GetIndex();
            idx -= patchOffset;
            patchRegion.SetIndex(idx);
            patchRegion.SetSize(patchSize);
            RealImageIteratorType patchIter(image, patchRegion);
            int i = 0;
            for (patchIter.GoToBegin(); !patchIter.IsAtEnd(); ++patchIter) {
                patch[i] = patchIter.Get();
                i++;
            }
            outputImage->SetPixel(idx, patch);
        }
        return outputImage;
    }

    void PatchCompare::buildGradientPatchImage(RealImage::Pointer image, const int radius, PatchImageVector& output) {

        LocalPatch zeroVector;
        zeroVector.Fill(0);

        ImageIO<PatchImage> io;
        for (int i = 0; i < RealImage::ImageDimension; i++) {
            PatchImage::Pointer patch = io.NewImageS<RealImage>(image);
            patch->FillBuffer(zeroVector);
            output.push_back(patch);
        }

        ImageProcessing proc;
        GradientImage::Pointer gradientImage = proc.ComputeGaussianGradient(image, 0.1);

        // set region to contain full patch
        RealImage::RegionType region = image->GetBufferedRegion();
        RealImage::SizeType imageSize = region.GetSize();

        RealImage::IndexType regionIdx;
        regionIdx[0] = radius/2 + 1;
        regionIdx[1] = radius/2 + 1;

        RealImage::SizeType regionSize;
        regionSize[0] = imageSize[0] - radius;
        regionSize[1] = imageSize[1] - radius;

        region.SetIndex(regionIdx);
        region.SetSize(regionSize);

        // iterate over the region
        itk::ImageRegionIteratorWithIndex<GradientImage> iter(gradientImage, region);
        RealImage::RegionType patchRegion;
        RealImage::OffsetType patchOffset;
        patchOffset.Fill(radius/2);

        RealImage::SizeType patchSize;
        patchSize.Fill(radius);

        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            LocalPatch patch;
            // extract patch
            RealImage::IndexType idx = iter.GetIndex();
            idx -= patchOffset;
            patchRegion.SetIndex(idx);
            patchRegion.SetSize(patchSize);
            itk::ImageRegionIteratorWithIndex<GradientImage> patchIter(gradientImage, patchRegion);
            int i = 0;
            for (patchIter.GoToBegin(); !patchIter.IsAtEnd(); ++patchIter) {
                GradientImage::PixelType grad = patchIter.Get();
                for (int j = 0; j < RealImage::ImageDimension; j++) {
                    output[j]->GetPixel(idx)[i] = grad[j];
                }
                i++;
            }
        }
    }


    DisplacementFieldType::Pointer PatchCompare::computeOpticalFlow(RealImage::Pointer fixed, RealImage::Pointer moving, double dt) {
        ImageIO<DisplacementFieldType> io;

        PatchImage::Pointer fPatch = buildPatchImage(fixed);
        PatchImage::Pointer mPatch = buildPatchImage(moving);

        PatchImageVector gradientPatchImages;
        buildGradientPatchImage(moving, 5, gradientPatchImages);

        DisplacementFieldType::Pointer flowOutput = io.NewImageS<RealImage>(fixed);
        DisplacementFieldType::PixelType zeroFlow;
        zeroFlow.Fill(0);
        flowOutput->FillBuffer(zeroFlow);

        PatchImage::PixelType* fBuf = fPatch->GetBufferPointer();
        PatchImage::PixelType* mBuf = mPatch->GetBufferPointer();
        PatchImage::PixelType* gBuf[RealImage::ImageDimension];
        fordim (k) {
            gBuf[k] = gradientPatchImages[k]->GetBufferPointer();
        }

        DisplacementFieldType::PixelType* oBuf = flowOutput->GetBufferPointer();

        const int nPatchElems = gBuf[0]->Size();
        const int nImageSize = fixed->GetPixelContainer()->Size();

        VNLDoubleMatrix A(__Dim, __Dim);
        VNLDoubleVector b(__Dim);
        for (int i = 0; i < nImageSize; i++) {
            A.fill(0);
            b.fill(0);
            for (int j = 0; j < nPatchElems; j++) {
                double gx = (*gBuf[0])[j];
                double gy = (*gBuf[1])[j];
                double It  = (*fBuf)[j] - (*mBuf)[j];
                A[0][0] += (gx*gx);
                A[0][1] += (gx*gy);
                A[1][1] += (gy*gy);
                b[0] -= gx*It;
                b[1] -= gy*It;
            }
            A[1][0] = A[0][1];

            if (abs(vnl_det(A[0], A[1])) > 1e-3) {
                VNLDoubleMatrix inv = vnl_matrix_inverse<double>(A);
                VNLDoubleVector sol = inv * b;
                fordim (k) {
                    (*oBuf)[k] = - (dt * sol[k]);
                }
            }

            fBuf ++;
            mBuf ++;
            oBuf ++;
            fordim (k) {
                gBuf[k]++;
            }
        }

        return flowOutput;
    }

    DisplacementFieldType::Pointer PatchCompare::computeDemonsFlow(RealImage::Pointer fixed, RealImage::Pointer moving, double dt) {
        ImageIO<DisplacementFieldType> io;

        // allocate output
        DisplacementFieldType::Pointer flowOutput = io.NewImageS<RealImage>(fixed);

        // compute gradient image
        ImageProcessing proc;
        RealImage::SpacingType spacing = fixed->GetSpacing();
        GradientImage::Pointer gF = proc.ComputeGaussianGradient(fixed, spacing[0]/2.0);
        GradientImage::Pointer gM = proc.ComputeGaussianGradient(moving, spacing[0]/2.0);

        // assign pointers
        RealImage::PixelType* fBuf = fixed->GetBufferPointer();
        RealImage::PixelType* mBuf = moving->GetBufferPointer();

        GradientImage::PixelType* gFv = gF->GetBufferPointer();
        GradientImage::PixelType* gMv = gM->GetBufferPointer();

        DisplacementFieldType::PixelType* oBuf = flowOutput->GetBufferPointer();
        const int nPixels = gF->GetPixelContainer()->Size();
        for (int i = 0; i < nPixels; i++) {
            double speed = (*fBuf - *mBuf);
            double gx = (*gMv)[0];
            double gy = (*gMv)[1];
            double gm = gx*gx + gy*gy;
            double denom = gm + speed*speed;
            if (denom > 1e-2) {
                (*oBuf)[0] = dt*speed/denom*gx;
                (*oBuf)[1] = dt*speed/denom*gy;
            }
            fBuf++;
            mBuf++;
            gFv++;
            gMv++;
            oBuf++;
        }

        return flowOutput;
    }



    void transferLabelsWithPatch(StringVector& args, std::string output, int searchRadius, int kNearest) {
        int nAtlas = (args.size() - 2) / 2;
        if (args.size() == 0 || args.size() != (nAtlas+1)*2 || nAtlas < 1) {
            cout << "LabelTransfer: target-image target-roi atlas-patch-#1 atlas-patch-#2 ... atlas-label-#1 atlas-label-#2 ..." << endl;
            return;
        }

        ImageIO<RealImage> realIO;
        ImageIO<LabelImage> labelIO;
        ImageIO<PatchImage> patchIO;

        RealImage::Pointer targetImage = realIO.ReadCastedImage(args[0]);
        LabelImage::Pointer targetROI = labelIO.ReadImage(args[1]);

        cout << "# of atlases: " << nAtlas << endl;
        PatchImageVector atlasPatches;
        for (int i = 2; i < nAtlas + 2; i++) {
            PatchImage::Pointer patchAtlas = patchIO.ReadImage(args[i]);
            cout << "Reading patch #" << (i-1) << endl;
            atlasPatches.push_back(patchAtlas);
        }

        LabelImageVector atlasLabels;
        for (int i = nAtlas + 2; i < nAtlas * 2 + 2; i++) {
            LabelImage::Pointer labelAtlas = labelIO.ReadImage(args[i]);
            cout << "Reading label #" << (i-nAtlas-1) << endl;
            atlasLabels.push_back(labelAtlas);
        }

        PatchCompare patchRun;
        patchRun.setAtlasImages(atlasPatches);
        patchRun.setAtlasLabels(atlasLabels);
        patchRun.setTargetImage(targetImage);
        patchRun.setTargetROI(targetROI);
        patchRun.estimateLabel(searchRadius, kNearest);

        labelIO.WriteImage(output, patchRun.getTargetLabel());
    }



    /**
     * ITK v4 framework
    int PatchCompare::performDeformableRegistration(PatchImage::Pointer fixedImage, PatchImage::Pointer movingImage, RealImage::Pointer movingSource) {
        typedef itk::TranslationTransform<double, __Dim> TransformType;
        typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
        typedef itk::MeanSquaresImageToImageMetric<PatchImage, PatchImage> MetricType;
        typedef itk::LinearInterpolateImageFunction<PatchImage, double> InterpolatorType;
        typedef itk::ImageRegistrationMethod<PatchImage, PatchImage> RegistrationType;
        typedef itk::ResampleImageFilter<RealImage, RealImage> ResampleFilterType;

        // Create components
        MetricType::Pointer         metric        = MetricType::New();
        TransformType::Pointer      transform     = TransformType::New();
        OptimizerType::Pointer      optimizer     = OptimizerType::New();
        InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
        RegistrationType::Pointer   registration  = RegistrationType::New();

        // Each component is now connected to the instance of the registration method.
        registration->SetMetric(        metric        );
        registration->SetOptimizer(     optimizer     );
        registration->SetTransform(     transform     );
        registration->SetInterpolator(  interpolator  );

        // Set the registration inputs
        registration->SetFixedImage(fixedImage);
        registration->SetMovingImage(movingImage);

        registration->SetFixedImageRegion(fixedImage->GetLargestPossibleRegion());

        //  Initialize the transform
        typedef RegistrationType::ParametersType ParametersType;
        ParametersType initialParameters( transform->GetNumberOfParameters() );

        initialParameters[0] = 0.0;  // Initial offset along X
        initialParameters[1] = 0.0;  // Initial offset along Y

        registration->SetInitialTransformParameters( initialParameters );

        optimizer->SetMaximumStepLength( 4.00 );
        optimizer->SetMinimumStepLength( 0.01 );

        // Set a stopping criterion
        optimizer->SetNumberOfIterations( 200 );

        // Connect an observer
        //CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
        //optimizer->AddObserver( itk::IterationEvent(), observer );

        try
        {
          registration->Update();
        }
        catch( itk::ExceptionObject & err )
        {
          std::cerr << "ExceptionObject caught !" << std::endl;
          std::cerr << err << std::endl;
          return EXIT_FAILURE;
        }

        //  The result of the registration process is an array of parameters that
        //  defines the spatial transformation in an unique way. This final result is
        //  obtained using the \code{GetLastTransformParameters()} method.

        ParametersType finalParameters = registration->GetLastTransformParameters();

        //  In the case of the \doxygen{TranslationTransform}, there is a
        //  straightforward interpretation of the parameters.  Each element of the
        //  array corresponds to a translation along one spatial dimension.

        const double TranslationAlongX = finalParameters[0];
        const double TranslationAlongY = finalParameters[1];

        //  The optimizer can be queried for the actual number of iterations
        //  performed to reach convergence.  The \code{GetCurrentIteration()}
        //  method returns this value. A large number of iterations may be an
        //  indication that the maximum step length has been set too small, which
        //  is undesirable since it results in long computational times.

        const unsigned int numberOfIterations = optimizer->GetCurrentIteration();

        //  The value of the image metric corresponding to the last set of parameters
        //  can be obtained with the \code{GetValue()} method of the optimizer.

        const double bestValue = optimizer->GetValue();

        // Print out results
        //
        std::cout << "Result = " << std::endl;
        std::cout << " Translation X = " << TranslationAlongX  << std::endl;
        std::cout << " Translation Y = " << TranslationAlongY  << std::endl;
        std::cout << " Iterations    = " << numberOfIterations << std::endl;
        std::cout << " Metric value  = " << bestValue          << std::endl;


        //  A resampling filter is created and the moving image is connected as  its input.
        ResampleFilterType::Pointer resampler = ResampleFilterType::New();
        resampler->SetInput(movingSource);
        resampler->SetTransform(registration->GetOutput()->Get());
        resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
        resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
        resampler->SetOutputSpacing( fixedImage->GetSpacing() );
        resampler->SetOutputDirection( fixedImage->GetDirection() );
        resampler->SetDefaultPixelValue( 100 );
        resampler->Update();

        pi::ImageIO<RealImage> realIO;
        realIO.WriteImage("/tmpfs/resample_output.nrrd", resampler->GetOutput());
    }
     */

    class PatchTracker {
    public:
        typedef itk::VectorLinearInterpolateImageFunction<PatchImage> PatchInterpolator;

        PatchTracker(PatchImage::Pointer f, PatchImage::Pointer m) {
            fixedImage = f;
            movingImage = m;
            fixedInterpolator = PatchInterpolator::New();
            movingInterpolator = PatchInterpolator::New();
            fixedInterpolator->SetInputImage(f);
            movingInterpolator->SetInputImage(m);
        }
        void SetMovingPoint(PatchImage::PointType point) {
            movingPoint = point;
            movingPixel = movingInterpolator->Evaluate(movingPoint);

        }
        double operator()(int n, double *x) {
            PatchImage::PointType point;
            point[0] = x[0];
            point[1] = x[1];
            PatchImage::PixelType f = fixedInterpolator->Evaluate(point);
            PatchImage::PixelType::ValueType* fBuf = f.GetDataPointer();
            PatchImage::PixelType::ValueType* mBuf = movingPixel.GetDataPointer();
            const int size = f.Size();
            double ssd = 0;
            for (int i = 0; i < size; i++) {
                ssd += (*fBuf - *mBuf)*(*fBuf - *mBuf);
                fBuf++;
                mBuf++;
            }
            ssd = std::sqrt(ssd);
            cout << x[0] << ", " << x[1] << ", " << ssd << endl;
            return ssd;
        }

    private:
        PatchImage::PointType movingPoint;
        PatchImage::PixelType movingPixel;
        PatchImage::Pointer fixedImage;
        PatchImage::Pointer movingImage;
        PatchInterpolator::Pointer fixedInterpolator;
        PatchInterpolator::Pointer movingInterpolator;
    };

    DisplacementFieldType::Pointer PatchCompare::performDenseMapping(PatchImage::Pointer fixed, PatchImage::Pointer moving, PatchImage::RegionType activeRegion) {

        typedef itk::ImageRegionConstIteratorWithIndex<PatchImage> PatchIteratorType;

        PatchImage::PointType fixedPoint, movingPoint;
        PatchIteratorType movingIter(moving, activeRegion);

        DisplacementFieldType::Pointer deformField = DisplacementFieldType::New();
        DisplacementFieldType::PixelType zeroDisplacement;
        zeroDisplacement.Fill(0);
        deformField->SetSpacing(moving->GetSpacing());
        deformField->SetRegions(moving->GetBufferedRegion());
        deformField->SetOrigin(moving->GetOrigin());
        deformField->SetDirection(moving->GetDirection());
        deformField->Allocate();
        deformField->FillBuffer(zeroDisplacement);

        PatchTracker tracker(fixed, moving);
        for (movingIter.GoToBegin(); !movingIter.IsAtEnd(); ++movingIter) {
            // compute index point
            PatchImage::IndexType idx = movingIter.GetIndex();
            fixed->TransformIndexToPhysicalPoint(idx, fixedPoint);
            moving->TransformIndexToPhysicalPoint(idx, movingPoint);

            // from this point, compute closest point
            PowellOpti<PatchTracker> optimizer;

            PowellParams initial(2);
            initial[0] = fixedPoint[0];
            initial[1] = fixedPoint[1];

            tracker.SetMovingPoint(movingPoint);
            optimizer.minimizeNEWUOA(tracker, initial, 1);

            // set as displacement
            DisplacementFieldType::PixelType displacement;
            fordim (k) {
                displacement[k] = initial[k] - movingPoint[k];
            }
            cout << displacement[0] << "," << displacement[1] << endl;
            deformField->SetPixel(idx, displacement);
        }

        return deformField;
    }

    int PatchCompare::performDemonsRegistration(RealImage::Pointer fixedImage, RealImage::Pointer movingImage, std::string outputImage) {

        typedef itk::DemonsRegistrationFilter<RealImage, RealImage, DisplacementFieldType> DemonsFilterType;

        DemonsFilterType::Pointer demonsFilter = DemonsFilterType::New();
        demonsFilter->SetFixedImage(fixedImage);
        demonsFilter->SetMovingImage(movingImage);
        demonsFilter->SetNumberOfIterations(1000);
        demonsFilter->SmoothDisplacementFieldOn();
        demonsFilter->SmoothUpdateFieldOn();
//        demonsFilter->SetUpdateFieldStandardDeviations(0.1);
//        demonsFilter->SetStandardDeviations(0.1);
        demonsFilter->Update();


        typedef itk::DisplacementFieldTransform<double, __Dim> DisplacementTransformType;
        DisplacementTransformType::Pointer displacementTransform = DisplacementTransformType::New();
        displacementTransform->SetDisplacementField(demonsFilter->GetDisplacementField());

        typedef itk::ResampleImageFilter<RealImage, RealImage> ResampleFilterType;
        ResampleFilterType::Pointer resampler = ResampleFilterType::New();
        resampler->SetInput(movingImage);
        resampler->SetTransform(displacementTransform);
        resampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
        resampler->SetOutputOrigin(fixedImage->GetOrigin());
        resampler->SetOutputSpacing(fixedImage->GetSpacing());
        resampler->SetOutputDirection(fixedImage->GetDirection());
        resampler->SetDefaultPixelValue(0);
        resampler->Update();

        ImageIO<RealImage> io;
        io.WriteImage(outputImage, resampler->GetOutput());

        return EXIT_SUCCESS;
    }

#pragma mark -
    /**
     * Constructor
     */
    template< class TInputImage, class TCoordRep, class TOutputType >
    CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >
    ::CentralDifferenceImageFunction()
    {
        this->m_UseImageDirection = true;

        /* Interpolator. Default to linear. */
        typedef itk::LinearInterpolateImageFunction< TInputImage, TCoordRep >
        LinearInterpolatorType;
        this->m_Interpolator = LinearInterpolatorType::New();
    }

    /**
     *
     */
    template< class TInputImage, class TCoordRep, class TOutputType >
    void
    CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >
    ::SetInputImage(const TInputImage *inputData)
    {
        if ( inputData != this->m_Image )
        {
            Superclass::SetInputImage( inputData );
            this->m_Interpolator->SetInputImage( inputData );

            // Verify the output vector is the right size.
            // OutputType of VariablelengthVector will have size 0 until allocated, so this
            // case can't be tested.
            if( inputData != NULL )
            {
                itk::SizeValueType nComponents = OutputConvertType::GetNumberOfComponents();
                if( nComponents > 0 )
                {
                    if( nComponents != inputData->GetNumberOfComponentsPerPixel() * TInputImage::ImageDimension )
                    {
                        itkExceptionMacro("The OutputType is not the right size (" << nComponents << ") for the given pixel size ("
                                          << inputData->GetNumberOfComponentsPerPixel() << ") and image dimension (" << TInputImage::ImageDimension << ").")
                    }
                }
            }
            this->Modified();
        }
    }

    /**
     *
     */
    template< class TInputImage, class TCoordRep, class TOutputType >
    void
    CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >
    ::SetInterpolator(InterpolatorType *interpolator )
    {
        if ( interpolator != this->m_Interpolator )
        {
            this->m_Interpolator = interpolator;
            if( this->GetInputImage() != NULL )
            {
                this->m_Interpolator->SetInputImage( this->GetInputImage() );
            }
            this->Modified();
        }
    }

    /**
     *
     */
    template< class TInputImage, class TCoordRep, class TOutputType >
    void
    CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >
    ::PrintSelf(std::ostream & os, itk::Indent indent) const
    {
        this->Superclass::PrintSelf(os, indent);
        os << indent << "UseImageDirection = " << this->m_UseImageDirection << std::endl;
    }

    /**
     * EvaluateAtIndex
     */
    template< class TInputImage, class TCoordRep, class TOutputType >
    typename CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >::OutputType
    CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >
    ::EvaluateAtIndex(const IndexType & index) const
    {
        OutputType derivative;

        // When ScalarDerivativeType is the same as OutputType, this calls
        // the version specialized for scalar pixels since in that case,
        // the two vector types are the same.
        EvaluateAtIndexSpecialized<ScalarDerivativeType>( index, derivative, OutputTypeSpecializationStructType<ScalarDerivativeType>() );

        return derivative;
    }

    /*
     * Specialized for scalar pixels
     */
    template< class TInputImage, class TCoordRep, class TOutputType >
    template< class Type >
    void
    CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >
    ::EvaluateAtIndexSpecialized(const IndexType & index, OutputType & orientedDerivative, OutputTypeSpecializationStructType<OutputType>) const
    {
        OutputType derivative;

        IndexType neighIndex = index;

        const InputImageType *inputImage = this->GetInputImage();

        const typename InputImageType::RegionType & region =
        inputImage->GetBufferedRegion();

        const typename InputImageType::SizeType & size   = region.GetSize();
        const typename InputImageType::IndexType & start = region.GetIndex();

        const unsigned int MaxDims = Self::ImageDimension;
        for ( unsigned int dim = 0; dim < MaxDims; dim++ )
        {
            // bounds checking
            // checks for index either on the boundary or out of bounds.
            // note that the documentation says this method assumes the index
            // is in-bounds, so we don't do anything else if the point is out of bounds.
            if ( index[dim] < start[dim] + 1 || index[dim] > ( start[dim] + static_cast< itk::OffsetValueType >( size[dim] ) - 2 ) )
            {
                derivative[dim] = itk::NumericTraits<OutputValueType>::Zero;
                continue;
            }

            // compute derivative
            neighIndex[dim] += 1;
            derivative[dim] = inputImage->GetPixel(neighIndex);

            neighIndex[dim] -= 2;
            derivative[dim] -= inputImage->GetPixel(neighIndex);

            derivative[dim] *= static_cast<OutputValueType>(0.5) / inputImage->GetSpacing()[dim];
            neighIndex[dim] += 1;
        }

        if ( this->m_UseImageDirection )
        {
            inputImage->TransformLocalVectorToPhysicalVector(derivative, orientedDerivative);
        }
        else
        {
            orientedDerivative = derivative;
        }

    }

    /*
     * Specialized for vector pixels
     */
    template< class TInputImage, class TCoordRep, class TOutputType >
    template< class Type >
    void
    CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >
    ::EvaluateAtIndexSpecialized(const IndexType & index, OutputType & derivative, OutputTypeSpecializationStructType<Type>) const
    {
        const InputImageType *inputImage = this->GetInputImage();
        const unsigned int numberComponents = this->GetInputImage()->GetNumberOfComponentsPerPixel();

        IndexType neighIndex = index;

        const typename InputImageType::RegionType & region = inputImage->GetBufferedRegion();
        const typename InputImageType::SizeType & size     = region.GetSize();
        const typename InputImageType::IndexType & start   = region.GetIndex();

        typedef typename InputImageType::PixelType PixelType;
        const PixelType * neighPixels[Self::ImageDimension][2];
        const PixelType zeroPixel = itk::NumericTraits<PixelType>::ZeroValue();
        const unsigned int MaxDims = Self::ImageDimension;
        bool  dimOutOfBounds[Self::ImageDimension];

        for ( unsigned int dim = 0; dim < MaxDims; dim++ )
        {
            // initialize to quiet compiler warnings
            neighPixels[dim][0] = &zeroPixel;
            neighPixels[dim][1] = &zeroPixel;

            // cached bounds checking
            dimOutOfBounds[dim] = ( ( index[dim] < (start[dim] + 1) ) || index[dim] > ( start[dim] + static_cast< itk::OffsetValueType >( size[dim] ) - 2 ) );
        }

        for ( unsigned int nc = 0; nc < numberComponents; nc++)
        {
            ScalarDerivativeType componentDerivative;

            for ( unsigned int dim = 0; dim < MaxDims; dim++ )
            {
                // bounds checking
                if( dimOutOfBounds[dim] )
                {
                    componentDerivative[dim] = itk::NumericTraits<OutputValueType>::ZeroValue();
                    continue;
                }

                // get pixels
                if( nc == 0 )
                {
                    neighIndex[dim] += 1;
                    neighPixels[dim][0] = &( inputImage->GetPixel(neighIndex) );
                    neighIndex[dim] -= 2;
                    neighPixels[dim][1] = &( inputImage->GetPixel(neighIndex) );
                    neighIndex[dim] += 1;
                }

                // compute derivative
                componentDerivative[dim] = InputPixelConvertType::GetNthComponent( nc, *neighPixels[dim][0] );
                componentDerivative[dim] -= InputPixelConvertType::GetNthComponent( nc, *neighPixels[dim][1] );
                componentDerivative[dim] *= 0.5 / inputImage->GetSpacing()[dim];
            }

            if ( this->m_UseImageDirection )
            {
                ScalarDerivativeType componentDerivativeOut;
                inputImage->TransformLocalVectorToPhysicalVector(componentDerivative, componentDerivativeOut);
                for ( unsigned int dim = 0; dim < MaxDims; dim++ )
                {
                    OutputConvertType::SetNthComponent( nc * MaxDims + dim, derivative, componentDerivativeOut[dim] );
                }
            }
            else
            {
                for ( unsigned int dim = 0; dim < MaxDims; dim++ )
                {
                    OutputConvertType::SetNthComponent( nc * MaxDims + dim, derivative, componentDerivative[dim] );
                }
            }
        }
    }

    /**
     *
     */
    template< class TInputImage, class TCoordRep, class TOutputType >
    typename CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >::OutputType
    CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >
    ::Evaluate(const PointType & point) const
    {
        OutputType derivative;

        // When ScalarDerivativeType is the same as OutputType, this calls
        // the version specialized for scalar pixels since in that case,
        // the two vector types are the same.
        EvaluateSpecialized<ScalarDerivativeType>( point, derivative, OutputTypeSpecializationStructType<ScalarDerivativeType>() );

        return derivative;
    }

    /*
     * Specialized for scalar pixels
     */
    template< class TInputImage, class TCoordRep, class TOutputType >
    template< class Type >
    void
    CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >
    ::EvaluateSpecialized(const PointType & point, OutputType & orientedDerivative, OutputTypeSpecializationStructType<OutputType>) const
    {
        typedef typename PointType::ValueType           PointValueType;
        typedef typename OutputType::ValueType          DerivativeValueType;
        typedef typename ContinuousIndexType::ValueType ContinuousIndexValueType;

        PointType neighPoint1 = point;
        PointType neighPoint2 = point;

        const InputImageType *inputImage = this->GetInputImage();

        const SpacingType & spacing = inputImage->GetSpacing();

        const unsigned int MaxDims = Self::ImageDimension;
        for ( unsigned int dim = 0; dim < MaxDims; dim++ )
        {
            PointValueType offset = static_cast<PointValueType>(0.5) * spacing[dim];
            // Check the bounds using the point because the image direction may swap dimensions,
            // making checks in index space inaccurate.
            // If on a boundary, we set the derivative to zero. This is done to match the behavior
            // of EvaluateAtIndex. Another approach is to calculate the 1-sided difference.
            neighPoint1[dim] = point[dim] - offset;
            if( ! this->IsInsideBuffer( neighPoint1 ) )
            {
                orientedDerivative[dim] = itk::NumericTraits<DerivativeValueType>::Zero;
                neighPoint1[dim] = point[dim];
                neighPoint2[dim] = point[dim];
                continue;
            }
            neighPoint2[dim] = point[dim] + offset;
            if( ! this->IsInsideBuffer( neighPoint2 ) )
            {
                orientedDerivative[dim] = itk::NumericTraits<DerivativeValueType>::Zero;
                neighPoint1[dim] = point[dim];
                neighPoint2[dim] = point[dim];
                continue;
            }

            PointValueType delta = neighPoint2[dim] - neighPoint1[dim];
            if( delta > 10.0 * itk::NumericTraits<PointValueType>::epsilon() )
            {
                orientedDerivative[dim] = ( this->m_Interpolator->Evaluate( neighPoint2 ) - this->m_Interpolator->Evaluate( neighPoint1 ) ) / delta;
            }
            else
            {
                orientedDerivative[dim] = static_cast<DerivativeValueType>(0.0);
            }

            neighPoint1[dim] = point[dim];
            neighPoint2[dim] = point[dim];
        }

        // Since we've implicitly calculated the derivative with respect to image
        // direction, we need to reorient into index-space if the user desires.
        if ( ! this->m_UseImageDirection )
        {
            OutputType derivative;
            inputImage->TransformPhysicalVectorToLocalVector( orientedDerivative, derivative );
            orientedDerivative = derivative;
        }
    }

    /*
     * Specialized for vector pixels
     */
    template< class TInputImage, class TCoordRep, class TOutputType >
    template< class Type >
    void
    CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >
    ::EvaluateSpecialized(const PointType & point, OutputType & derivative, OutputTypeSpecializationStructType<Type>) const
    {
        typedef typename PointType::ValueType           PointValueType;
        typedef typename OutputType::ValueType          DerivativeValueType;
        typedef typename ContinuousIndexType::ValueType ContinuousIndexValueType;

        const InputImageType *inputImage = this->GetInputImage();
        const unsigned int numberComponents = inputImage->GetNumberOfComponentsPerPixel();

        PointType neighPoint1 = point;
        PointType neighPoint2 = point;

        const SpacingType & spacing = inputImage->GetSpacing();

        typedef typename InputImageType::PixelType PixelType;
        PixelType neighPixels[Self::ImageDimension][2];
        bool  dimOutOfBounds[Self::ImageDimension];
        const unsigned int MaxDims = Self::ImageDimension;
        PointValueType delta[Self::ImageDimension];
        PixelType zeroPixel = itk::NumericTraits<PixelType>::ZeroValue();

        ScalarDerivativeType componentDerivativeOut;
        ScalarDerivativeType componentDerivative;
        componentDerivative.Fill( itk::NumericTraits<OutputValueType>::Zero );

        for ( unsigned int dim = 0; dim < Self::ImageDimension; dim++ )
        {
            // initialize to quiet compiler warnings
            neighPixels[dim][0] = zeroPixel;
            neighPixels[dim][1] = zeroPixel;
            delta[dim] = itk::NumericTraits<PointValueType>::ZeroValue();
            dimOutOfBounds[dim] = true;
        }

        for ( unsigned int nc = 0; nc < numberComponents; nc++ )
        {
            for ( unsigned int dim = 0; dim < MaxDims; dim++ )
            {
                // Initialize values that only depend on dimension and not component number.
                if( nc == 0 )
                {
                    // Check the bounds using the point because the image direction may swap dimensions,
                    // making checks in index space inaccurate.
                    // If on a boundary, we set the derivative to zero. This is done to match the behavior
                    // of EvaluateAtIndex. Another approach is to calculate the 1-sided difference.
                    PointValueType offset = static_cast<PointValueType>(0.5) * spacing[dim];
                    neighPoint1[dim] = point[dim] - offset;
                    neighPoint2[dim] = point[dim] + offset;
                    dimOutOfBounds[dim] = ( ! this->IsInsideBuffer( neighPoint1 ) || ! this->IsInsideBuffer( neighPoint2 ) );

                    if( dimOutOfBounds[dim] )
                    {
                        componentDerivative[dim] = itk::NumericTraits<OutputValueType>::Zero;
                        neighPoint1[dim] = point[dim];
                        neighPoint2[dim] = point[dim];
                        continue;
                    }

                    neighPixels[dim][0] = this->m_Interpolator->Evaluate( neighPoint2 );
                    neighPixels[dim][1] = this->m_Interpolator->Evaluate( neighPoint1 );

                    delta[dim] = neighPoint2[dim] - neighPoint1[dim];

                    neighPoint1[dim] = point[dim];
                    neighPoint2[dim] = point[dim];
                }
                else
                {
                    if( dimOutOfBounds[dim] )
                    {
                        continue;
                    }
                }

                if( delta[dim] > 10.0 * itk::NumericTraits<PointValueType>::epsilon() )
                {
                    OutputValueType left = InputPixelConvertType::GetNthComponent( nc, neighPixels[dim][0] );
                    OutputValueType right = InputPixelConvertType::GetNthComponent( nc, neighPixels[dim][1] );
                    componentDerivative[dim] = (left - right) / delta[dim];
                }
                else
                {
                    componentDerivative[dim] = itk::NumericTraits<OutputValueType>::ZeroValue();
                }
            }

            // Since we've implicitly calculated the derivative with respect to image
            // direction, we need to reorient into index-space if the user
            // desires.
            if ( ! this->m_UseImageDirection )
            {
                inputImage->TransformPhysicalVectorToLocalVector(componentDerivative, componentDerivativeOut);
                for ( unsigned int dim = 0; dim < MaxDims; dim++ )
                {
                    OutputConvertType::SetNthComponent( nc * MaxDims + dim, derivative, componentDerivativeOut[dim] );
                }
            }
            else
            {
                for ( unsigned int dim = 0; dim < MaxDims; dim++ )
                {
                    OutputConvertType::SetNthComponent( nc * MaxDims + dim, derivative, componentDerivative[dim] );
                }
            }
        }
    }

    /**
     * EvaluateAtContinuousIndex
     */
    template< class TInputImage, class TCoordRep, class TOutputType >
    typename CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >::OutputType
    CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >
    ::EvaluateAtContinuousIndex(const ContinuousIndexType & cindex) const
    {
        OutputType derivative;
        // When ScalarDerivativeType is the same as OutputType, this calls
        // the version specialized for scalar pixels since in that case,
        // the two vector types are the same.
        this->EvaluateAtContinuousIndexSpecialized<ScalarDerivativeType>( cindex, derivative, OutputTypeSpecializationStructType<ScalarDerivativeType>() );
        return derivative;
    }

    /*
     * Specialized for scalar pixels
     */
    template< class TInputImage, class TCoordRep, class TOutputType >
    template< class Type >
    void
    CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >
    ::EvaluateAtContinuousIndexSpecialized(const ContinuousIndexType & cindex, OutputType & orientedDerivative, OutputTypeSpecializationStructType<OutputType>) const
    {
        typedef typename OutputType::ValueType          DerivativeValueType;
        typedef typename ContinuousIndexType::ValueType ContinuousIndexValueType;

        OutputType derivative;

        ContinuousIndexType neighIndex = cindex;

        const InputImageType *inputImage = this->GetInputImage();

        const typename InputImageType::RegionType & region = inputImage->GetBufferedRegion();

        const typename InputImageType::SizeType & size   = region.GetSize();
        const typename InputImageType::IndexType & start = region.GetIndex();

        const unsigned int MaxDims = Self::ImageDimension;
        for ( unsigned int dim = 0; dim < MaxDims; dim++ )
        {
            // bounds checking
            if ( cindex[dim] < static_cast<ContinuousIndexValueType>(start[dim] + 1)
                || cindex[dim] > static_cast<ContinuousIndexValueType>
                ( start[dim] + static_cast< itk::OffsetValueType >( size[dim] ) - 2 ) )
            {
                derivative[dim] = itk::NumericTraits<DerivativeValueType>::Zero;
                continue;
            }

            // compute derivative
            neighIndex[dim] += static_cast<ContinuousIndexValueType>(1.0);
            derivative[dim] = this->m_Interpolator->EvaluateAtContinuousIndex(neighIndex);

            neighIndex[dim] -= static_cast<ContinuousIndexValueType>(2.0);
            derivative[dim] -= this->m_Interpolator->EvaluateAtContinuousIndex(neighIndex);

            derivative[dim] *= static_cast<ContinuousIndexValueType>(0.5) / inputImage->GetSpacing()[dim];
            neighIndex[dim] += static_cast<ContinuousIndexValueType>(1.0);
        }

        if ( this->m_UseImageDirection )
        {
            inputImage->TransformLocalVectorToPhysicalVector(derivative, orientedDerivative);
        }
        else
        {
            orientedDerivative = derivative;
        }
    }
    
    /*
     * Specialized for vector pixels
     */
    template< class TInputImage, class TCoordRep, class TOutputType >
    template< class Type >
    void
    CentralDifferenceImageFunction< TInputImage, TCoordRep, TOutputType >
    ::EvaluateAtContinuousIndexSpecialized(const ContinuousIndexType & cindex, OutputType & derivative, OutputTypeSpecializationStructType<Type>) const
    {
        typedef typename OutputType::ValueType          DerivativeValueType;
        typedef typename ContinuousIndexType::ValueType ContinuousIndexValueType;
        
        const InputImageType *inputImage = this->GetInputImage();
        const unsigned int numberComponents = inputImage->GetNumberOfComponentsPerPixel();
        
        ContinuousIndexType neighIndex = cindex;
        const typename InputImageType::RegionType & region = inputImage->GetBufferedRegion();
        
        const typename InputImageType::SizeType & size   = region.GetSize();
        const typename InputImageType::IndexType & start = region.GetIndex();
        
        typedef typename InputImageType::PixelType PixelType;
        PixelType neighPixels[Self::ImageDimension][2];
        bool  dimOutOfBounds[Self::ImageDimension];
        const unsigned int MaxDims = Self::ImageDimension;
        PixelType zeroPixel = itk::NumericTraits<PixelType>::ZeroValue();
        
        for ( unsigned int dim = 0; dim < MaxDims; dim++ )
        {
            // initialize to quiet compiler warnings
            neighPixels[dim][0] = zeroPixel;
            neighPixels[dim][1] = zeroPixel;
            
            // bounds checking
            dimOutOfBounds[dim] = ( ( cindex[dim] < static_cast<ContinuousIndexValueType>(start[dim] + 1) )
                                   || cindex[dim] > static_cast<ContinuousIndexValueType> ( start[dim] + static_cast< itk::OffsetValueType >( size[dim] ) - 2 ) );
        }
        
        for ( unsigned int nc = 0; nc < numberComponents; nc++)
        {
            ScalarDerivativeType componentDerivative;
            ScalarDerivativeType componentDerivativeOut;
            
            for ( unsigned int dim = 0; dim < MaxDims; dim++ )
            {
                if( dimOutOfBounds[dim] )
                {
                    componentDerivative[dim] = itk::NumericTraits<DerivativeValueType>::ZeroValue();
                    continue;
                }
                
                // get pixels
                if( nc == 0 )
                {
                    neighIndex[dim] += static_cast<ContinuousIndexValueType>(1.0);
                    neighPixels[dim][0] = this->m_Interpolator->EvaluateAtContinuousIndex(neighIndex);
                    neighIndex[dim] -= static_cast<ContinuousIndexValueType>(2.0);
                    neighPixels[dim][1] = this->m_Interpolator->EvaluateAtContinuousIndex(neighIndex);
                    neighIndex[dim] += static_cast<ContinuousIndexValueType>(1.0);
                }
                
                // compute derivative
                componentDerivative[dim] = InputPixelConvertType::GetNthComponent(nc, neighPixels[dim][0] );
                componentDerivative[dim] -= InputPixelConvertType::GetNthComponent(nc, neighPixels[dim][1] );
                componentDerivative[dim] *= static_cast<ContinuousIndexValueType>(0.5) / inputImage->GetSpacing()[dim];
            }
            
            if ( this->m_UseImageDirection )
            {
                inputImage->TransformLocalVectorToPhysicalVector(componentDerivative, componentDerivativeOut);
                for ( unsigned int dim = 0; dim < MaxDims; dim++ )
                {
                    OutputConvertType::SetNthComponent( nc * MaxDims + dim, derivative, componentDerivativeOut[dim] );
                }
            }
            else
            {
                for ( unsigned int dim = 0; dim < MaxDims; dim++ )
                {
                    OutputConvertType::SetNthComponent( nc * MaxDims + dim, derivative, componentDerivative[dim] );
                }
            }
        }
    }

#pragma mark -
    /**
     * Default constructor
     */
    template< class TFixedImage, class TMovingImage, class TDisplacementField >
    DemonsRegistrationFunction< TFixedImage, TMovingImage, TDisplacementField >
    ::DemonsRegistrationFunction()
    {
        RadiusType   r;
        unsigned int j;

        for ( j = 0; j < ImageDimension; j++ )
        {
            r[j] = 0;
        }
        this->SetRadius(r);

        m_TimeStep = 1.0;
        m_DenominatorThreshold = 1e-9;
        m_IntensityDifferenceThreshold = 0.001;
        this->SetMovingImage(NULL);
        this->SetFixedImage(NULL);
        //m_FixedImageSpacing.Fill( 1.0 );
        //m_FixedImageOrigin.Fill( 0.0 );
        m_Normalizer = 1.0;
        m_FixedImageGradientCalculator = GradientCalculatorType::New();

        typename DefaultInterpolatorType::Pointer interp =
        DefaultInterpolatorType::New();

        m_MovingImageInterpolator = static_cast< InterpolatorType * >(
                                                                      interp.GetPointer() );

        m_Metric = itk::NumericTraits< double >::max();
        m_SumOfSquaredDifference = 0.0;
        m_NumberOfPixelsProcessed = 0L;
        m_RMSChange = itk::NumericTraits< double >::max();
        m_SumOfSquaredChange = 0.0;

        m_MovingImageGradientCalculator = MovingImageGradientCalculatorType::New();
        m_UseMovingImageGradient = false;
    }

    /**
     * Standard "PrintSelf" method.
     */
    template< class TFixedImage, class TMovingImage, class TDisplacementField >
    void
    DemonsRegistrationFunction< TFixedImage, TMovingImage, TDisplacementField >
    ::PrintSelf(std::ostream & os, itk::Indent indent) const
    {
        Superclass::PrintSelf(os, indent);

        os << indent << "MovingImageIterpolator: ";
        os << m_MovingImageInterpolator.GetPointer() << std::endl;
        os << indent << "FixedImageGradientCalculator: ";
        os << m_FixedImageGradientCalculator.GetPointer() << std::endl;
        os << indent << "DenominatorThreshold: ";
        os << m_DenominatorThreshold << std::endl;
        os << indent << "IntensityDifferenceThreshold: ";
        os << m_IntensityDifferenceThreshold << std::endl;

        os << indent << "UseMovingImageGradient: ";
        os << m_UseMovingImageGradient << std::endl;

        os << indent << "Metric: ";
        os << m_Metric << std::endl;
        os << indent << "SumOfSquaredDifference: ";
        os << m_SumOfSquaredDifference << std::endl;
        os << indent << "NumberOfPixelsProcessed: ";
        os << m_NumberOfPixelsProcessed << std::endl;
        os << indent << "RMSChange: ";
        os << m_RMSChange << std::endl;
        os << indent << "SumOfSquaredChange: ";
        os << m_SumOfSquaredChange << std::endl;
    }

    /**
     *
     */
    template< class TFixedImage, class TMovingImage, class TDisplacementField >
    void
    DemonsRegistrationFunction< TFixedImage, TMovingImage, TDisplacementField >
    ::SetIntensityDifferenceThreshold(double threshold)
    {
        m_IntensityDifferenceThreshold = threshold;
    }

    /**
     *
     */
    template< class TFixedImage, class TMovingImage, class TDisplacementField >
    double
    DemonsRegistrationFunction< TFixedImage, TMovingImage, TDisplacementField >
    ::GetIntensityDifferenceThreshold() const
    {
        return m_IntensityDifferenceThreshold;
    }

    /**
     * Set the function state values before each iteration
     */
    template< class TFixedImage, class TMovingImage, class TDisplacementField >
    void
    DemonsRegistrationFunction< TFixedImage, TMovingImage, TDisplacementField >
    ::InitializeIteration()
    {
        if ( !this->GetMovingImage() || !this->GetFixedImage() || !m_MovingImageInterpolator )
        {
            itkExceptionMacro(<< "MovingImage, FixedImage and/or Interpolator not set");
        }

        // cache fixed image information
        SpacingType fixedImageSpacing    = this->GetFixedImage()->GetSpacing();
        m_ZeroUpdateReturn.Fill(0.0);

        // compute the normalizer
        m_Normalizer      = 0.0;
        for ( unsigned int k = 0; k < ImageDimension; k++ )
        {
            m_Normalizer += fixedImageSpacing[k] * fixedImageSpacing[k];
        }
        m_Normalizer /= static_cast< double >( ImageDimension );

        // setup gradient calculator
        m_FixedImageGradientCalculator->SetInputImage( this->GetFixedImage() );
        m_MovingImageGradientCalculator->SetInputImage( this->GetMovingImage() );

        // setup moving image interpolator
        m_MovingImageInterpolator->SetInputImage( this->GetMovingImage() );

        // initialize metric computation variables
        m_SumOfSquaredDifference  = 0.0;
        m_NumberOfPixelsProcessed = 0L;
        m_SumOfSquaredChange      = 0.0;
    }

    /**
     * Compute update at a specify neighbourhood
     */
    template< class TFixedImage, class TMovingImage, class TDisplacementField >
    typename DemonsRegistrationFunction< TFixedImage, TMovingImage, TDisplacementField >
    ::PixelType
    DemonsRegistrationFunction< TFixedImage, TMovingImage, TDisplacementField >
    ::ComputeUpdate( const NeighborhoodType & it, void *gd,
                    const FloatOffsetType & itkNotUsed(offset) )
    {
        // Get fixed image related information
        // Note: no need to check the index is within
        // fixed image buffer. This is done by the external filter.
        const IndexType index = it.GetIndex();
        const typename TFixedImage::PixelType fixedValue = this->GetFixedImage()->GetPixel(index);

        // Get moving image related information
        PointType mappedPoint;

        this->GetFixedImage()->TransformIndexToPhysicalPoint(index, mappedPoint);
        for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
            mappedPoint[j] += it.GetCenterPixel()[j];
        }

        typename TMovingImage::PixelType movingValue;
        if ( m_MovingImageInterpolator->IsInsideBuffer(mappedPoint) )
        {
            movingValue = m_MovingImageInterpolator->Evaluate(mappedPoint);
        }
        else
        {
            return m_ZeroUpdateReturn;
        }

        CovariantVectorType gradient;
        // Compute the gradient of either fixed or moving image
        if ( !m_UseMovingImageGradient )
        {
            gradient = m_FixedImageGradientCalculator->EvaluateAtIndex(index);
        }
        else
        {
            gradient = m_MovingImageGradientCalculator->Evaluate(mappedPoint);
        }

        double gradientSquaredMagnitude = 0;
        for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
            gradientSquaredMagnitude += vnl_math_sqr(gradient[j]);
        }

        /**
         * Compute Update.
         * In the original equation the denominator is defined as (g-f)^2 + grad_mag^2.
         * However there is a mismatch in units between the two terms.
         * The units for the second term is intensity^2/mm^2 while the
         * units for the first term is intensity^2. This mismatch is particularly
         * problematic when the fixed image does not have unit spacing.
         * In this implemenation, we normalize the first term by a factor K,
         * such that denominator = (g-f)^2/K + grad_mag^2
         * where K = mean square spacing to compensate for the mismatch in units.
         */

        typename TMovingImage::PixelType diffValue = fixedValue - movingValue;

        // estimate speedValue from patch difference norm
        const double speedValue = diffValue.GetNorm();
        const double sqr_speedValue = vnl_math_sqr(speedValue);

        // update the metric
        GlobalDataStruct *globalData = (GlobalDataStruct *)gd;
        if ( globalData )
        {
            globalData->m_SumOfSquaredDifference += sqr_speedValue;
            globalData->m_NumberOfPixelsProcessed += 1;
        }

        const double denominator = sqr_speedValue / m_Normalizer
        + gradientSquaredMagnitude;

        if ( vnl_math_abs(speedValue) < m_IntensityDifferenceThreshold
            || denominator < m_DenominatorThreshold )
        {
            return m_ZeroUpdateReturn;
        }
        
        PixelType update;
        for ( unsigned int j = 0; j < ImageDimension; j++ )
        {
            update[j] = speedValue * gradient[j] / denominator;
            if ( globalData )
            {
                globalData->m_SumOfSquaredChange += vnl_math_sqr(update[j]);
            }
        }
        return update;
    }
    
    /**
     * Update the metric and release the per-thread-global data.
     */
    template< class TFixedImage, class TMovingImage, class TDisplacementField >
    void
    DemonsRegistrationFunction< TFixedImage, TMovingImage, TDisplacementField >
    ::ReleaseGlobalDataPointer(void *gd) const
    {
        GlobalDataStruct *globalData = (GlobalDataStruct *)gd;
        
        m_MetricCalculationLock.Lock();
        m_SumOfSquaredDifference += globalData->m_SumOfSquaredDifference;
        m_NumberOfPixelsProcessed += globalData->m_NumberOfPixelsProcessed;
        m_SumOfSquaredChange += globalData->m_SumOfSquaredChange;
        if ( m_NumberOfPixelsProcessed )
        {
            m_Metric = m_SumOfSquaredDifference
            / static_cast< double >( m_NumberOfPixelsProcessed );
            m_RMSChange = vcl_sqrt( m_SumOfSquaredChange
                                   / static_cast< double >( m_NumberOfPixelsProcessed ) );
        }
        m_MetricCalculationLock.Unlock();
        
        delete globalData;
    }

#pragma mark -

    /**
     * Default constructor
     */
    template< class TFixedImage, class TMovingImage, class TDisplacementField >
    DemonsRegistrationFilter< TFixedImage, TMovingImage, TDisplacementField >
    ::DemonsRegistrationFilter()
    {
        typename DemonsRegistrationFunctionType::Pointer drfp;
        drfp = DemonsRegistrationFunctionType::New();

        this->SetDifferenceFunction( static_cast< FiniteDifferenceFunctionType * >(
                                                                                   drfp.GetPointer() ) );

        m_UseMovingImageGradient = false;
    }

    template< class TFixedImage, class TMovingImage, class TDisplacementField >
    void
    DemonsRegistrationFilter< TFixedImage, TMovingImage, TDisplacementField >
    ::PrintSelf(std::ostream & os, itk::Indent indent) const
    {
        Superclass::PrintSelf(os, indent);
        os << indent << "UseMovingImageGradient: ";
        os << m_UseMovingImageGradient << std::endl;
        os << indent << "Intensity difference threshold: "
        << this->GetIntensityDifferenceThreshold() << std::endl;
    }

    /*
     * Set the function state values before each iteration
     */
    template< class TFixedImage, class TMovingImage, class TDisplacementField >
    void
    DemonsRegistrationFilter< TFixedImage, TMovingImage, TDisplacementField >
    ::InitializeIteration()
    {
        // call the superclass  implementation
        Superclass::InitializeIteration();

        // set the gradient selection flag
        DemonsRegistrationFunctionType *drfp =
        dynamic_cast< DemonsRegistrationFunctionType * >
        ( this->GetDifferenceFunction().GetPointer() );

        if ( !drfp )
        {
            itkExceptionMacro(
                              << "Could not cast difference function to DemonsRegistrationFunction");
        }

        drfp->SetUseMovingImageGradient(m_UseMovingImageGradient);

        /**
         * Smooth the deformation field
         */
        if ( this->GetSmoothDisplacementField() )
        {
            this->SmoothDisplacementField();
        }
    }

    /**
     * Get the metric value from the difference function
     */
    template< class TFixedImage, class TMovingImage, class TDisplacementField >
    double
    DemonsRegistrationFilter< TFixedImage, TMovingImage, TDisplacementField >
    ::GetMetric() const
    {
        DemonsRegistrationFunctionType *drfp =
        dynamic_cast< DemonsRegistrationFunctionType * >
        ( this->GetDifferenceFunction().GetPointer() );

        if ( !drfp )
        {
            itkExceptionMacro(
                              << "Could not cast difference function to DemonsRegistrationFunction");
        }

        return drfp->GetMetric();
    }

    /**
     *
     */
    template< class TFixedImage, class TMovingImage, class TDisplacementField >
    double
    DemonsRegistrationFilter< TFixedImage, TMovingImage, TDisplacementField >
    ::GetIntensityDifferenceThreshold() const
    {
        DemonsRegistrationFunctionType *drfp =
        dynamic_cast< DemonsRegistrationFunctionType * >
        ( this->GetDifferenceFunction().GetPointer() );

        if ( !drfp )
        {
            itkExceptionMacro(
                              << "Could not cast difference function to DemonsRegistrationFunction");
        }

        return drfp->GetIntensityDifferenceThreshold();
    }

    /**
     *
     */
    template< class TFixedImage, class TMovingImage, class TDisplacementField >
    void
    DemonsRegistrationFilter< TFixedImage, TMovingImage, TDisplacementField >
    ::SetIntensityDifferenceThreshold(double threshold)
    {
        DemonsRegistrationFunctionType *drfp =
        dynamic_cast< DemonsRegistrationFunctionType * >
        ( this->GetDifferenceFunction().GetPointer() );

        if ( !drfp )
        {
            itkExceptionMacro(
                              << "Could not cast difference function to DemonsRegistrationFunction");
        }

        drfp->SetIntensityDifferenceThreshold(threshold);
    }

    /**
     * Get the metric value from the difference function
     */
    template< class TFixedImage, class TMovingImage, class TDisplacementField >
    void
    DemonsRegistrationFilter< TFixedImage, TMovingImage, TDisplacementField >
    ::ApplyUpdate(const TimeStepType& dt)
    {
        // If we smooth the update buffer before applying it, then the are
        // approximating a viscuous problem as opposed to an elastic problem
        if ( this->GetSmoothUpdateField() )
        {
            this->SmoothUpdateField();
        }
        
        this->Superclass::ApplyUpdate(dt);
        
        DemonsRegistrationFunctionType *drfp =
        dynamic_cast< DemonsRegistrationFunctionType * >
        ( this->GetDifferenceFunction().GetPointer() );
        
        if ( !drfp )
        {
            itkExceptionMacro(
                              << "Could not cast difference function to DemonsRegistrationFunction");
        }
        
        this->SetRMSChange( drfp->GetRMSChange() );
    }

}

