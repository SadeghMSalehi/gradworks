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

        PatchImage::Pointer outputImage = io.NewImageT(imageSize);
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


    void PatchCompare::performDeformableRegistration(PatchImage::Pointer fixedImage, PatchImage::Pointer movingImage, RealImage::Pointer movingSource) {
        typedef itk::TranslationTransform<double, __Dim> TransformType;
        typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
        typedef itk::MeanSquaresImageToImageMetric<PatchImage, PatchImage> MetricType;
        typedef itk:: LinearInterpolateImageFunction<PatchImage, double> InterpolatorType;
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
}
