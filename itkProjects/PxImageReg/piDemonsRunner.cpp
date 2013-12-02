//
//  piDemonsRunner.cpp
//  PxImageReg
//
//  Created by Joohwi Lee on 11/20/13.
//
//

#include "piDemonsRunner.h"
#include "piOptions.h"
#include "piImageIO.h"
#include "piImageProc.h"

#include <numeric>
#include <algorithm>
#include <itkResampleImageFilter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkVectorResampleImageFilter.h>
#include <itkWarpImageFilter.h>
#include <itkWarpVectorImageFilter.h>

using namespace std;
using namespace libconfig;


namespace pi {
    void executeDemonsRunner(pi::Options &opts, StringVector &args) {
        DemonsRunner runner(opts, args);
        runner.computeOpticalFlowMapping();

        /*
        if (args.size() > 0) {
            if (args[0] == "optical-flow-mapping") {
            } else if (args[0] == "patch-mapping") {
                runner.computePatchMapping();
            }
        }
         */
    }



    DemonsRunner::DemonsRunner(pi::Options &opts, StringVector &args) {
        string file = opts.GetString("--config");
        if (file != "") {
            _config.readFile(file.c_str());
        }
        cout << "config file: " << file << endl;
    }

    /*
    void DemonsRunner::computePatchMapping() {
        if (!_config.exists("dense-patch-mapping")) {
            return;
        }

        Setting& config = _config["dense-patch-mapping"];
        int fixedImageIdx = config["fixed-image-idx"];
        int movingImageIdx = config["moving-image-idx"];
        PatchImage::Pointer fixedImage = __patchImages[fixedImageIdx];
        PatchImage::Pointer movingImage = __patchImages[movingImageIdx];

        PatchImage::RegionType sourceRegion = fixedImage->GetBufferedRegion();
        PatchImage::RegionType activeRegion = _config.offsetRegion(sourceRegion, "dense-patch-mapping.region-offset");

        PatchCompare patchMaker;
        DisplacementFieldType::Pointer deformationField = patchMaker.performDenseMapping(fixedImage, movingImage, activeRegion);

        // deformation field output
        ImageIO<DisplacementFieldType> io;
        string deformationFieldFile = config["displacement-field-output"];
        io.WriteImage(deformationFieldFile, deformationField);

        // warped image output
        string warpedImage = config["warped-image-output"];
        io.WriteImageS<RealImage>(warpedImage, deformImage(_config.image(movingImageIdx), deformationField, _config.image(fixedImageIdx)));
    }
     */

    DisplacementFieldType::Pointer DemonsRunner::computeDemonsFlow(RealImage::Pointer fixed, RealImage::Pointer moving, double dt) {
        ImageIO<DisplacementFieldType> io;

        // allocate output
        DisplacementFieldType::Pointer flowOutput = io.NewImageS<RealImage>(fixed);

        // compute gradient image
        RealImage::SpacingType spacing = fixed->GetSpacing();
        GradientImage::Pointer gF = ComputeGaussianGradient(fixed, spacing[0]/2.0);
        GradientImage::Pointer gM = ComputeGaussianGradient(moving, spacing[0]/2.0);

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
            } else {
                for (int j = 0; j < DisplacementFieldType::ImageDimension; j++) {
                    (*oBuf)[j] = 0;
                }
            }
            fBuf++;
            mBuf++;
            gFv++;
            gMv++;
            oBuf++;
        }

        return flowOutput;
    }


    /// Compute an intermediate image to record the transformed coordinate
    ///
    /// \param sourceImage an image defines the space of the fixed image
    ///
    DisplacementFieldType::Pointer DemonsRunner::computeCoordImage(RealImage::Pointer sourceImage) {
        ImageIO<DisplacementFieldType> io;
        DisplacementFieldType::Pointer coordImage = io.NewImageS<RealImage>(sourceImage);

        DisplacementFieldType::PixelType* oBuf = coordImage->GetBufferPointer();
        const int nPixels = coordImage->GetPixelContainer()->Size();

        itk::ImageRegionConstIteratorWithIndex<DisplacementFieldType> iter(coordImage, coordImage->GetBufferedRegion());
        for (int i = 0; i < nPixels; i++) {
            DisplacementFieldType::IndexType idx = iter.GetIndex();
            for (int j = 0; j < DisplacementFieldType::ImageDimension; j++) {
                oBuf[i][j] = idx[j];
            }
            ++iter;
        }

        return coordImage;
    }



    void DemonsRunner::computeOpticalFlowMapping() {
        Setting& files = _config.lookup("demons.files");
        Setting& param = _config.lookup("demons.params");

        double dt = 0.1;
        param.lookupValue("timestep", dt);

        int iter = 1;
        param.lookupValue("iteration", iter);

        double fieldSigma = -1;
        param.lookupValue("field-sigma", fieldSigma);

        double velocitySigma = -1;
        param.lookupValue("velocity-sigma", velocitySigma);

        cout << "dt = " << dt << endl;
        cout << "iterations  = " << iter << endl;
        cout << "field-sigma = " << fieldSigma << endl;
        cout << "velocity-sigma = " << velocitySigma << endl;

        ImageIO<DisplacementFieldType> fieldIO;
        ImageIO<RealImage> io;

        for (int i = 0; i < files.getLength(); i++) {
            string source = files[i][0];
            string target = files[i][1];
            string flowoutput = files[i][2];
            string warpoutput = files[i][3];

            cout << source << " => " << target << "; " << flowoutput << ", " << warpoutput << endl;

            RealImage::Pointer sourceImage = io.ReadCastedImage(source);
            RealImage::Pointer targetImage = io.ReadCastedImage(target);

            DisplacementFieldType::Pointer coordImage = computeCoordImage(sourceImage);

            DisplacementFieldType::Pointer velocityImage;


            // warped image at each iteration
            RealImage::Pointer tempImage = io.CopyImage(sourceImage);

            for (int j = 0; j < iter; j++) {
                velocityImage = computeDemonsFlow(targetImage, tempImage, dt);

                if (velocitySigma > 0) {
                    velocityImage = applyGaussian<DisplacementFieldType>(velocityImage, velocitySigma);
                }

                // compute approximated update field
                coordImage = deformVectorImage(coordImage, velocityImage, targetImage);
                if (fieldSigma > 0) {
                    coordImage = applyGaussian<DisplacementFieldType>(coordImage, fieldSigma);
                }

//                fieldIO.WriteImage("tempflow.nii.gz", coordImage);
                tempImage = resampleImage(sourceImage, coordImage);
            }
            io.WriteImage(warpoutput, tempImage);
            fieldIO.WriteImage(flowoutput, coordImage);
        }
    }


    /// create a warped image from sourceImage in according to the coordImage
    /// The pixel value of coordImage is the index before transformation
    /// Each voxel of outputImage will sample a pixel value from sourceImage at the coordImage designates
    ///
    RealImage::Pointer DemonsRunner::resampleImage(RealImage::Pointer sourceImage, DisplacementFieldType::Pointer coordImage) {
        ImageIO<RealImage> io;
        RealImage::Pointer outputImage = io.NewImage(sourceImage);

        // assume that the size of sourceImage and coordImage is same
        DisplacementFieldType::PixelType* dBuf = coordImage->GetBufferPointer();
        RealImage::PixelType* oBuf = outputImage->GetBufferPointer();
        int nPixels = outputImage->GetPixelContainer()->Size();


        // use linear interpolation to sample
        typedef itk::LinearInterpolateImageFunction<RealImage> ImageInterpolator;
        ImageInterpolator::Pointer sampler = ImageInterpolator::New();
        sampler->SetInputImage(sourceImage);

        // for each pixel of outputImage
        for (int i = 0; i < nPixels; i++) {
            // determine where to sample
            itk::ContinuousIndex<double, RealImage::ImageDimension> idx;
            for (int j = 0; j < RealImage::ImageDimension; j++) {
                idx[j] = dBuf[i][j];
            }
            // read a pixel value from the sourceImage at idx
            oBuf[i] = sampler->EvaluateAtContinuousIndex(idx);
        }
        return outputImage;
    }


    /// As described in the header file, the output pixels are resampled from the input image in according to the given displacement image.
    ///
    DisplacementFieldType::Pointer DemonsRunner::deformVectorImage(DisplacementFieldType::Pointer input, DisplacementFieldType::Pointer displacement, RealImage::Pointer refImage) {
        typedef itk::WarpVectorImageFilter<DisplacementFieldType, DisplacementFieldType, DisplacementFieldType> WarpImageFilter;

        // create an instance of warp image filter
        WarpImageFilter::Pointer warpFilter = WarpImageFilter::New();
        warpFilter->SetInput(input);
        warpFilter->SetDisplacementField(displacement);

        // the output space should be defined in according to the refImage
        warpFilter->SetOutputDirection(refImage->GetDirection());
        warpFilter->SetOutputOrigin(refImage->GetOrigin());
        warpFilter->SetOutputSpacing(refImage->GetSpacing());

        // The LargetPossibleRegion for the output is inherited from the input displacement field.
        // There is no need to assign output size because it is determined by the displacement field.
        // This is inconsistent with WarpImageFilter
        // warpFilter->SetOutputSize(refImage->GetBufferedRegion().GetSize());
        warpFilter->Update();

        // return the produced resampled image
        DisplacementFieldType::Pointer warpedInput = warpFilter->GetOutput();
        warpedInput->DisconnectPipeline();
        return warpedInput;
    }


    RealImage::Pointer DemonsRunner::deformImage(RealImage::Pointer input, DisplacementFieldType::Pointer displacement, RealImage::Pointer refImage) {
        FieldTransformType::Pointer transform = FieldTransformType::New();
        transform->SetDisplacementField(displacement);

        typedef itk::ResampleImageFilter<RealImage, RealImage> ResampleFilterType;

        ResampleFilterType::Pointer resampler = ResampleFilterType::New();
        resampler->SetInput(input);
        resampler->SetUseReferenceImage(true);
        resampler->SetReferenceImage(refImage);
        resampler->SetTransform(transform);
        resampler->Update();
        RealImage::Pointer warpedImage = resampler->GetOutput();
        warpedImage->DisconnectPipeline();
        return warpedImage;
    }



    // resample displaecement field
    DisplacementFieldType::Pointer DemonsRunner::resampleField(DisplacementFieldType::Pointer currentField, DisplacementFieldType::Pointer resamplingField) {
        typedef itk::VectorLinearInterpolateImageFunction<DisplacementFieldType> InterpolatorType;
        InterpolatorType::Pointer interpolator = InterpolatorType::New();

        FieldTransformType::Pointer transform = FieldTransformType::New();
        transform->SetDisplacementField(currentField);
        transform->SetInterpolator(interpolator);

        typedef itk::VectorResampleImageFilter<DisplacementFieldType, DisplacementFieldType> ResampleFilter;
        ResampleFilter::Pointer resampler = ResampleFilter::New();

        resampler->SetInput(resamplingField);
        resampler->SetSize(currentField->GetBufferedRegion().GetSize());
        resampler->SetOutputOrigin(currentField->GetOrigin());
        resampler->SetOutputSpacing(currentField->GetSpacing());
        resampler->SetOutputDirection(currentField->GetDirection());
        resampler->SetTransform(transform);

        resampler->Update();

        DisplacementFieldType::Pointer resampledField = resampler->GetOutput();
        resampledField->DisconnectPipeline();
        return resampledField;
    }




//    void DemonsRunner::buildPatches(libconfig::Setting &setting) {
//        bool checkFiles = setting["check-files"];
//        PatchCompare patchMaker;
//        string inputImageTag = setting["input"];
//        string patchImageTag = "patch-images";
//        for (int i = 0; i < _config[patchImageTag].getLength(); i++) {
//            string output = _config[patchImageTag][i];
//            if (!checkFiles || !io.FileExists(output.c_str())) {
//                string input = _config[inputImageTag][i];
//                RealImage::Pointer inputImage = io.ReadCastedImage(input);
//                PatchImage::Pointer patchImage = patchMaker.buildPatchImage(inputImage);
//                io.WriteImageS<PatchImage>(output, patchImage);
//            }
//        }
//    }




}
