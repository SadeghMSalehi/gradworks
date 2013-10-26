//
//  piParticleRunner.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/25/13.
//
//

#include "piImageDef.h"
#include "piParticleRunner.h"
#include "piPatchCompare.h"
#include "piImageIO.h"
#include "piImageProcessing.h"

#include <itkResampleImageFilter.h>

using namespace libconfig;
using namespace std;

namespace pi {
    typedef itk::ResampleImageFilter<RealImage, RealImage> ResampleImageFilterType;

    static ImageIO<RealImage> io;
    vector<PatchImage::Pointer> __patchImages;

    void executeParticleRunner(pi::Options &opts, StringVector &args) {
        ParticleRunner runner;
        runner.main(opts, args);
    }

    void ParticleRunner::buildPatches(libconfig::Setting &setting) {
        bool checkFiles = setting["check-files"];
        PatchCompare patchMaker;
        string inputImageTag = setting["input"];
        string patchImageTag = "patch-images";
        for (int i = 0; i < _config[patchImageTag].getLength(); i++) {
            string output = _config[patchImageTag][i];
            if (!checkFiles || !io.FileExists(output.c_str())) {
                string input = _config[inputImageTag][i];
                RealImage::Pointer inputImage = io.ReadCastedImage(input);
                PatchImage::Pointer patchImage = patchMaker.buildPatchImage(inputImage);
                io.WriteImageS<PatchImage>(output, patchImage);
            }
        }
    }

    void ParticleRunner::initialize(Options& opts, StringVector& args) {
        _config.load(opts.GetConfigFile());

        if (_config.exists("build-patch")) {
            buildPatches(_config["build-patch"]);
        }

        _config.readImages<PatchImage>(__patchImages, "patch-images");
    }

    void ParticleRunner::computePatchMapping() {
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
        FieldTransformType::Pointer transform = FieldTransformType::New();
        transform->SetDisplacementField(deformationField);

        ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
        resampler->SetInput(_config.image(movingImageIdx));
        resampler->SetUseReferenceImage(true);
        resampler->SetReferenceImage(_config.image(fixedImageIdx));
        resampler->SetTransform(transform);
        resampler->Update();
        io.WriteImageS<RealImage>(warpedImage, resampler->GetOutput());
    }

    void ParticleRunner::main(pi::Options &opts, StringVector &args) {
        initialize(opts, args);
        computePatchMapping();
    }

    void ParticleRunner::print() {

    }
}