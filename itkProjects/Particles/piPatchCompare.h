//
//  piPatchCompare.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/13/13.
//
//

#ifndef __ParticleGuidedRegistration__piPatchCompare__
#define __ParticleGuidedRegistration__piPatchCompare__

#include <iostream>
#include "piImageDef.h"
#include "piParticleSystem.h"

#define PATCH_SIZE 5

namespace pi {

    typedef itk::Vector<ImageReal,PATCH_SIZE*PATCH_SIZE> LocalPatch;
    typedef itk::Image<LocalPatch,__Dim> PatchImage;
    typedef std::vector<PatchImage::Pointer> PatchImageVector;

    class PatchCompare {
    public:
        PatchCompare() {}
        virtual ~PatchCompare() {}

        void setTargetRadius(int r);
        void setPatchRadius(int r);

        void setAtlasImages(PatchImageVector atlasImages);
        void setAtlasLabels(LabelImageVector atlasLabels);
        void setTargetImage(RealImage::Pointer target);
        void setTargetROI(LabelImage::Pointer target);
        void setParticleSystem(ParticleSystem* system);

        // return the estimated target label
        LabelImage::Pointer getTargetLabel();

        void estimateLabel(int searchRadius, int kNearest);

        PatchImage::Pointer buildPatchImage(RealImage::Pointer image);
    private:
        ParticleSystem* _system;
        RealImage::Pointer _targetImage;
        LabelImage::Pointer _targetROI;
        LabelImage::Pointer _targetLabel;
        PatchImageVector _atlasImages;
        LabelImageVector _atlasLabels;

        int _targetRadius;
        int _patchRadius;
    };


    /**
     * Compute label transfer from string arguments
     * \param args a list of string arguments
     * \param output filename for estimated label transfer
     */
    void transferLabelsWithPatch(StringVector& args, std::string output, int searchRadius = 3, int kNearest = 3);
}

#endif /* defined(__ParticleGuidedRegistration__piPatchCompare__) */
