//
//  piContourSystem.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/6/13.
//
//

#ifndef __ParticleGuidedRegistration__piContourSystem__
#define __ParticleGuidedRegistration__piContourSystem__

#include <iostream>
#include "piParticle.h"
#include "piParticleCore.h"
#include "itkImage.h"

namespace pi {
    typedef itk::Image<DataReal,2> RealSlice;
    typedef itk::Image<LabelPixel,2> LabelSlice;

    typedef std::vector<RealSlice::Pointer> RealSliceVector;
    typedef std::vector<LabelSlice::Pointer> LabelSliceVector;

    template <class T>
    class EntropyComputer;

    class ContourSystem {
    public:
        const int NAttrs;
    private:
        ParticleSystem _system;
        ParticleVector _particles;

        RealSliceVector _sliceImages;
        RealSliceVector _attrImages;
        EntropyComputer<float>* _attrs;

        int _width;
        int _height;
        int _sliceIdx;

    public:
        ContourSystem(): NAttrs(5), _width(0), _height(0), _sliceIdx(0) {}
        void SetInitialParticles(ParticleVector& initial) { _particles = initial; }
        void SetSlices(RealSliceVector& vector) { _sliceImages = vector; }
        void SetAttributeDimensions(int w, int h) { _width = w; _height = h; }

        void AllocateAttributeBuffer();
        void ExtractAttributes();
        void ExtractParticleAttributes(RealSlice::Pointer sliceImg,
                                       Particle& par, float* attrOut, int nWidth, int nHeight);
        void Track(int nSlices);
    };
}
#endif /* defined(__ParticleGuidedRegistration__piContourSystem__) */