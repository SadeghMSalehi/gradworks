//
//  piPlutoCore.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/2/13.
//
//

#include "piPlutoCore.h"

namespace pi {
    ImageIO<RealImage3> __real3IO;
    ImageIO<RealImage2> __real2IO;




    template <>
    class ImageSamplesMSE<RealImage2, GradientImage2> {
    public:
        void computeGradient(RealSamples2& model, RealSamples2& test) {
            const int nPoints = model.n;
            float* i1 = model._values;
            float* i0 = test._values;
            float* g = test._gradients;

            for (int i = 0; i < nPoints; i++) {
                Particle& p = test._particles[0][i];
                for (int k = 0; k < 2; k++) {
                    p.f[k] = 0;
                }
                for (int j = 0; j < model.s; j++) {
                    float delta = (*i1 - *i0);
                    ++i1;
                    ++i0;

                    for (int k = 0; k < 2; k++) {
                        p.f[k] += (-delta * (*g));
                        ++g;
                    }
                }
                p.f[0] /= test.s;
                p.f[1] /= test.s;
            }
        }

        void updateParticles(RealSamples2& test, double dt) {
            for (int i = 0; i < test.n; i++) {
                Particle& p = test._particles[0][i];
                for (int k = 0; k < 2; k++) {
                    p.x[k] += dt * p.f[k];
                }
            }
            cout << test._particles[0][0] << endl;
        }
    };
    typedef ImageSamplesMSE<RealImage2, GradientImage2> RealSamplesMSE;
    

    PlutoCore::PlutoCore(QObject* parent): QObject(parent) {
        _samples = NULL;
    }

    PlutoCore::~PlutoCore() {
        if (_samples) {
            delete[] _samples;
            _samples = NULL;
        }
    }

    void PlutoCore::setImages(RealImage2Vector images) {
        if (_samples) {
            delete[] _samples;
        }
        _images = images;
        _gradientImages = __realImageTools.computeGaussianGradient(_images, 0.5);
        _samples = new RealSamples2[_images.size()];
    }
    
    void PlutoCore::setInitialParticles(ParticleVector initialParticles) {
        _initialParticles = initialParticles;
    }
    

    std::vector<ParticleVector>& PlutoCore::getParticles() {
        return _particles;
    }

    void PlutoCore::run() {
        
    }
    
    void PlutoCore::initialize() {
        const int radius = 3;
        const int size = 7;
        
        RealImage2::RegionType region;
        region.SetIndex(0, -radius);
        region.SetIndex(1, -radius);
        region.SetSize(0, size);
        region.SetSize(1, size);
        
        _particles.resize(_images.size());
        for (int i = 0; i < _images.size(); i++) {
            // set up each particle per image
            _particles[i] = _initialParticles;
            
            // interpolator creation
            typedef itk::LinearInterpolateImageFunction<RealImage2> InterpolatorType;
            InterpolatorType::Pointer interp = InterpolatorType::New();
            interp->SetInputImage(_images[i]);
            
            Gradient2InterpolatorType::Pointer gradientInterp = Gradient2InterpolatorType::New();
            gradientInterp->SetInputImage(_gradientImages[i]);
            
            // sampler initialization
            RealSamples2& samples = _samples[i];
            samples.setSampleRegion(region);
            samples.addInterpolator(interp);
            samples.addGradientInterpolator(gradientInterp);
            samples.addParticles(&_particles[i][0], _particles[i].size());
            samples.allocateValues();
            samples.allocateGradients();
            samples.sampleValues();
            samples.sampleGradients();
        }
    }


    void PlutoCore::track(int model, int test) {
        _samples[model].sampleValues();
        _samples[model].sampleGradients();


        RealSamplesMSE mse;

        // 4) Repeat until converge
        for (int i = 0; i < 10; i++) {
            // 1) Compute data
            _samples[test].sampleValues();
            _samples[test].sampleGradients();

            // 2) Compute gradient from current data

            mse.computeGradient(_samples[model], _samples[test]);

            // 3) Update particle positions
            mse.updateParticles(_samples[test], 0.05);
        }
    }

}