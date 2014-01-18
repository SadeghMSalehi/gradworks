//
//  piParticleRunner.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/25/13.
//
//

#ifndef __ParticleGuidedRegistration__piParticleRunner__
#define __ParticleGuidedRegistration__piParticleRunner__

#include <iostream>
#include <cmath>
#include <itkVectorLinearInterpolateImageFunction.h>

#include "piOptions.h"
#include "piConfigFile.h"
#include "piPx.h"
#include "piTools.h"



namespace pi {
#pragma mark ParticleRunner
    const int __PatchSize = 5;

    /// globally affected parameters
    class PxGlobal {
    public:
        int nsubjs;
        int nlabels;
        int totalParticles;
        bool useLocalRepulsion;
        bool useAffineTransform;
        bool useEnsembleForce;

        PxGlobal(): nsubjs(0), nlabels(0), totalParticles(0), useLocalRepulsion(false), useAffineTransform(false), useEnsembleForce(true) {}
        std::vector<int> numParticles;
        std::vector<double> cutoffParams;
        std::vector<double> sigmaParams;
        std::vector<double> repulsionCoeff;
        std::vector<double> ensembleCoeff;
        std::vector<double> imageCoeff;

        typedef std::vector<IntVector> Neighbors;
        Neighbors neighbors;

        void load(ConfigFile& config);
    };


    class PxI {
    public:
        RealImage::Pointer image;
        GradientImage::Pointer gradient;

        LinearImageInterpolatorType::Pointer imageSampler;
        GradientInterpolatorType::Pointer gradientSampler;

        PxI();
        void load(std::string file);
    };

    class PxR {
    public:
        typedef std::vector<PxR> Vector;
        int subjId;
        LabelImage::Pointer labelmap;
        VectorImage::Pointer distmap;
        GradientImage::Pointer gradmap;

        typedef itk::NearestNeighborInterpolateImageFunction<LabelImage> LabelImageInterpolatorType;
        LabelImageInterpolatorType::Pointer lablIntp;

        LinearVectorImageInterpolatorType::Pointer distIntp;

        typedef itk::VectorLinearInterpolateImageFunction<GradientImage> GradientImageInterpolatorType;
        GradientImageInterpolatorType::Pointer gradIntp;

        bool isIn(Px& p);
        void computeNormal(Px& p, Px& nOut);

        bool projectParticle(Px& p);
        bool normalForce(Px& p, Px& f, Px& fout);

        bool load(libconfig::Setting& cache, bool error = true);
        void deriveImages(libconfig::Setting& cache, LabelImage::Pointer label, bool saveLabel);

    };



    class PxSubj;

    class PxAffine {
    public:
        VNLDoubleMatrix r;

        PxAffine(Px::Vector* p, Px::Vector* q): _p(p), _q(q) {
            r.set_size(__Dim, __Dim);
            r.set_identity();
        }
        void estimateAffineTransform(PxSubj* b);
        void transformVector(Px::Elem* f, Px::Elem *fOut);

    private:
        friend class PxSubj;

        PxGlobal* global;
        Px::Vector* _p;
        Px::Vector* _q;
    };


    /// subject for particle registration
    class PxSubj {
    private:
        PxGlobal* global;

    public:
        typedef std::vector<PxSubj> Vector;

        PxI image;
        Px::Vector particles;
        Px::Vector affineAligned;
        Px::Vector forces;
        Px::Vector ensembleForces;
        Px::Vector imageForces;
        PxA::Vector attrs;
        PxR::Vector regions;
        PxAffine affineTxf;

        PxSubj(): global(NULL), affineTxf(&particles, &affineAligned) {}

        inline void setGlobalConfig(PxGlobal* g) {
            global = g;
            affineTxf.global = g;
        }

        inline const Px& operator[](int ix) const { return particles[ix]; }
        inline Px& operator[](int ix) { return particles[ix]; }

        inline int size() const {
            return particles.size();
        }
        inline void resize(int n) {
            particles.resize(n);
            affineAligned.resize(n);
            forces.resize(n);
            ensembleForces.resize(n);
            imageForces.resize(n);
            attrs.resize(n);
        }
        inline void clearForce() {
            std::fill(forces.begin(), forces.end(), 0.0);
            std::fill(ensembleForces.begin(), ensembleForces.end(), 0.0);
            std::fill(imageForces.begin(), imageForces.end(), 0.0);

            
        }
        inline void clearVector(Px::Vector& v) {
            std::fill(v.begin(), v.end(), 0.0);
        }

        LabelImage::IndexType getIndex(int i);

        void sampleParticles(std::vector<int>& numParticles);

        void computeRepulsion(bool ignoreNeighbors = false);
        void constrainParticles();
        void constrainForces();
        void updateSystem(double dt);

        void save(ostream& os);
        bool load(std::string filename);


    };
    std::ostream& operator<<(std::ostream& os, const PxSubj& par);

    class PxEnsemble {
    public:
        PxGlobal* global;
        PxEnsemble(): global(NULL) {}

        void computeAttraction(PxSubj::Vector& subjs);
    };


    class NeighborSampler;

    class PxImageTerm {
    public:
        typedef std::vector<float> PixelVector;

        PxGlobal* global;
        PxImageTerm(): global(NULL) {}

        void computePixelEntropy(int i, int npx, int nsx, int nex, PxSubj::Vector& subjs, NeighborSampler* sampler);
        void computeImageTerm(PxSubj::Vector& subjs, int w);

    private:
        RealImage::RegionType patchRegion;

    };

    class PxSystem {
    public:

        PxGlobal global;

        // what is sampler?
        // maybe the intersection?
        PxSubj sampler;
        PxSubj::Vector subjs;
        PxEnsemble ensemble;

        void main(Options& opts, StringVector& args);

    private:
        double t0, dt, t1;
        ConfigFile _config;

//        void clearForces();



        void initialize(Options& opts, StringVector& args);
        void loadSystem(ConfigFile& config);
        bool loadSampler(ConfigFile& config);
        bool loadParticles(ConfigFile& config);
        bool saveParticles(ConfigFile& config, std::string outputName);

        // create sampler
        void computeIntersection(LabelImageVector& regions);
        void initialLoop();

        /// construct particle neighbor structure
        void setupNeighbors(const IntVector& numPx, PxSubj& subj);

        /// duplicate particles from the sampler to subjects
        void transferParticles();


        /// main particle registration loop
        void loop();

        /// affine transform particles
        void affineTransformParticles();

        /// warp result images
        void warpLabels(libconfig::Setting& data);
        void warpImages(libconfig::Setting& data);


        // auxilirary (unimportnat) member functions
        void print();
        void saveAdjacencyMatrix(std::string file);


    };

    class ParticleRunner {
    public:
        void main(Options& opts, StringVector& args);

    private:
        void printHelp();
        PxSystem _system;
    };

}

#endif /* defined(__ParticleGuidedRegistration__piParticleRunner__) */
