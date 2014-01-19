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

    /// A class stores configuration used as global
    /// @property nsubjs the number of subjects
    /// @property nlabels the number of labels. This should be the same across subjects
    /// @property totalParticles the total number of particles of a subject
    /// @property useLocalRepulsion
    /// @property useAffineTransform
    /// @property useEnsembleForce
    /// @property numParticles a vector stores the number of particles for each label
    /// @property cutoffParams the cut-off parameters used in the sampling process
    /// @property sigmaParams the sigma parameters used in the sampling process for controlling spacing
    /// @property repulsionCoeff the distribution weight between terms
    /// @property ensembleCoeff the ensemble weight between terms
    /// @property imageCoeff the image term weight
    /// @property Neighbors ???
    /// @method load to load a config file
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

    /// A class represents an image to be registered.
    /// @property image an intensity image
    /// @property gradient the gradient of the image
    /// @method load(std::string file)

    class PxI {
    public:
        RealImage::Pointer image;
        GradientImage::Pointer gradient;

        LinearImageInterpolatorType::Pointer imageSampler;
        GradientInterpolatorType::Pointer gradientSampler;

        PxI();
        void load(std::string file);
    };


    /// A class to represent a label region
    /// A subset of particles is allowed move freely only in a given region
    /// This class abstracts such a region and provides information to constrain particles.
    /// For the constraint, it needs a label defined in the image domain,
    /// a distance map to relocate a particle into the closest boundary,
    /// and the gradient map to identify the boundary of the label.
    /// @property subjId
    /// @property label  which label is this
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



    // forward declaration for class
    class PxSubj;

    /// A class to estimate affine transformation between subjects
    /// @property r a matrix store the transformation
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


    /// This class represents a subject to store image and particle information
    /// @brief A class encapsulates a subject information
    class PxSubj {
    private:
        PxGlobal* global;

    public:
        typedef std::vector<PxSubj> Vector;

        PxI image;

        Px::Vector particles;           /// a vector of particles in the subject space
        Px::Vector affineAligned;       /// a vector of particles after the affine alignment
        Px::Vector forces;              /// a vector of net forces applied to each particle
        Px::Vector ensembleForces;      /// a vector of forces from the ensemble term
        Px::Vector imageForces;         /// a vector of forces from the image term
        PxA::Vector attrs;              /// a vector of attributes for particles
        PxR::Vector regions;            /// a vector of regions for labels
        PxAffine affineTxf;             /// affine transform for this subject

        PxSubj(): global(NULL), affineTxf(&particles, &affineAligned) {}

        /// setting the global configuration
        /// @param g a pointer to PxGlobal
        inline void setGlobalConfig(PxGlobal* g) {
            global = g;
            affineTxf.global = g;
        }

        /// access a particle in the subject space (const version)
        inline const Px& operator[](int ix) const { return particles[ix]; }

        inline Px& operator[](int ix) { return particles[ix]; }

        /// return the number of particles
        inline int size() const {
            return particles.size();
        }

        /// change the number of particles
        inline void resize(int n) {
            particles.resize(n);
            affineAligned.resize(n);
            forces.resize(n);
            ensembleForces.resize(n);
            imageForces.resize(n);
            attrs.resize(n);
        }

        /// clear forces used
        inline void clearForce() {
            std::fill(forces.begin(), forces.end(), 0.0);
            std::fill(ensembleForces.begin(), ensembleForces.end(), 0.0);
            std::fill(imageForces.begin(), imageForces.end(), 0.0);

            
        }

        /// a utility function to clear a particle vector
        inline void clearVector(Px::Vector& v) {
            std::fill(v.begin(), v.end(), 0.0);
        }

        LabelImage::IndexType getIndex(int i);

        /// extract particles with a list of numbers
        void sampleParticles(std::vector<int>& numParticles);

        /// compute the distribution of particles
        void computeRepulsion(bool ignoreNeighbors = false);

        /// constrain particles into a label
        void constrainParticles();

        /// constrain forces not to drive a particle toward outside
        void constrainForces();

        /// update the system with the amount of time
        void updateSystem(double dt);

        /// save this system into a file
        void save(ostream& os);

        /// load a system from a file
        bool load(std::string filename);


    };
    std::ostream& operator<<(std::ostream& os, const PxSubj& par);


    /// A class to compute the ensemble term
    class PxEnsemble {
    public:
        PxGlobal* global;
        PxEnsemble(): global(NULL) {}

        void computeAttraction(PxSubj::Vector& subjs);
    };


    // forward declaration of a class to sample a local intensity patch
    class NeighborSampler;


    /// A class to compute the image term in the equation. This class iterates over all particles and compute the gradient to minimize the intensity entropy term. Since this class requires the image intensity and the gradient of intensities, the information of images should be provided.
    class PxImageTerm {
    public:
        /// a data type to store pixel values
        typedef std::vector<float> PixelVector;

        /// a point to the global configuration
        PxGlobal* global;

        /// constructor
        PxImageTerm(): global(NULL) {}

        /// @brief Compute an entropy value for each particle. this function is repeated for every pixel in an image.
        /// @param i the index of a particle in a subject
        /// @param npx the number of particles
        /// @param nsx the number of subjects
        /// @param nex the number of elements in a local intensity patch
        /// @param subjs the vector of subjects
        /// @param sampler the intensity sampler which has the size of a patch
        void computePixelEntropy(int i, int npx, int nsx, int nex, PxSubj::Vector& subjs, NeighborSampler* sampler);

        /// @brief Compute the image term for all subjects
        /// @param subjs the vector of subjects
        /// @param w the width of a local intensity patch
        void computeImageTerm(PxSubj::Vector& subjs, int w);

    private:
        RealImage::RegionType patchRegion;

    };


    /// A class represents the particle system
    class PxSystem {
    public:

        PxGlobal global; /// global configuration

        PxSubj sampler;         /// used in the initialization step
        PxSubj::Vector subjs;   /// the vector of subjects
        PxEnsemble ensemble;    /// a class to compute the ensemble term


        /// the main entry point
        /// @param opts Options
        /// @param args Arguments
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
