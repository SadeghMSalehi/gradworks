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

#include "piOptions.h"
#include "piConfigFile.h"
#include <cmath>


namespace pi {
#pragma mark DemonsRunner
    class DemonsRunner {
    public:
        DemonsRunner(Options& opts, StringVector& args);

        void buildPatches(libconfig::Setting& setting);
        void computePatchMapping();
        void computeOpticalFlowMapping();

        RealImage::Pointer deformImage(RealImage::Pointer input, DisplacementFieldType::Pointer displacement, RealImage::Pointer refImage);

        DisplacementFieldType::Pointer resampleField(DisplacementFieldType::Pointer currentField, DisplacementFieldType::Pointer resamplingField);

    private:
        ConfigFile _config;
    };

#pragma mark ParticleRunner

    class PxTools {
    public:
        bool checkFile(std::string filename) {
            ifstream i(filename.c_str());
            bool file = i.is_open();
            i.close();
            return file;
        }

        template<class T>
        void readImages(std::vector<typename T::Pointer>& data, libconfig::Setting& files) {
            for (int i = 0; i < files.getLength(); i++) {
                std::string in = files[i];
                ImageIO<T> io;
                data.push_back(io.ReadCastedImage(in));
            }
        }

        template<class T>
        void writeImages(std::vector<typename T::Pointer>& data, libconfig::Setting& files) {
            for (int i = 0; i < files.getLength(); i++) {
                std::string in = files[i];
                ImageIO<T> io;
                io.WriteImage(in, data[i]);
            }
        }
    };

    class Px {
    public:
        double x[__Dim];
        Px() {}
        Px(double v) {
            fordim (k) x[k] = v;
        }
        inline bool isnan() {
            fordim (k) {
                if (x[k] != x[k]) {
                    return true;
                }
            }
            return false;
        }
        /*
        inline Px& operator=(int v) {
            fordim (k) x[k] = v;
            return (*this);
        }
         */
        inline double& operator[](int i) {
            return x[i];
        }
        inline const double& operator[](int i) const {
            return x[i];
        }
        inline Px& operator=(const double v) {
            fordim (k) x[k] = v;
            return (*this);
        }
        inline Px& operator+=(const Px& p) {
            fordim (k) x[k] += p.x[k];
            return (*this);
        }
        inline Px& operator-=(const Px& p) {
            fordim (k) x[k] -= p.x[k];
            return (*this);
        }
        inline double operator*(const Px& p) {
            double d = 0;
            fordim (k) d += x[k] += p.x[k];
            return d;
        }
        inline double dist(const Px& o) {
            return sqrt(dist2(o));
        }
        inline double dist2(const Px& o) {
            double s = 0;
            for (int k = 0; k < __Dim; k++) {
                s += ((x[k] - o.x[k])*(x[k] - o.x[k]));
            }
            return s;
        }
        void dot(const Px& o, const Px& f) {
            fordim (k) x[k] = o[k] * f[k];
        }
        typedef double Elem;
        typedef std::vector<Px> Vector;
    };

    class PxA {
    public:
        typedef std::vector<PxA> Vector;

        int label;
        bool bound;

        PxA() { label = 0; }
    };

    // utility operator overloading
    std::ostream& operator<<(std::ostream& os, const Px& par);
    std::ostream& operator<<(std::ostream& os, const Px::Vector& par);


    class PxR {
    public:
        typedef std::vector<PxR> Vector;
        int subjId;
        LabelImage::Pointer labelmap;
        VectorImage::Pointer distmap;
        GradientImage::Pointer gradmap;
        double repulsionParams[3] = { 1, 3, 5 };

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

    /// globally affected parameters
    class PxGlobal {
    public:
        int nsubjs;
        int nlabels;
        int totalParticles;
        bool useLocalRepulsion;

        PxGlobal(): nsubjs(0), nlabels(0), totalParticles(0), useLocalRepulsion(false) {}
        std::vector<int> numParticles;
        std::vector<double> cutoffParams;
        std::vector<double> sigmaParams;
        std::vector<double> repulsionCoeff;
        std::vector<double> ensembleCoeff;

        typedef std::vector<IntVector> Neighbors;
        Neighbors neighbors;

        void load(ConfigFile& config);
    };


    class PxSubj;

    class PxAffine {
    public:
        VNLDoubleMatrix r;

        PxAffine(Px::Vector& p, Px::Vector& q): _p(p), _q(q) {
            r.set_size(__Dim, __Dim);
            r.set_identity();
        }
        void estimateAffineTransform(PxSubj* b);
        void transformVector(Px::Elem* f, Px::Elem *fOut);

    private:
        friend class PxSubj;

        PxGlobal* global;
        Px::Vector& _p;
        Px::Vector& _q;
    };


    /// subject for particle registration
    class PxSubj {
    private:
        PxGlobal* global;

    public:
        typedef std::vector<PxSubj> Vector;

        Px::Vector particles;
        Px::Vector affineAligned;
        Px::Vector forces;
        Px::Vector ensembleForces;
        PxA::Vector attrs;
        PxR::Vector regions;
        PxAffine affineTxf;

        PxSubj(): global(NULL), affineTxf(particles, affineAligned) {}

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
            attrs.resize(n);
        }
        inline void clearForce() {
            std::fill(forces.begin(), forces.end(), 0.0);
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

    class PxSystem {
    public:

        PxGlobal global;

        PxSubj sampler;
        PxSubj::Vector subjs;
        PxEnsemble ensemble;

        void main(Options& opts, StringVector& args);

    private:
        double t0, dt, t1;
        ConfigFile _config;

        void clearForces();



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
        void print();
        PxSystem _system;
    };

}

#endif /* defined(__ParticleGuidedRegistration__piParticleRunner__) */
