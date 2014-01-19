//
//  piParticleRunner.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/25/13.
//
//

#include "piImageDef.h"
#include "piParticleRunner.h"
#include "piImageIO.h"
#include "piParticleWarp.h"
#include "piEntropyComputer.h"
#include "piImageProc.h"

#include <numeric>
#include <algorithm>
#include <itkResampleImageFilter.h>
#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkVectorResampleImageFilter.h>

using namespace libconfig;
using namespace std;

namespace pi {
    typedef itk::ResampleImageFilter<RealImage, RealImage> ResampleImageFilterType;

    static PxTools sys;
    static ImageIO<RealImage> io;
    static ImageIO<LabelImage> labelIO;

    void executeParticleRunner(pi::Options &opts, StringVector &args) {
        ParticleRunner runner;
        runner.main(opts, args);
    }


    ostream& operator<<(ostream& os, const PxSubj& p) {
        for (int i = 0; i < p.size(); i++) {
            cout << p.particles[i] << endl;
        }
        return os;
    }


    PxI::PxI() {
        imageSampler = LinearImageInterpolatorType::New();
        gradientSampler = GradientInterpolatorType::New();
    }

    void PxI::load(std::string file) {
        this->image = io.ReadCastedImage(file);
        if (this->image.IsNull()) {
            cout << "can't read " << file << endl;
            exit(0);
        }
        this->gradient = ComputeGaussianGradient(image, image->GetSpacing()[0] / 2.0);

        this->imageSampler->SetInputImage(this->image);
        this->gradientSampler->SetInputImage(this->gradient);
    }

#pragma mark PxR implementations
    // load label map, distance map, gradient map for each region
    bool PxR::load(libconfig::Setting& cache, bool error) {
        string labelfile = cache[0];
        string distfile = cache[1];
        string gradfile = cache[2];


        if (sys.checkFile(labelfile)) {
            labelmap = io.ReadImageS<LabelImage>(labelfile);
        } else {
            if (error) {
                cout << "can't read " << labelfile << endl;
            }
            return false;
        }
        if (sys.checkFile(distfile)) {
            distmap = io.ReadImageS<VectorImage>(distfile);
        }
        if (sys.checkFile(gradfile)) {
            gradmap = io.ReadImageS<GradientImage>(gradfile);
        }

        deriveImages(cache, labelmap, false);
        return true;
    }

    void PxR::deriveImages(libconfig::Setting& cache, LabelImage::Pointer label, bool saveLabel) {
        string labelfile = cache[0];
        string distfile = cache[1];
        string gradfile = cache[2];

        if (saveLabel) {
            labelmap = label;
            io.WriteImageS<LabelImage>(labelfile, label);
        }
        if (distmap.IsNull()) {
            distmap = ComputeDistanceMap(label);
            io.WriteImageS<VectorImage>(distfile, distmap);
        }
        if (gradmap.IsNull()) {
            LabelImage::SpacingType spacing = label->GetSpacing();
            gradmap = ComputeGaussianGradient(label, spacing[0] / 2.0);
            io.WriteImageS<GradientImage>(gradfile, gradmap);
        }

        lablIntp = LabelImageInterpolatorType::New();
        lablIntp->SetInputImage(labelmap);

        distIntp = LinearVectorImageInterpolatorType::New();
        distIntp->SetInputImage(distmap);

        gradIntp = GradientImageInterpolatorType::New();
        gradIntp->SetInputImage(gradmap);

    }

    void PxR::computeNormal(Px& p, Px& nOut) {
        GradientImageInterpolatorType::InputType pi;
        fordim (k) {
            pi[k] = p[k];
        }
        GradientImageInterpolatorType::OutputType nx = gradIntp->Evaluate(pi);
        fordim (k) {
            nOut[k] = nx[k];
        }
    }

    bool PxR::isIn(Px& p) {
        LabelImageInterpolatorType::InputType pi;
        fordim (k) {
            pi[k] = p[k];
        }
        if (lablIntp->IsInsideBuffer(pi) && lablIntp->Evaluate(pi) > 0) {
            return true;
        }
        return false;
    }

    bool PxR::projectParticle(Px& p) {
        if (p.isnan()) {
            return false;
        }
        LabelImage::PointType px;
        fordim (k) {
            px[k] = p[k];
        }
        LabelImage::IndexType ix;
        labelmap->TransformPhysicalPointToIndex(px, ix);
        if (!lablIntp->IsInsideBuffer(ix)) {
            return false;
        }
        if (isIn(p)) {
            return false;
        }
        LinearVectorImageInterpolatorType::InputType pi;
        LinearVectorImageInterpolatorType::OutputType dx = distIntp->EvaluateAtIndex(ix);
        fordim (k) {
            ix[k] += dx[k];
        }
        labelmap->TransformIndexToPhysicalPoint(ix, px);
        fordim (k) {
            p[k] = px[k];
        }
        return true;
    }

    bool PxR::normalForce(Px& p, Px& f, Px& fout) {
        if (isIn(p)) {
            return false;
        }

        Px n = 0;
        computeNormal(p, n);

        fordim (k) {
            fout[k] = n[k] * abs(f[k]);
        }
        return true;
    }


#pragma mark PxSubj implementations
    LabelImage::IndexType PxSubj::getIndex(int i) {
        LabelImage::Pointer regionLabel = regions[attrs[i].label].labelmap;
        LabelImage::PointType px;
        fordim (k) {
            px[k] = particles[i][k];
        }
        LabelImage::IndexType ix;
        regionLabel->TransformPhysicalPointToIndex(px, ix);
        return ix;
    }



    /// inline function
    inline void computeRepulsionWeight(VNLMatrix& weights, int i, int j, Px& pi, Px& pj, double cutoff, double sigma2, double kappa) {
        // compute weight
        double dij = pi.dist(pj);
        if (dij < cutoff) {
            if (dij > 0) {
                weights[i][j] = exp(-dij*dij*kappa/(sigma2)) / dij;
            } else if (dij == 0) {
                cout << "particle overlap!" << endl;
            }
        } else {
            weights[i][j] = 0;
        }
    }

    // compute repulsion force
    void PxSubj::computeRepulsion(bool ignoreNeighbors) {
        const int npx = size();
        VNLMatrix weights(npx, npx);
        weights.fill(0);

        // loop over all particle pairs
        for (int i = 0; i < npx; i++) {
            Px& pi = particles[i];
            weights[i][i] = 0;
            double kappa = 1;

            const int label = attrs[i].label;
            double sigma = global->sigmaParams[label];
            double cutoff = global->cutoffParams[label];
            const double sigma2 = sigma * sigma;

            // loop over neighbors
            bool useNeighbors = global->neighbors.size() == npx && global->useLocalRepulsion;
            if (!useNeighbors || ignoreNeighbors) {
                // neighbor structure isn't yet built (sampler case)
                for (int jj = i + 1; jj < npx; jj++) {
                    Px& pj = particles[jj];
                    computeRepulsionWeight(weights, i, jj, pi, pj, cutoff, sigma2, kappa);
                    weights[jj][i] = weights[i][jj];
                }
            } else {
                for (int j = 0; j < global->neighbors[i].size(); j++) {
                    int jj = global->neighbors[i][j];
                    Px& pj = particles[jj];
                    computeRepulsionWeight(weights, i, jj, pi, pj, cutoff, sigma2, kappa);
                }
            }
        }

        // compute force
        VNLVector wsum(npx);
        for (int i = 0; i < npx; i++) {
            wsum[i] = 0;
            for (int j = 0; j < npx; j++) {
                wsum[i] += weights[i][j];
            }
        }
        for (int i = 0; i < npx; i++) {
            Px& fi = forces[i];
            Px& pi = particles[i];
            const int label = attrs[i].label;
            double coeff = global->repulsionCoeff[label];
            for (int j = i + 1; j < npx; j++) {
                Px& pj = particles[j];
                Px& fj = forces[j];
                double c = coeff * weights[i][j];
                if (wsum[i] > 1e-5 && wsum[j] > 1e-5) {
                    fordim (k) {
                        double fk = c * (pi.x[k] - pj.x[k]);
                        fi.x[k] += fk / wsum[i];
                        fj.x[k] -= fk / wsum[j];
                    }
                }
            }
        }
    }


    // update the particle dynamic equation
    void PxSubj::updateSystem(double dt) {
        Px::Vector::iterator p = particles.begin();
        Px::Vector::iterator f = forces.begin();
        Px::Vector::iterator e = ensembleForces.begin();
        Px::Vector::iterator i = imageForces.begin();

        for (; p != particles.end(); p++, f++, e++, i++) {
            fordim (k) {
                p->x[k] += (dt * (f->x[k] + e->x[k] + i->x[k]));
                assert(!std::isnan(p->x[k]));
            }
        }
    }

    // iterate all particles and move inside of the object
    void PxSubj::constrainParticles() {
        Px::Vector::iterator p = particles.begin();
        PxA::Vector::iterator a = attrs.begin();
        for (; p != particles.end(); a++, p++) {
            a->bound = regions[a->label].projectParticle(*p);
        }
    }

    // iterate all particles and remove boundary normal component of forces
    void PxSubj::constrainForces() {
        Px::Vector::iterator p = particles.begin();
        Px::Vector::iterator f = forces.begin();
        PxA::Vector::iterator a = attrs.begin();
        Px fn = 0;
        for (; f != forces.end(); p++, a++, f++) {
            if (a->bound) {
                if (regions[a->label].normalForce(*p, *f, fn)) {
                    assert(!fn.isnan());
                    assert(!f->isnan());
                    *f += fn;
                }
            }
        }
    }

    // random sample particles from each region
    // this function is only for the sampler
    void PxSubj::sampleParticles(std::vector<int>& numParticles) {
        // allocate particles
        int totalParticles = std::accumulate(numParticles.begin(), numParticles.end(), 0);
        resize(totalParticles);

        // assign labels
        int total = 0, label = 0;
        for (int i = 0; i < numParticles.size(); i++) {
            for (int j = 0; j < numParticles[i]; j++) {
                attrs[total].label = label;
                total ++;
            }
            label ++;
        }

        // loop over labels
        total = 0;
        for (unsigned int i = 0; i < regions.size(); i++) {
            // identify avaiable pixels
            std::vector<LabelImage::IndexType> pixelIds;
            LabelImageIteratorType iter(regions[i].labelmap, regions[i].labelmap->GetBufferedRegion());

            for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
                if (iter.Get() > 0) {
                    pixelIds.push_back(iter.GetIndex());
                }
            }

            // shuffle avaiable pixels and choose some
            const int npx = numParticles[i];
            std::random_shuffle(pixelIds.begin(), pixelIds.end());

            for (int j = 0; j < npx; j++) {
                LabelImage::PointType pts;
                regions[i].labelmap->TransformIndexToPhysicalPoint(pixelIds[j], pts);

                attrs[total].label = i;
                particles[total].x[0] = pts[0];
                particles[total].x[1] = pts[1];

                total ++;
            }
        }
    }


    // load particle coordinates and labels from libconfig format
    bool PxSubj::load(string filename) {
        if (!sys.checkFile(filename)) {
            return false;
        }

        libconfig::Config config;
        config.readFile(filename.c_str());

        Setting& labels = config.lookup("labels");
        Setting& data = config.lookup("particles");

        if (global->totalParticles != labels.getLength() || labels.getLength() != data[0].getLength()) {
            return false;
        }

        int npx = labels.getLength();
        resize(npx);
        for (int i = 0; i < npx; i++) {
            attrs[i].label = labels[i];
            fordim (k) {
                particles[i][k] = data[k][i];
            }
        }
        return true;
    }

    // save particle coordinates and labels in libconfig format
    void PxSubj::save(ostream &os) {
        os << "labels = [" << attrs[0].label;
        for (int i = 1; i < attrs.size(); i++) {
            os << "," << attrs[i].label;
        }
        os << "]" << endl;

        os.setf(ios::fixed, ios::floatfield);
        os << "particles = (";

        fordim (k) {
            if (k > 0) {
                os << ",";
            }
            os << "[" << particles[0][k];
            for (int i = 1; i < attrs.size(); i++) {
                os << "," << particles[i][k];
            }
            os << "]" << endl;
        }
        os << ")" << endl;
    }

#pragma mark PxEnsemble Implementations
    void PxAffine::estimateAffineTransform(pi::PxSubj *b) {
        VNLDoubleMatrix mat(_p->size(), __Dim);
        VNLDoubleMatrix mbt(_p->size(), __Dim);
        for (int i = 0; i < mat.rows(); i++) {
            fordim (k) {
                mbt[i][k] = (*_p)[i][k] - 51;
                mat[i][k] = b->particles[i][k] - 51;
            }
        }
        double ma = 0, mb = 0;
        for (int i = 0; i < mat.rows(); i++) {
            ma += mat[i][0]*mat[i][0] + mat[i][1]*mat[i][1];
            mb += mbt[i][0]*mbt[i][0] + mbt[i][1]*mbt[i][1];
        }
        ma /= mat.rows();
        mb /= mbt.rows();
        for (int i = 0; i < mat.rows(); i++) {
            fordim (k) {
                mat[i][k] /= ma;
                mbt[i][k] /= mb;
            }
        }
        VNLDoubleMatrix m = mbt.transpose() * mat;

        vnl_svd<double> svd(m);
        VNLDoubleMatrix v = svd.V();
        VNLDoubleMatrix u = svd.U();
        this->r = v.transpose() * u;

        double rm = sqrt(r[0][0]*r[0][0] + r[0][1]*r[0][1]);
        r /= rm;

        for (int i = 0; i < mat.rows(); i++) {
            (*_q)[i][0] = ma * (r[0][0] * mat[i][0] + r[0][1] * mat[i][1]) + 51;
            (*_q)[i][1] = ma * (r[1][0] * mat[i][0] + r[1][1] * mat[i][1]) + 51;
        }
    }

    void PxAffine::transformVector(Px::Elem *f, Px::Elem *fOut) {
        fordim (i) {
            fOut[i] = 0;
            fordim (j) {
                fOut[i] += r[i][j] * f[j];
            }
        }
    }

    void PxEnsemble::computeAttraction(PxSubj::Vector &subjs) {
        // compute entropy for each particle
        EntropyComputer<double> comp(global->nsubjs, global->totalParticles, __Dim);
        // copy data into the entropy computer
        comp.dataIter.FirstData();
        for (int i = 0; i < global->nsubjs; i++) {
            for (int j = 0; j < global->totalParticles; j++) {
                fordim (k) {
                    comp.dataIter.sample[k] = subjs[i].affineAligned[j][k];
                }
                comp.dataIter.NextSample();
            }
            comp.dataIter.NextData();
        }
        // move to center (subtract mean from each sample)
        comp.MoveToCenter();
        // compute covariance matrix
        if (!comp.ComputeCovariance()) {
            cout << "Wrong covariance!!" << endl;
            for (int i = 0; i < subjs.size(); i++) {
                cout << "SUBJ[" << i << "]" << endl;
                cout << subjs[i].particles << endl;
            }
        }
        // compute inverse covariance
        comp.ComputeGradient();

        EntropyComputer<double>::Iterator gradIter(comp.gradient.data_block(), global->totalParticles, __Dim);
        gradIter.FirstData();
        for (int i = 1; i < global->nsubjs; i++) {
            for (int j = 0; j < global->totalParticles; j++) {
                double coeff = global->ensembleCoeff[subjs[i].attrs[j].label];
                Px ensembleForce;
                subjs[i].affineTxf.transformVector(gradIter.sample, ensembleForce.x);
                fordim (k) {
                    subjs[i].ensembleForces[j][k] += coeff*ensembleForce[k];
                }
                gradIter.NextSample();
            }
            gradIter.NextData();
        }
    }

#pragma mark Neighbor Intensity Sampler implementations
    class NeighborSampler {
    public:
        typedef RealImage::PixelType PixelType;
        typedef RealImage::IndexType IndexType;
        typedef RealImage::PointType PointType;
        typedef RealImage::RegionType RegionType;
        typedef RealImage::SizeType SizeType;

        NeighborSampler(RegionType region, RealImage::Pointer image);
        ~NeighborSampler();

        void setSampleRegion(RegionType& region, RealImage* reference);

        /// sample image intensities from a given interpolator
        /// @interp ImageInterpolator
        /// @px particle locations at the original space
        /// @pz particle locations at the transformed space (actual sampling)
        void sampleValues(LinearImageInterpolatorType* interp, Px& px, Px& pz, double* pxValue);

        /// sample image gradients from a given interpolator
        /// @interp GradientInterpolator
        /// @px particle locations at the original space
        /// @pz particle locations at the transformed space (actual sampling)
        void sampleGradients(GradientInterpolatorType* interp, Px& px, Px& pz, Px::Vector& pxGrad);

        void createSampleIndexes2(IndexType& startIdx);
        void createSampleIndexes3(IndexType& startIdx);

    private:
        RealImage::SizeType _regionSize;

        int _numberOfSamples;
        std::vector<IndexType> _indexes;
        std::vector<PointType> _points;
    };

    NeighborSampler::NeighborSampler(RegionType region, RealImage::Pointer image) {
        setSampleRegion(region, image);
    }

    NeighborSampler::~NeighborSampler() {

    }

    void NeighborSampler::sampleValues(LinearImageInterpolatorType* interp, Px& px, Px& pz, double* pxValue) {
        // sample intensity values from images
        for (int k = 0; k < _numberOfSamples; k++) {
            PointType samplePoint = _points[k];
            fordim (l) {
                // translate to the particle
                samplePoint[l] += pz.x[l];
            }
            if (interp->IsInsideBuffer(samplePoint)) {
                pxValue[k] = interp->Evaluate(samplePoint);
            } else {
                IntIndex idx;
                interp->GetInputImage()->TransformPhysicalPointToIndex(samplePoint, idx);
                cout << "Out of region: " << idx << "; " << px.x[0] << ", " << px.x[1] << endl;
                pxValue[k] = 0;
            }
        }
    }

    void NeighborSampler::sampleGradients(GradientInterpolatorType* interp, Px& px, Px& pz, Px::Vector& pxGrad) {
        // sample intensity values from images
        for (int k = 0; k < _numberOfSamples; k++) {
            PointType samplePoint = _points[k];
            fordim (l) {
                // translate to the particle
                samplePoint[l] += pz[l];
            }
            GradientPixel pix = interp->Evaluate(samplePoint);
            fordim (l) {
                pxGrad[k][l] = pix[l];
            }
        }
    }

    // should have s sample index
    void NeighborSampler::setSampleRegion(RegionType& region, RealImage* reference) {
        _regionSize = region.GetSize();

        _indexes.clear();
        _points.clear();

        IndexType startIdx = region.GetIndex();
        if (RegionType::ImageDimension == 2) {
            createSampleIndexes2(startIdx);
        } else if (RegionType::ImageDimension == 3){
            createSampleIndexes3(startIdx);
        }
        _numberOfSamples = _indexes.size();

        _points.resize(_numberOfSamples);
        for (int i = 0; i < _numberOfSamples; i++) {
            reference->TransformIndexToPhysicalPoint(_indexes[i], _points[i]);
        }
    }

    void NeighborSampler::createSampleIndexes2(IndexType& startIdx) {
        IndexType idx = startIdx;
        _indexes.reserve(_regionSize[1]*_regionSize[0]);
        for (int j = 0; j < _regionSize[1]; j++) {
            idx[0] = startIdx[0];
            for (int i = 0; i < _regionSize[0]; i++) {
                _indexes.push_back(idx);
                idx[0] ++;
            }
            idx[1] ++;
        }
    }

    void NeighborSampler::createSampleIndexes3(IndexType& startIdx) {
        IndexType idx = startIdx;
        for (int k = 0; k < _regionSize[2]; k++) {
            idx[1] = startIdx[1];
            for (int j = 0; j < _regionSize[1]; j++) {
                idx[0] = startIdx[0];
                for (int i = 0; i < _regionSize[0]; i++) {
                    _indexes.push_back(idx);
                    idx[0] ++;
                }
                idx[1] ++;
            }
            idx[2] ++;
        }
    }

#pragma mark PxImageTerm implementations
    void PxImageTerm::computePixelEntropy(int i, int npx, int nsx, int nex, PxSubj::Vector &subjs, NeighborSampler* sampler) {

        /// ---
        /// ### The entropy computation for image patch
        /// * Prepare the entropy computer
        EntropyComputer<double> comp(nsx, npx, 1);
        comp.dataIter.FirstData();

        /// * Iterate every subject with *j*
        ///   * Sample intensity values
        for (int j = 0; j < nsx; j++) {
            sampler->sampleValues(subjs[j].image.imageSampler, subjs[j].particles[i], subjs[j].affineAligned[i], comp.dataIter.sample);
            comp.dataIter.NextData();
        }

        /// * Move to the center of the intensity values \f$\bar{X} = \sum X_i \f$
        comp.MoveToCenter();

        /// * Compute covariance matrix
        comp.ComputeCovariance();

        /// * Compute inverse covariance
        comp.ComputeGradient();

        /// * Compute image gradient
        Px::Vector imageGradient(nex);

        /// * Iterate every subject with *j*
        ///   * Sample the gradient of each image by the chain rule
        ///   * Compute the gradient of entropy with respect to a particle \f$dH/dx\f$
        for (int j = 0; j < nsx; j++) {
            sampler->sampleGradients(subjs[j].image.gradientSampler, subjs[j].particles[i], subjs[j].affineAligned[i], imageGradient);
            const double coeff = global->imageCoeff[subjs[j].attrs[i].label];
            fordim (d) {
                subjs[j].imageForces[i][d] = 0;
                for (int k = 0; k < nex; k++) {
                    subjs[j].imageForces[i][d] += comp.gradient[j][k] * imageGradient[k][d];
                }
                subjs[j].imageForces[i][d] *= coeff;
            }
        }
    }

    void PxImageTerm::computeImageTerm(PxSubj::Vector& subjs, int w) {
        /// ---
        const int npx = global->totalParticles;
        const int nsx = global->nsubjs;
        const int nex = (__Dim == 2) ? w * w : w * w * w;


        fordim (j) {
            patchRegion.SetIndex(j, -w/2);
            patchRegion.SetSize(j, w);
        }

        // sample intensity and gradients
        NeighborSampler sampler(patchRegion, subjs[0].image.image);

        /// Iterate every particle with *i*, and call computePixelEntropy()
        for (int i = 0; i < npx; i++) {
            computePixelEntropy(i, npx, nsx, nex, subjs, &sampler);
        }
    }

#pragma mark PxSystem Implementations
    void PxGlobal::load(ConfigFile& config) {
        nsubjs = config["particles.number-of-subjects"];
        nlabels = config["particles.number-of-labels"];

        useLocalRepulsion = false;
        config["particles"].lookupValue("use-local-repulsion", useLocalRepulsion);

        useAffineTransform = false;
        config["particles"].lookupValue("use-local-repulsion", useAffineTransform);

        useEnsembleForce = true;
        config["particles"].lookupValue("use-ensemble-force", useEnsembleForce);
        Setting& settingNumParticles = config["particles.number-of-particles"];
        for (int i = 0; i < settingNumParticles.getLength(); i++) {
            numParticles.push_back(settingNumParticles[i]);
        }
        totalParticles = std::accumulate(numParticles.begin(), numParticles.end(), 0);

        Setting& params = config["particles.parameters"];
        repulsionCoeff.resize(nlabels);
        ensembleCoeff.resize(nlabels);
        sigmaParams.resize(nlabels);
        cutoffParams.resize(nlabels);
        imageCoeff.resize(nlabels);
        for (int i = 0; i < params.getLength(); i++) {
            repulsionCoeff[i] = params[i][0];
            ensembleCoeff[i] = params[i][1];
            sigmaParams[i] = params[i][2];
            cutoffParams[i] = params[i][3];
            imageCoeff[i] = params[i][4];
        }
    }

    void PxSystem::setupNeighbors(const IntVector& numPx, PxSubj& subj) {
        global.neighbors.resize(subj.size());
        assert(numPx.size() == global.nlabels);

        // construct delaunay triangulations
        // implicitly assume that the index of numPx represents the label index
        for (int i = 0; i < numPx.size(); i++) {
            ParticleMesh mesh;
            mesh.constructNeighbors(i, numPx[i], subj, global.cutoffParams[i], global.neighbors);
        }
    }


//    void PxSystem::clearForces() {
//        for (int i = 0; i < subjs.size(); i++) {
//            std::fill(subjs[i].forces.begin(), subjs[i].forces.end(), 0);
//            std::fill(subjs[i].ensembleForces.begin(), subjs[i].ensembleForces.end(), 0);
//            std::fill(subjs[i].imageForces.begin(), subjs[i].imageForces.end(), 0);
//        }
//    }

    // load data for particle system
    void PxSystem::loadSystem(ConfigFile& config) {
        global.load(config);

        sampler.setGlobalConfig(&global);
        subjs.resize(global.nsubjs);
        Setting& subjconfig = config["particles.subjects"];
        for (int i = 0; i < global.nsubjs; i++) {
            subjs[i].setGlobalConfig(&global);
            Setting& subjdata = subjconfig[i];
            subjs[i].regions.resize(global.nlabels);
            for (int j = 0; j < global.nlabels; j++) {
                subjs[i].regions[j].load(subjdata["labels"][j]);
            }
            subjs[i].image.load(subjdata["images"][0]);
        }
    }


    // iterate subjects to load particle data
    // only when ignore-particle-input = false
    bool PxSystem::loadParticles(ConfigFile& config) {
        if (config["particles.ignore-particle-input"]) {
            return false;
        }

        Setting& subjconfig = config["particles.subjects"];

        bool ok = true;
        for (int i = 0; i < subjconfig.getLength() && ok; i++) {
            string f = subjconfig[i]["particle-input"];
            ok = subjs[i].load(f);
        }

        return ok;
    }

    // initialize sampler on first run
    bool PxSystem::loadSampler(ConfigFile& config) {
        // sampler needs to load region information
        // if sampler-input doesn't exist!!
        if (!config["particles.ignore-sampler-input"]) {
            string samplerCache = config["particles.sampler-cache"];
            if (sys.checkFile(samplerCache)) {
                if (sampler.load(samplerCache)) {
                    return true;
                }
            }
        }

        // main configuration setting
        Setting& mainConfig = config["particles"];

        // load each label information and create intersection
        Setting& regionList = config["particles.sampler.labels"];
        if (global.nlabels != regionList.getLength()) {
            cout << "sampler has different number of labels" << endl;
            exit(0);
        }

        /// global configuration setting
        sampler.setGlobalConfig(&global);

        // sampler region load
        sampler.regions.resize(global.nlabels);
        for (int i = 0; i < regionList.getLength(); i++) {
            // if there exists intersection cache
            if (!sampler.regions[i].load(regionList[i], false)) {
                // allocate new image buffer
                LabelImage::Pointer intersection = labelIO.NewImage(subjs[0].regions[i].labelmap);

                int nPixels = intersection->GetPixelContainer()->Size();
                LabelImage::PixelType* outputBuff = intersection->GetBufferPointer();

                // initialize label pointers for fast access
                std::vector<LabelImage::PixelType*> ptrs;
                ptrs.resize(global.nsubjs);

                for (int j = 0; j < global.nsubjs; j++) {
                    ptrs[j] = subjs[j].regions[i].labelmap->GetBufferPointer();
                }

                // loop over pixels to check if the pixel is intersection
                for (int ii = 0; ii < nPixels; ii++) {
                    // test if every labels have label
                    bool all = true;
                    for (unsigned int j = 0; all && j < ptrs.size(); j++) {
                        all = all && (*ptrs[j] > 0);
                    }
                    if (all) {
                        *outputBuff = 1;
                    }

                    // increment pixel pointer
                    for (unsigned int j = 0; j < ptrs.size(); j++) {
                        ptrs[j] ++;
                    }
                    outputBuff ++;
                }
                sampler.regions[i].labelmap = intersection;
                sampler.regions[i].deriveImages(regionList[i], intersection, true);
            }
        }

        // sample particles with no cache
        sampler.sampleParticles(global.numParticles);

        // save cache
        string file = "";
        if (mainConfig.lookupValue("sampler-initial-cache", file)) {
            ofstream of(file.c_str());
            sampler.save(of);
            of.close();
        }
        return false;
    }

    // iterate subjects to save particle data
    bool PxSystem::saveParticles(ConfigFile& config, string outputName) {
        Setting& subjconfig = config["particles.subjects"];

        bool ok = true;
        for (int i = 0; i < global.nsubjs && ok; i++) {
            string f = subjconfig[i][outputName];
            ofstream of(f.c_str());
            subjs[i].save(of);
            of.close();
        }
        return ok;
    }

    // duplicate particles from sampler to each subject
    void PxSystem::transferParticles() {
        for (int i = 0; i < global.nsubjs; i++) {
            subjs[i].setGlobalConfig(&global);
            subjs[i].resize(sampler.size());
            subjs[i].attrs = sampler.attrs;
            subjs[i].particles = sampler.particles;
        }
    }

    /**
     * 1) load configuration
     * 2) load particle data
     * 3) load sampler and run preprocessing if necessary
     * 4) save initialized particles
     */
    void PxSystem::initialize(Options& opts, StringVector& args) {
        _config.load(opts.GetConfigFile());
        loadSystem(_config);
        if (!loadParticles(_config)) {
            if (!loadSampler(_config)) {
                initialLoop();
                string samplerCache = _config["particles.sampler-cache"];
                ofstream of(samplerCache.c_str());
                sampler.save(of);
                of.close();
            }
            transferParticles();
            saveParticles(_config, "particle-input");
        }
    }

    // initial loop to distribute particles inside region
    void PxSystem::initialLoop() {
        double t0 = _config["particles.sampler-time-steps.[0]"];
        double dt = _config["particles.sampler-time-steps.[1]"];
        double t1 = _config["particles.sampler-time-steps.[2]"];

        ImageIO<LabelImage3> io;
        LabelImage::Pointer refImage = sampler.regions[0].labelmap;
        LabelImage3::Pointer tracker = CreateImage3(refImage, (t1 - t0) / dt + 1);

        cout << sampler << endl;

        int m = 0;
        for (double t = t0; t <= t1; t += dt, m++) {
            cout << "t = " << t << endl;
            sampler.clearForce();
            sampler.constrainParticles();

            // how to set repulsion paramter per region?
            sampler.computeRepulsion(true);
            sampler.constrainForces();
            sampler.updateSystem(dt);

            for (int i = 0; i < sampler.particles.size(); i++) {
                LabelImage::IndexType idx = sampler.getIndex(i);
                LabelImage3::IndexType tx;
                fordim (k) {
                    tx[k] = idx[k];
                }
                tx[2] = m;
                if (tracker->GetBufferedRegion().IsInside(tx)) {
                    tracker->SetPixel(tx, 1);
                }
            }
        }
        cout << sampler << endl;

//        printAdjacencyMatrix();

        io.WriteImage("initialLoop.nrrd", tracker);
    }


    /**
     * registration loop
     *
     */
    void PxSystem::loop() {
        double t0 = _config["particles.time-steps.[0]"];
        double dt = _config["particles.time-steps.[1]"];
        double t1 = _config["particles.time-steps.[2]"];

        ensemble.global = &global;

        ImageIO<LabelImage3> io;
        LabelImage::Pointer refImage = subjs[0].regions[0].labelmap;
        std::vector<LabelImage3::Pointer> trackers;

        for (int i = 0; i < global.nsubjs; i++) {
            LabelImage3::Pointer tracker = CreateImage3(refImage, (t1 - t0) / dt + 1);
            trackers.push_back(tracker);
        }


        setupNeighbors(global.numParticles, subjs[0]);

        int m = 0;
        int s = 0;
        for (double t = t0; t <= t1; t += dt, m++) {
            cout << "t = " << t << endl;

            if (m % 10 == 0) {
                // construct neighbors
                setupNeighbors(global.numParticles, subjs[++s%2]);
            }

            // compute internal forces
            for (int i = 0; i < global.nsubjs; i++) {
                subjs[i].clearForce();
                subjs[i].constrainParticles();
                // how to set repulsion paramter per region?
                subjs[i].computeRepulsion();
            }

            // apply least-square based affine transform
            affineTransformParticles();

            // now compute ensemble and intensity forces
            if (global.useEnsembleForce) {
                ensemble.computeAttraction(subjs);
            }


            // compute pixel entropy
            
//            cout << subjs[0].ensembleForces << endl;

            // update system
            for (int i = 0; i < global.nsubjs; i++ ) {
                subjs[i].constrainForces();
                subjs[i].updateSystem(dt);
                subjs[i].constrainParticles();
            }



            // mark each particle location at a tracker image
            for (int i = 0; i < global.nsubjs; i++) {
                for (int j = 0; j < subjs[0].particles.size(); j++) {
                    LabelImage::IndexType idx = subjs[i].getIndex(j);
                    LabelImage3::IndexType tx;
                    fordim (k) {
                        tx[k] = idx[k];
                    }
                    tx[2] = m;
                    if (trackers[i]->GetBufferedRegion().IsInside(tx)) {
                        trackers[i]->SetPixel(tx, 1);
                    }
                }
            }
        }


        for (int i = 0; i < global.nsubjs; i++) {
            string file = _config["particles.particle-images"][i];
            io.WriteImage(file, trackers[i]);
        }
    }


    void PxSystem::affineTransformParticles() {
        // copy from particle space to affineAligned space
        for (int i = 0; i < global.nsubjs; i++) {
            std::copy(subjs[i].particles.begin(), subjs[i].particles.end(), subjs[i].affineAligned.begin());
        }

        // estimate affine transform only if it is allowed
        if (global.useAffineTransform) {
            subjs[1].affineTxf.estimateAffineTransform(&subjs[0]);
        }
    }

    /**
     * warp label images in according to its setting
     * setting = [ srcIdx, dstIdx, inputFile, outputFile ]
     * control point spacing and other bspline parameters may be cached to the system
     */
    void PxSystem::warpLabels(libconfig::Setting& warpedLabels) {
        for (int i = 0; i < warpedLabels.getLength(); i++) {
            int src = warpedLabels[i][0];
            int dst = warpedLabels[i][1];
            string input = warpedLabels[i][2];
            string output = warpedLabels[i][3];

            // if src-index and dst-index are not correct
            if (src < 0 || src >= global.nsubjs) {
                return;
            }
            if (dst < 0 || dst >= global.nsubjs) {
                return;
            }

            cout << "Warping: " << src << " => " << dst << " : " << input << " => " << output << endl;
            LabelImage::Pointer inputImage = labelIO.ReadCastedImage(input);
            ParticleWarp warp;
            warp.setParameters(_config);
            warp.reference = inputImage;
            warp.estimateBsplineWarp(subjs[src].particles, subjs[dst].particles);
            LabelImage::Pointer outputImage = warp.warpLabel(inputImage);
            labelIO.WriteImage(output, outputImage);
        }
    }

    /**
     * warp intensity images in according to its setting
     * setting = [ srcIdx, dstIdx, inputFile, outputFile ]
     * control point spacing and other bspline parameters may be cached to the system
     */
    void PxSystem::warpImages(libconfig::Setting& data) {
        for (int i = 0; i < data.getLength(); i++) {
            int src = data[i][0];
            int dst = data[i][1];
            string input = data[i][2];
            string output = data[i][3];

            // if src-index and dst-index are not correct
            if (src < 0 || src >= global.nsubjs) {
                return;
            }
            if (dst < 0 || dst >= global.nsubjs) {
                return;
            }

            // warp output
            cout << "Warping: " << src << " => " << dst << " : " << input << " => " << output << endl;
            RealImage::Pointer inputImage = io.ReadCastedImage(input);
            ParticleWarp warp;
            warp.reference = subjs[0].regions[0].labelmap;
            warp.estimateBsplineWarp(subjs[dst].particles, subjs[src].particles);
            RealImage::Pointer outputImage = warp.warpImage(inputImage);
            io.WriteImage(output, outputImage);
        }
    }

    void PxSystem::main(pi::Options &opts, StringVector &args) {
        initialize(opts, args);
        loop();
        saveParticles(_config, "particle-output");

        warpLabels(_config["particles/warped-labels"]);
        warpImages(_config["particles/warped-images"]);

        print();
        saveAdjacencyMatrix("adj.txt");
    }

    void PxSystem::print() {
        for (int i = 0; i < sampler.particles.size(); i++) {
            cout << sampler.particles[i] << endl;
        }
    }

    void PxSystem::saveAdjacencyMatrix(std::string file) {
        ofstream of(file.c_str());

        vnl_matrix<int> adjmat(global.totalParticles, global.totalParticles);
        adjmat.fill(0);

        for (int i = 0; i < global.neighbors.size(); i++) {
            for (int j = 0; j < global.neighbors[i].size(); j++) {
                adjmat[i][ global.neighbors[i][j] ] = 1;
            }
        }
        of << adjmat;
        of.close();
    }

#pragma mark ParticleRunner implementations
    void ParticleRunner::main(pi::Options &opts, StringVector &args) {
        try {
            if (opts.GetBool("--help")) {
                printHelp();
                return;
            }
            _system.main(opts, args);
        } catch (ParseException & ex) {
            cout << "Parsing Error" << endl;
            cout << "File: " << ex.getFile() << endl;
            cout << "Line: " << ex.getLine() << endl;
            cout << "Error: " << ex.getError() << endl;
        } catch (SettingNotFoundException& ex) {
            cout << "Setting Not Found Error" << endl;
            cout << "Path: " << ex.getPath() << endl;
            cout << "What: " << ex.what() << endl;
        } catch (SettingTypeException& ex) {
            cout << "Setting Type Error" << endl;
            cout << "Path: " << ex.getPath() << endl;
        } catch (FileIOException& ex) {
            cout << "File I/O Exception" << endl;
            cout << "Path: " << ex.what() << endl;
        }
    }

    void ParticleRunner::printHelp() {
        cout << "particle image registration: built (" << __DATE__ << ")" << endl;
    }
}