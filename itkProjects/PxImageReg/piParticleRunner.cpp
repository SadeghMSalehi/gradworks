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

//    void executeParticleRunner(pi::Options &opts, StringVector &args);
//    void executeEntropyImage(Options& opts, StringVector& args);

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



#pragma mark PxGlobal Implementation
    void PxGlobal::setNumberOfLabels(int nLabels) {
        numParticles.resize(nLabels);
        cutoffParams.resize(nLabels);
        sigmaParams.resize(nLabels);
        repulsionCoeff.resize(nLabels);
        ensembleCoeff.resize(nLabels);
        imageCoeff.resize(nLabels);
    }

#pragma mark PxI implementation

    PxI::PxI() {
        imageSampler = LinearImageInterpolatorType::New();
        gradientSampler = GradientInterpolatorType::New();
    }

    void PxI::load(std::string file) {
        /// ---
        /// If the file exists and is loaded correctily, set the image
        RealImage::Pointer image = io.ReadCastedImage(file);
        if (image.IsNull()) {
            cout << "can't read image: " << file << endl;
            exit(0);
        }
        setImage(image);
    }


    void PxI::setImage(RealImage::Pointer image) {
        /// ---
        this->image = image;
        /// * The gradient is computed after applying the Gaussian filter with the sigma of the half of the spacing.
        this->gradient = ComputeGaussianGradient(image, image->GetSpacing()[0] / 2.0);

        this->imageSampler->SetInputImage(this->image);
        this->gradientSampler->SetInputImage(this->gradient);
    }


    void PxI::write(std::string file) {
        ImageIO<RealImage> io;
        io.WriteImage(file, this->image);
    }

    void PxI::writeGradient(std::string file) {
        ImageIO<RealImage> io;
        io.WriteImage(file.c_str(), ComputeGradientMagnitude(gradient));
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
            distmap = ComputeDistanceMap(label, "");
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
        GradientImageInterpolatorType::PointType pi;
        fordim (k) {
            pi[k] = p[k];
        }
        GradientImageInterpolatorType::OutputType nx = gradIntp->Evaluate(pi);
        fordim (k) {
            nOut[k] = nx[k];
        }
    }

    bool PxR::isIn(Px& p, bool& isInsideBuffer) {
        LabelImageInterpolatorType::PointType pi;
        fordim (k) {
            pi[k] = p[k];
        }

        isInsideBuffer = lablIntp->IsInsideBuffer(pi);
        if (isInsideBuffer && lablIntp->Evaluate(pi) > 0) {
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

        bool isInsideBuffer = false;
        if (isIn(p, isInsideBuffer)) {
            return false;
        }
//        LinearVectorImageInterpolatorType::InputType pi;

        LabelImage::IndexType newIx = ix;
        LinearVectorImageInterpolatorType::OutputType dx = distIntp->EvaluateAtIndex(ix);
        fordim (k) {
            newIx[k] += dx[k];
        }
        labelmap->TransformIndexToPhysicalPoint(newIx, px);
        fordim (k) {
            p[k] = px[k];
        }
        return true;
    }

    bool PxR::normalForce(Px& p, Px& f, Px& fout) {
        bool isInsideBuffer = false;
        if (isIn(p, isInsideBuffer)) {
            return false;
        }

        if (!isInsideBuffer) {
            LabelImage::IndexType pi;
            LabelImage::PointType px;
            fordim (k) {
                px[k] = p[k];
            }
            cout << this->labelmap->GetBufferedRegion().GetSize() << endl;
            this->labelmap->TransformPhysicalPointToIndex(px, pi);
            cout << "A particle is located out of the image buffer: " << p << "; " << pi << endl;
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
    void PxSubj::computeSamplingTerm(const bool isInitial, bool ignoreNeighbors) {
        const int npx = size();
        VNLMatrix weights(npx, npx);
        weights.fill(0);

        // loop over all particle pairs
        for (int i = 0; i < npx; i++) {
            Px& pi = particles[i];
            weights[i][i] = 0;
            double kappa = 1;

            const int label = attrs[i].label;
            double sigma = 0;
            double cutoff = 0;
            if (isInitial) {
                sigma = global->initialSigmaParams[label];
                cutoff = global->initialCutoffParams[label];
            } else {
                sigma = global->sigmaParams[label];
                cutoff = global->cutoffParams[label];
            }

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


    void PxSubj::updateSystem(double dt) {
        /// ---
        Px::Vector::iterator p = particles.begin();
        Px::Vector::iterator f = forces.begin();
        Px::Vector::iterator e = ensembleForces.begin();
        Px::Vector::iterator i = imageForces.begin();

        for (int j = 0; p != particles.end(); p++, f++, e++, i++, j++) {
            fordim (k) {
                p->x[k] += (dt * (f->x[k] + e->x[k] + i->x[k]));
                assert(!std::isnan(p->x[k]));
            }
            if (j == 0) {
                cout << (*p) << "; " << (*f) << "; " << (*e) << "; " << (*i) << endl;
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
        for (int i = 0; f != forces.end(); p++, a++, f++, i++) {
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
                fordim (d) {
                    particles[total].x[d] = pts[d];
                }

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

        /// @brief sample image intensities from a given interpolator
        /// @param interp A image interpolator
        /// @param px A particle to sample
        /// @param pxValue Output variable for the sampled intensity vectors
        void sampleValues(LinearImageInterpolatorType* interp, Px& px,  double* pxValue);

        /// @brief Sample image gradients of a local patch neary by a particle *px*
        /// @param interp A gradient interpolator
        /// @param px A particle to sample
        /// @param pxGrad Output variable for the sampled gradient vectors
        void sampleGradients(GradientInterpolatorType* interp, Px& px, Px::Vector& pxGrad);

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

    void NeighborSampler::sampleValues(LinearImageInterpolatorType* interp, Px& px, double* pxValue) {
        /// ---
        /// * Iterate over sample points with *k*
        for (int k = 0; k < _numberOfSamples; k++) {
            /// * For each displacement vector at *px*
            PointType samplePoint = _points[k];
            fordim (l) {
                // translate to the particle
                samplePoint[l] += px.x[l];
            }

            /// * Evaluate the intensity; This doesn't check if the point is inside. Should be careful if the patch is fully inside an image buffer.
            if (true) {
                if (interp->IsInsideBuffer(samplePoint)) {
                    pxValue[k] = interp->Evaluate(samplePoint);
                    if (pxValue[k] > 1e5) {
                        cout << "Error" << endl;
                    }
                } else {
                    /// * If the point is outside a buffer, then set to zero.
                    IntIndex idx;
                    interp->GetInputImage()->TransformPhysicalPointToIndex(samplePoint, idx);
                    cout << "Out of region: " << idx << "; " << px.x[0] << ", " << px.x[1] << endl;
                    pxValue[k] = 0;
                }
            }
        }
    }

    void NeighborSampler::sampleGradients(GradientInterpolatorType* interp, Px& px, Px::Vector& pxGrad) {
        /// ---
        /// * For each sampling point, sample the gradient values from images
        /// * This does not check if the point is inside a buffer. Please be careful to use
        for (int k = 0; k < _numberOfSamples; k++) {
            PointType samplePoint = _points[k];
            fordim (l) {
                // translate to the particle
                samplePoint[l] += px[l];
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
    double PxImageTerm::computePixelEntropy(int i, int npx, int nsx, int nex, PxSubj::Vector &subjs, NeighborSampler* sampler) {

        /// ---
        /// ### The entropy computation for image patch
        /// * Prepare the entropy computer _(nSubj x nPatchElems)_
        EntropyComputer<double> comp(nsx, nex, 1);
        comp.dataIter.FirstData();

        /// * Iterate every subject with *j*
        ///   * Sample intensity values; The intensity sampled from the warped image
        for (int j = 0; j < nsx; j++) {
            sampler->sampleValues(subjs[j].warpedImage.imageSampler, subjs[j].particles[i], comp.dataIter.sample);
            comp.dataIter.NextData();
        }

        /// * Move to the center of the intensity values \f$\bar{X} = \sum X_i \f$
        comp.MoveToCenter();

        /// * Compute covariance matrix
        comp.ComputeCovariance();

        /// * Compute inverse covariance
        comp.ComputeGradient();


        /// * Compute the entropy for this particle
        const double entropy = comp.ComputeEntropy();

        /// * Compute image gradient
        Px::Vector imageGradient(nex);

        /// * Iterate every subject with *j*
        ///   * Sample the gradient of each image by the chain rule
        ///   * Compute the gradient of entropy with respect to a particle \f$dH/dx\f$
        for (int j = 0; j < nsx; j++) {
            sampler->sampleGradients(subjs[j].warpedImage.gradientSampler, subjs[j].particles[i], imageGradient);
            const double coeff = global->imageCoeff[subjs[j].attrs[i].label];
            fordim (d) {
                subjs[j].imageForces[i][d] = 0;
                for (int k = 0; k < nex; k++) {
                    subjs[j].imageForces[i][d] += comp.gradient[j][k] * imageGradient[k][d];
                }
                subjs[j].imageForces[i][d] *= coeff;
            }
        }

        return entropy;
    }

    void PxImageTerm::computeImageTerm(PxSubj::Vector& subjs, int w, DoubleVector& entropy) {
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
            entropy[i] = computePixelEntropy(i, npx, nsx, nex, subjs, &sampler);
        }
    }


#pragma mark PxBsplineDeformation Implementation
    void PxBsplineDeformation::computeDeformationToAverage(PxSubj::Vector &subjs) {
        /// ---
        /// 1. Compute the average of particles
        Px::Vector average;
        average.resize(global->totalParticles);
        for (int j = 0; j < global->totalParticles; j++) {
            forfill(average[j].x, 0);
            for (int i = 0; i < global->nsubjs; i++) {
                fordim(d) {
                    average[j].x[d] += subjs[i].particles[j].x[d];
                }
            }
            fordim (d) {
                average[j].x[d] /= global->nsubjs;
            }
        }

        for (int i = 0; i < global->nsubjs; i++) {
            /// 2. Compute the warp using ParticleWarp
            ParticleWarp warp;
            warp.controlSpacing = global->controlPointSpacing;
            warp.reference = global->referenceImage;
            warp.estimateBsplineWarp(subjs[i].particles, average);

            /// 3. Compute the deformed images and its gradient
            subjs[i].warpedImage.setImage(warp.warpImage(subjs[i].image.image));
        }

        const bool debugWarpedImages = false;
        if (debugWarpedImages) {
            subjs[0].warpedImage.write("warpedImage.0.nrrd");
            subjs[0].warpedImage.writeGradient("gradientImage.0.nrrd");
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

        useImageTerm = true;
        config["particles"].lookupValue("use-image-term", useImageTerm);

        Setting& settingNumParticles = config["particles.number-of-particles"];
        for (int i = 0; i < settingNumParticles.getLength(); i++) {
            numParticles.push_back(settingNumParticles[i]);
        }
        totalParticles = std::accumulate(numParticles.begin(), numParticles.end(), 0);


        /// Set parameters related optimization function
        /// - particles.parameters
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

        /// Set parameters related to the initial sampling process
        /// - particles.initial-parameters
        Setting& initialParams = config["particles.initial-parameters"];
        initialSigmaParams.resize(nlabels);
        initialCutoffParams.resize(nlabels);
        for (int i = 0; i < initialParams.getLength(); i++) {
            initialSigmaParams[i] = initialParams[i][0];
            initialCutoffParams[i] = initialParams[i][1];
        }


        /// Set B-spline grid related parameters
        /// - particles.referenceImage
        /// - particles.controlPointSpacing
        ImageIO<LabelImage> io;
        referenceImage = io.ReadCastedImage(config["particles.bspline-reference-grid"]);
        controlPointSpacing = config["particles.bspline-control-point-spacing"];
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

    void PxSystem::loadSystem(ConfigFile& config) {
        /// ---
        /// 1. Load global parameters
        global.load(config);

        /// 2. Initialization setup
        sampler.setGlobalConfig(&global);
        subjs.resize(global.nsubjs);

        /// 3. Load subject parameters
        Setting& subjconfig = config["particles.subjects"];
        for (int i = 0; i < global.nsubjs; i++) {
            subjs[i].setGlobalConfig(&global);
            Setting& subjdata = subjconfig[i];
            subjs[i].regions.resize(global.nlabels);
            /// 4. Load each label
            for (int j = 0; j < global.nlabels; j++) {
                subjs[i].regions[j].load(subjdata["labels"][j]);
            }
            /// 5. Load a intensity image
            cout << "Loading image ..." << endl;
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

    void PxSystem::initialLoop() {
        /// The initialization is important to place particles near optimal. Assuming affine registrations are performed all subjects, each particle should be placed near the optimal only bearing non-linear deformation. This sampling step is iterated with the time steps given as __particles.initial-time-step__.
        double t0 = _config["particles.initial-time-steps.[0]"];
        double dt = _config["particles.initial-time-steps.[1]"];
        double t1 = _config["particles.initial-time-steps.[2]"];

        ImageIO<LabelImage3> io;
        LabelImage::Pointer refImage = sampler.regions[0].labelmap;
        LabelImage3::Pointer tracker = CreateImage3(refImage, (t1 - t0) / dt + 1);

        cout << sampler << endl;

        int m = 0;
        /// #### Iteration process
        for (double t = t0; t <= t1; t += dt, m++) {
            cout << "t = " << t << endl;
            /// - Clear previous forces
            sampler.clearForce();
            /// - Move particles back inside the region
            sampler.constrainParticles();

            /// - Compute the sampling term
            sampler.computeSamplingTerm(true, false);
            sampler.constrainForces();
            sampler.updateSystem(dt);

            if (__Dim == 2) {
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
        }
        cout << sampler << endl;

//        printAdjacencyMatrix();

        io.WriteImage("initialLoop.nrrd", tracker);
    }



    void PxSystem::loop() {
        /// Perform the optimization from *t0* to *t1* with the timestep *dt* as given in the configuration file at __particles.time-steps.[]__. The current implementation assumes regestering two dimensional images, and the particles at each time step is written into 3D volumes given in __particles.particle-images__. For the construction of neighborhood information, the adjacency matrix is constructed for every 10 time steps.
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

            /// At each time step, it computes the gradient of _1) the sampling term, 2) the ensemble term, and 3) image term_.
            /// #### Sampling Term
            /// - Iterative every subject
            /// - Clear the gradient values
            /// - Move particles outside back to the inside of the region
            /// - Evaluate the sampling gradient
            for (int i = 0; i < global.nsubjs; i++) {
                subjs[i].clearForce();
                subjs[i].constrainParticles();
                subjs[i].computeSamplingTerm(false);
            }

            /// Then, apply least-square based affine transform
            affineTransformParticles();

            /// #### Ensemble Term
            /// The ensemble term is evaluated by PxEnsemble.
            if (global.useEnsembleForce) {
                ensemble.computeAttraction(subjs);
            }

            if (global.useImageTerm) {
                /// #### Image Deformation
                /// The deformation is performed by PxBsplineDeformation. This image term can be turned off using 'particles.use-image-term = false'.
                PxBsplineDeformation deformation(&global);
                deformation.computeDeformationToAverage(subjs);


                /// #### Image Term
                /// The image iterm is evaluated by PxImageTerm after the B-spline deformation.
                PxImageTerm imageTerm(&global);
                DoubleVector entropyVector;
                entropyVector.resize(global.totalParticles);
                imageTerm.computeImageTerm(subjs, 17, entropyVector);
            }

            /// #### Update the system
            /// The gradient evaluated is added with different weighting.
            /// - Update forces to direct the inside
            /// - Update the particle position
            /// - Update particles outside into the inside
            for (int i = 0; i < global.nsubjs; i++ ) {
                subjs[i].constrainForces();
                subjs[i].updateSystem(dt);
                subjs[i].constrainParticles();
            }

            /// #### Debugging
            /// Mark each particle location at a tracker image
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
            warp.controlSpacing = global.controlPointSpacing;
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


#pragma mark --entropyImage implementation
    void executeEntropyImage(Options& opts, StringVector& args) {
        /// * If it is in __test__ mode (`--test`), it automatically provides a set of test files.
        bool testMode = opts.GetBool("--test");
        string output = opts.GetString("-o");

        if (testMode) {
            args.clear();
            args.resize(4);
            args[0] = "/NIRAL/work/joohwi/c57/SliceData/07_slice.nrrd";
            args[1] = "/NIRAL/work/joohwi/c57/SliceData/11_slice.nrrd";
            args[2] = "/NIRAL/work/joohwi/c57/SliceData/12_slice.nrrd";
            args[3] = "/NIRAL/work/joohwi/c57/SliceData/15_slice.nrrd";
            output = "/NIRAL/work/joohwi/c57/entropy_test.nrrd";
        } else {
            if (output == "" || args.size() < 2) {
                cout << "Error: the output is not given or at least two images are not given." << endl;
                exit(0);
            }
        }

        /// * Create a vector of PxSubj list from arguments to feed into PxImageTerm::computeImageTerm()
        PxGlobal global;
        global.nsubjs = args.size();
        global.setNumberOfLabels(1);

        PxSubj::Vector subjs;
        subjs.resize(args.size());

        /// * Compute the image region by shrinking the region boundary
        // perform only for the first subject
        RealImage::RegionType region;
        for (int i = 0; i < args.size(); i++) {
            subjs[i].setGlobalConfig(&global);
            subjs[i].image.load(args[i]);

            /// * Create a list of particles to cover all voxels
            if (i == 0) {
                // shrink the region by the patch size
                region = subjs[i].image.image->GetBufferedRegion();
                region.ShrinkByRadius(5);

                global.totalParticles = region.GetNumberOfPixels();
                subjs[i].resize(global.totalParticles);

                RealImageIteratorType iter(subjs[i].image.image, region);
                for (int j = 0; !iter.IsAtEnd(); ++iter, j++) {
                    IntIndex idx = iter.GetIndex();
                    ImagePoint point;
                    subjs[i].image.image->TransformIndexToPhysicalPoint(idx, point);
                    for (int d = 0; d < __Dim; d++) {
                        subjs[i].particles[j].x[d] = point[d];
                    }
                }

                cout << "total # of points: " << global.totalParticles << endl;
            } else {
                // duplicate particles
                subjs[i].resize(global.totalParticles);
                subjs[i].particles = subjs[0].particles;
            }
        }

        // Prepare to store entropy values.
        DoubleVector entropyValues;
        entropyValues.resize(global.totalParticles);

        /// * Create an instance of PxImageTerm and call PxImageTerm::computeImageTerm
        PxImageTerm imageTerm(&global);
        imageTerm.global = &global;
        imageTerm.computeImageTerm(subjs, 5, entropyValues);

        /// * Convert the entropyValues vector into an image after scaling
        ImageIO<RealImage> io;
        RealImage::Pointer outputImage = io.NewImage(subjs[0].image.image);

        double minEntropy = *(std::min_element(entropyValues.begin(), entropyValues.end()));
        double maxEntropy = *(std::max_element(entropyValues.begin(), entropyValues.end()));

        RealImageIteratorType iter(outputImage, region);
        for (int i = 0; !iter.IsAtEnd(); ++iter, ++i) {
            RealImage::PixelType value = (10000*(entropyValues[i] - minEntropy) / maxEntropy);
            outputImage->SetPixel(iter.GetIndex(), value);
        }

        /// * Write the output image
        io.WriteImage(output, outputImage);
    }
}