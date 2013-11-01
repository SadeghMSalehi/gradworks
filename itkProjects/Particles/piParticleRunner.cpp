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
#include "piParticleWarp.h"

#include <numeric>
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

    std::vector<PatchImage::Pointer> __patchImages;

    template <typename T>
    typename T::Pointer applyGaussian(typename T::Pointer input, double sigma) {
        typedef itk::SmoothingRecursiveGaussianImageFilter<T> GaussianFilter;
        typename GaussianFilter::Pointer filter = GaussianFilter::New();
        filter->SetSigma(sigma);
        filter->SetInput(input);
        filter->Update();
        typename T::Pointer output = filter->GetOutput();
        output->DisconnectPipeline();
        return output;
    }

    LabelImage3::Pointer createImage3(LabelImage::Pointer refImage, int m) {
        ImageIO<LabelImage3> io;
        LabelImage::SizeType sz = refImage->GetBufferedRegion().GetSize();
        LabelImage3::Pointer tracker = io.NewImageT(sz[0], sz[1], m);
        LabelImage3::SpacingType spacing;
        LabelImage3::PointType origin;
        LabelImage3::DirectionType direction;
        direction.Fill(0);
        fordim (k) {
            spacing[k] = refImage->GetSpacing()[k];
            origin[k] = refImage->GetOrigin()[k];
            fordim (l) {
                direction[k][l] = refImage->GetDirection()[k][l];
            }
        }
        spacing[2] = spacing[0];
        origin[2] = 0;
        direction[2][2] = 1;

        tracker->FillBuffer(0);
        tracker->SetSpacing(spacing);
        tracker->SetOrigin(origin);
        tracker->SetDirection(direction);

        return tracker;
    }

    void executeParticleRunner(pi::Options &opts, StringVector &args) {
        ParticleRunner runner;
        runner.main(opts, args);
    }


    // utility operator overloading
    ostream& operator<<(ostream& os, const Px& par) {
        fordim(k) { os << par.x[k] << " "; }
        return os;
    }
    ostream& operator<<(ostream& os, const Px::Vector& par) {
        Px::Vector::const_iterator p = par.begin();

        if (p == par.end()) {
            return os;
        }

        os << *p;
        for (p++; p != par.end(); p++) {
            os << ";" << *p;
        }

        return os;
    }
    ostream& operator<<(ostream& os, const PxSubj& p) {
        for (int i = 0; i < p.size(); i++) {
            cout << p.particles[i] << endl;
        }
        return os;
    }

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

        ImageProcessing proc;
        if (saveLabel) {
            labelmap = label;
            io.WriteImageS<LabelImage>(labelfile, label);
        }
        if (distmap.IsNull()) {
            distmap = proc.ComputeDistanceMap(label);
            io.WriteImageS<VectorImage>(distfile, distmap);
        }
        if (gradmap.IsNull()) {
            LabelImage::SpacingType spacing = label->GetSpacing();
            gradmap = proc.ComputeGaussianGradient(label, spacing[0] / 2.0);
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

    void PxSubj::computeRepulsion(double coeff, double sigma, double cutoff) {
        const double sigma2 = sigma * sigma;

        const int npx = size();
        VNLMatrix weights(npx, npx);

        // loop over all particle pairs
        for (int i = 0; i < npx; i++) {
            Px& pi = particles[i];
            weights[i][i] = 0;
            double kappa = 1;
            for (int j = i + 1; j < npx; j++) {
                Px& pj = particles[j];

                // compute weight
                double dij = pi.dist(pj);
                if (dij < cutoff) {
                    weights[i][j] = weights[j][i] = exp(-dij*dij*kappa/(sigma2)) / dij;
                } else {
                    weights[i][j] = weights[j][i] = 0;
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

        for (; p != particles.end() && f != forces.end(); p++, f++) {
            fordim (k) { p->x[k] += (dt * f->x[k]); }
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

        if (labels.getLength() != data[0].getLength()) {
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


    void PxSystem::createSampler() {

    }

    void PxSystem::clearForces() {
        for (int i = 0; i < subjs.size(); i++) {
            std::fill(subjs[i].forces.begin(), subjs[i].forces.end(), 0);
        }
    }

    void PxSystem::constrainParticles() {

    }

    void PxSystem::computeForces() {

    }

    void PxSystem::projectParticles() {

    }

    void PxSystem::updateParticles() {

    }

    // load data for particle system
    void PxSystem::load(ConfigFile& config) {
        nsubjs = config["particles.number-of-subjects"];
        nlabels = config["particles.number-of-labels"];

        Setting& settingNumParticles = config["particles.number-of-particles"];
        for (int i = 0; i < settingNumParticles.getLength(); i++) {
            numParticles.push_back(settingNumParticles[i]);
        }
        totalParticles = std::accumulate(numParticles.begin(), numParticles.end(), 0);

        subjs.resize(nsubjs);
        Setting& subjconfig = config["particles.subjects"];
        for (int i = 0; i < nsubjs; i++) {
            Setting& subjdata = subjconfig[i];
            subjs[i].regions.resize(nlabels);
            for (int j = 0; j < nlabels; j++) {
                subjs[i].regions[j].load(subjdata["labels"][j]);
            }
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
            ok = subjs[i].load(f) && subjs[i].size() == totalParticles;
        }

        return ok;
    }

    // initialize sampler on first run
    bool PxSystem::loadSampler(ConfigFile& config) {
        if (config["particles.ignore-sampler-input"]) {
            return false;
        }

        string samplerCache = config["particles.sampler-cache"];
        if (sys.checkFile(samplerCache)) {
            if (sampler.load(samplerCache) && sampler.size() == totalParticles) {
                return true;
            }
        }
        Setting& regionList = config["particles.sampler.labels"];
        if (nlabels != regionList.getLength()) {
            cout << "sampler has different number of labels" << endl;
            exit(0);
        }

        sampler.regions.resize(nlabels);
        for (int i = 0; i < regionList.getLength(); i++) {
            // if there exists intersection cache
            if (!sampler.regions[i].load(regionList[i], false)) {
                // allocate new image buffer
                LabelImage::Pointer intersection = labelIO.NewImage(subjs[0].regions[i].labelmap);

                int nPixels = intersection->GetPixelContainer()->Size();
                LabelImage::PixelType* outputBuff = intersection->GetBufferPointer();

                // initialize label pointers for fast access
                std::vector<LabelImage::PixelType*> ptrs;
                ptrs.resize(nsubjs);

                for (int j = 0; j < nsubjs; j++) {
                    ptrs[j] = subjs[j].regions[i].labelmap->GetBufferPointer();
                }

                // loop over pixels to check if the pixel is intersection
                for (int i = 0; i < nPixels; i++) {
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
        std::vector<int> numParticles;
        Setting& numberOfParticleSetting = config["particles.number-of-particles"];
        for (int i = 0; i < numberOfParticleSetting.getLength(); i++) {
            numParticles.push_back(numberOfParticleSetting[i]);
        }
        sampler.sampleParticles(numParticles);
        cout << sampler.particles << endl;

        return false;
    }

    void PxSystem::sampleParticles() {

    }


    // iterate subjects to save particle data
    bool PxSystem::saveParticles(ConfigFile& config, string outputName) {
        Setting& subjconfig = config["particles.subjects"];

        bool ok = true;
        for (int i = 0; i < subjconfig.getLength() && ok; i++) {
            string f = subjconfig[i][outputName];
            ofstream of(f);
            subjs[i].save(of);
            of.close();
        }
        return ok;
    }

    // duplicate particles from sampler to each subject
    void PxSystem::duplicateParticles() {
        for (int i = 0; i < nsubjs; i++) {
            subjs[i].resize(sampler.size());
            subjs[i].attrs = sampler.attrs;
            subjs[i].particles = sampler.particles;
        }
    }

    void PxSystem::print() {
        for (int i = 0; i < sampler.particles.size(); i++) {
            cout << sampler.particles[i] << endl;
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
        load(_config);
        if (!loadParticles(_config)) {
            if (!loadSampler(_config)) {
                initialLoop();
                string samplerCache = _config["particles.sampler-cache"];
                cout << samplerCache << endl;
                ofstream of(samplerCache);
                sampler.save(of);
                of.close();
            }
            duplicateParticles();
            saveParticles(_config, "particle-input");
        }
    }

    // initial loop to distribute particles inside region
    void PxSystem::initialLoop() {
        double t0 = _config["particles.sampler-time-steps.[0]"];
        double dt = _config["particles.sampler-time-steps.[1]"];
        double t1 = _config["particles.sampler-time-steps.[2]"];
        double sigma = _config["particles.repulsion-parameters.[1]"];
        double cutoff = _config["particles.repulsion-parameters.[2]"];

        ImageIO<LabelImage3> io;
        LabelImage::Pointer refImage = sampler.regions[0].labelmap;
        LabelImage3::Pointer tracker = createImage3(refImage, (t1 - t0) / dt + 1);

        int m = 0;
        for (double t = t0; t <= t1; t += dt, m++) {
            cout << "t = " << t << endl;
            sampler.clearForce();
            sampler.constrainParticles();

            // how to set repulsion paramter per region?
            sampler.computeRepulsion(1, sigma, cutoff);
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


        ImageIO<LabelImage3> io;
        LabelImage::Pointer refImage = subjs[0].regions[0].labelmap;
        std::vector<LabelImage3::Pointer> trackers;

        for (int i = 0; i < nsubjs; i++) {
            LabelImage3::Pointer tracker = createImage3(refImage, (t1 - t0) / dt + 1);
            trackers.push_back(tracker);
        }

        int m = 0;
        for (double t = t0; t <= t1; t += dt, m++) {
            cout << "t = " << t << endl;

            // compute internal forces
            for (int i = 0; i < nsubjs; i++) {
                subjs[i].clearForce();
                subjs[i].constrainParticles();

                // how to set repulsion paramter per region?
                subjs[i].computeRepulsion(1, 0.3, 0.5);
                subjs[i].constrainForces();
                subjs[i].updateSystem(dt);
            }

            // mark each particle location at a tracker image
            for (int i = 0; i < nsubjs; i++) {
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

        for (int i = 0; i < nsubjs; i++) {
            string file = _config["particles.particle-images"][i];
            io.WriteImage(file, trackers[i]);
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
    }

    void ParticleRunner::main(pi::Options &opts, StringVector &args) {
        try {
            _system.main(opts, args);
        } catch (ParseException & ex) {
            cout << "File: " << ex.getFile() << endl;
            cout << "Line: " << ex.getLine() << endl;
            cout << "Error: " << ex.getError() << endl;
        }
    }

    void ParticleRunner::print() {
    }


#pragma mark DemonsRunner
    void executeDemonsRunner(pi::Options &opts, StringVector &args) {
        DemonsRunner runner(opts, args);
        if (args.size() > 0) {
            if (args[0] == "optical-flow-mapping") {
                runner.computeOpticalFlowMapping();
            } else if (args[0] == "patch-mapping") {
                runner.computePatchMapping();
            }
        }
    }



    DemonsRunner::DemonsRunner(pi::Options &opts, StringVector &args) {
        _config.load("config.txt");
        if (_config.exists("build-patch")) {
            buildPatches(_config["build-patch"]);
        }
        if (_config.exists("patch-images")) {
            _config.readImages<PatchImage>(__patchImages, "patch-images");
        }
    }

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

    void DemonsRunner::computeOpticalFlowMapping() {
        Setting& files = _config["optical-flow-mapping/files"];
        Setting& param = _config["optical-flow-mapping/params"];

        double dt = 0.1;
        param.lookupValue("timestep", dt);

        int iter = 1;
        param.lookupValue("iteration", iter);

        double fieldSigma = -1;
        param.lookupValue("field-sigma", fieldSigma);

        double velocitySigma = -1;
        param.lookupValue("velocity-sigma", velocitySigma);

        bool useDemonsFlow = false;
        param.lookupValue("demons-flow", useDemonsFlow);

        bool iterativeResampling = false;
        param.lookupValue("iterative-resampling", iterativeResampling);

        cout << "demons = " << useDemonsFlow << endl;
        cout << "iterative-resampling = " << iterativeResampling << endl;
        cout << "dt = " << dt << endl;
        cout << "iter  = " << iter << endl;
        cout << "field-sigma = " << fieldSigma << endl;
        cout << "velocity-sigma = " << velocitySigma << endl;

        PatchCompare patchTool;
        for (int i = 0; i < files.getLength(); i++) {
            string source = files[i][0];
            string target = files[i][1];
            string flowoutput = files[i][2];
            string warpoutput = files[i][3];

            cout << source << " => " << target << "; " << flowoutput << ", " << warpoutput << endl;

            RealImage::Pointer sourceImage = io.ReadCastedImage(source);
            RealImage::Pointer targetImage = io.ReadCastedImage(target);
            RealImage::Pointer initialImage = io.CopyImage(sourceImage);

            const int nPixels = sourceImage->GetPixelContainer()->Size();

            DisplacementFieldType::Pointer velocityImage;
            DisplacementFieldType::Pointer flowImage;

            ImageIO<DisplacementFieldType> flowIO;
            flowImage = flowIO.NewImageS<RealImage>(sourceImage);
            DisplacementFieldType::PixelType zeroFlow;
            zeroFlow.Fill(0);
            flowImage->FillBuffer(zeroFlow);

            for (int j = 0; j < iter; j++) {
                if (useDemonsFlow) {
                    velocityImage = patchTool.computeDemonsFlow(targetImage, sourceImage, dt);
                } else {
                    velocityImage = patchTool.computeOpticalFlow(targetImage, sourceImage, dt);
                }

                if (velocitySigma > 0) {
                    cout << "apply gaussian on velocity iter: " << j << endl;
                    velocityImage = applyGaussian<DisplacementFieldType>(velocityImage, velocitySigma);
                }

                if (iterativeResampling) {
                    sourceImage = deformImage(sourceImage, velocityImage, targetImage);
                } else {
                    // compute approximated update field
                    DisplacementFieldType::Pointer updateField = resampleField(flowImage, velocityImage);
                    DisplacementFieldType::PixelType* pFlow = flowImage->GetBufferPointer();
                    DisplacementFieldType::PixelType* pUpdateField = updateField->GetBufferPointer();

                    for (int k = 0; k < nPixels; k++) {
                        fordim (l) {
                            pFlow[k][l] += pUpdateField[k][l];
                        }
                    }
                    cout << endl;

                    if (fieldSigma > 0) {
                        cout << "apply gaussian on displacement iter: " << j << endl;
                        flowImage = applyGaussian<DisplacementFieldType>(flowImage, fieldSigma);
                    }
                    sourceImage = deformImage(initialImage, flowImage, targetImage);
                }

//                io.WriteImageS<DisplacementFieldType>("velocity.nrrd", velocityImage);
//                io.WriteImageS<DisplacementFieldType>("update.nrrd", updateField);

            }


            io.WriteImage(warpoutput, sourceImage);
            io.WriteImageS<DisplacementFieldType>(flowoutput, flowImage);
        }
    }

    RealImage::Pointer DemonsRunner::deformImage(RealImage::Pointer input, DisplacementFieldType::Pointer displacement, RealImage::Pointer refImage) {
        FieldTransformType::Pointer transform = FieldTransformType::New();
        transform->SetDisplacementField(displacement);

        ResampleImageFilterType::Pointer resampler = ResampleImageFilterType::New();
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
    



    void DemonsRunner::buildPatches(libconfig::Setting &setting) {
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

}