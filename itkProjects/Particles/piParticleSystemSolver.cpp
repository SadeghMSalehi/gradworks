//
//  piParticleSystemSolver.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/15/13.
//
//

#include <fstream>
#include <sstream>
#include <itkMacro.h>

#include "piParticleSystemSolver.h"
#include "piParticleCollision.h"
#include "piParticleForces.h"
#include "piOptions.h"
#include "piParticleTrace.h"
#include "piImageProcessing.h"
#include "piTimer.h"
#include "piImageIO.h"
#include "piParticleBSpline.h"

namespace pi {
    using namespace std;

    ParticleSystemSolver::ParticleSystemSolver() : verbose(true) {
        useAlignment = false;
    }

    bool ParticleSystemSolver::LoadConfig(const char* name) {
        ifstream in(name);
        if (!in.is_open()) {
            cout << "Can't open " << name << endl;
            return false;
        }
        cout << "loading " << name << " ..." << endl;
        return LoadConfig(in);
    }

    bool ParticleSystemSolver::LoadConfig(std::istream &in) {
        m_Options = Options();

        in >> m_Options;

        // debug
        cout << m_Options << endl;

        // this requires 'NumberOfParticles:', 'Subjects:'
        m_System = ParticleSystem();
        m_System.InitializeSystem(m_Options);
        
        m_ImageContext.Clear();
        StringVector& labelImages = m_Options.GetStringVector("LabelImages:");
        for (int i = 0; i < labelImages.size(); i++) {
            m_ImageContext.LoadLabel(labelImages[i]);
        }
        
        StringVector& realImages = m_Options.GetStringVector("RealImages:");
        for (int i = 0; i < realImages.size(); i++) {
            m_ImageContext.LoadRealImage(realImages[i]);
        }
        
        while (in.good()) {
            char buf[1024];
            in.getline(buf, sizeof(buf));
            if (in.good()) {
                stringstream ss(buf);
                string name;
                ss >> name;
                
                if (name == "Particles:") {
                    int subjId, nPoints;
                    ss >> subjId >> nPoints;
                    if (subjId < m_System.GetNumberOfSubjects()) {
                        if (m_Options.GetBool("load_position_only")) {
                            m_System[subjId].ReadParticlePositions(in, nPoints);
                        } else {
                            m_System[subjId].ReadParticles(in, nPoints);
                        }
                    }
                } else if (name == "InitialParticles:") {
                    int nPoints;
                    ss >> nPoints;
                    m_System.GetInitialSubject().ReadParticlePositions(in, nPoints);
                    cout << "Read " << nPoints << " initial particles" << endl;
                }
//                else if (name == "ParticleAlignment:") {
//                    int nsubj;
//                    ss >> nsubj;
//                    for (int i = 0; i < nsubj; i++) {
//                        in.getline(buf, sizeof(buf));
//                        string line(buf);
//                        stringstream ss(line);
//                        ss >> m_System[i].alignment;
//                    }
//                }
            }
        }
        return true;
    }
    
    bool ParticleSystemSolver::SaveConfig(const char* name) {
        ofstream out(name);
        if (!out.is_open()) {
            cout << "Can't write " << name << endl;
            return false;
        }
        out << m_Options << OPTION_END << endl << endl;
        const int nSubj = m_System.GetNumberOfSubjects();
        const int nParticles = m_System.GetNumberOfParticles();

//        if (nSubj > 0) {
//            out << "ParticleAlignment: " << nSubj << " " << nParticles << endl;
//        for (int i = 0; i < nSubj; i++) {
//            out << m_System[i].alignment << endl;
//        }
//        }
        const int nInitialPoints = m_System.GetInitialSubject().GetNumberOfPoints();
        if (nInitialPoints > 0) {
            out << "InitialParticles: " << nInitialPoints << endl;
            m_System.GetInitialSubject().WriteParticlePositions(out);
        }

        if (nSubj > 0) {

            for (int i = 0; i < nSubj; i++) {
                out << "Particles: " << i << " " << nParticles << endl;
                if (m_Options.GetBool("save_position_only")) {
                    m_System[i].WriteParticlePositions(out);
                } else {
                    m_System[i].WriteParticles(out);
                }
            }
        }
        return true;
    }
    
    void ParticleSystemSolver::Preprocessing() {
        string traceFile = m_Options.GetString("PreprocessingTrace:", "");
        const bool traceOn = traceFile != "";
        ParticleTrace trace;
        
        ParticleCollision boundary;
        boundary.applyMaskSmoothing = true;
        m_Options.GetStringTo("InitialIntersectionMaskCache:", boundary.binaryMaskCache);
        m_Options.GetStringTo("InitialIntersectionDistanceMapCache:", boundary.distanceMapCache);
        
        ParticleSubject& initial = m_System.GetInitialSubject();
        if (initial.GetNumberOfPoints() == 0) {
            
            // compute intersection in combination with particle boundary
            ImageIO<LabelImage> io;
            if (io.FileExists(boundary.binaryMaskCache.c_str())) {
                boundary.SetBinaryMask(io.ReadCastedImage(boundary.binaryMaskCache.c_str()));
            } else {
                int npixels = m_ImageContext.ComputeIntersection();
                if (npixels == 0) {
                    cout << "invalid intersection.. halt" << endl;
                    exit(0);
                    return;
                }
                boundary.SetLabelImage(m_ImageContext.GetIntersection());
            }
            
            boundary.UpdateImages();
            m_ImageContext.SetIntersection(boundary.GetBinaryMask());
            
            initial.m_Name = "Intersection";
            initial.NewParticles(m_Options.GetInt("NumberOfParticles:", 0));
            initial.InitializeRandomPoints(m_ImageContext.GetIntersection());
            string initialConfigOutput = m_Options.GetString("InitialConfigOutput:", "");
            if (initialConfigOutput != "") {
                SaveConfig(initialConfigOutput.c_str());
            }
        } else {
            cout << "Skip preprocessing and use pre-generated data..." << endl;
            return;
        }
        
        if (initial.GetNumberOfPoints() == 0) {
            cout << "Fail to initializing with random points" << endl;
            return;
        }
        
        DataReal t0 = m_Options.GetRealVectorValue("PreprocessingTimeRange:", 0);
        DataReal dt = m_Options.GetRealVectorValue("PreprocessingTimeRange:", 1);
        DataReal t1 = m_Options.GetRealVectorValue("PreprocessingTimeRange:", 2);
        
        EntropyInternalForce internalForce;
        m_Options.GetRealTo("InternalForceSigma:", internalForce.repulsionSigma);
        m_Options.GetRealTo("InternalForceCutoff:", internalForce.repulsionCutoff);
        
        // this must be off for preprocessing
        internalForce.useAdaptiveSampling = false;
        
        const int nPoints = initial.GetNumberOfPoints();

        int count = 0;
        Timer timer;
        timer.start();
        for (DataReal t = t0; t < t1; t += dt, count++) {
            ParticleSubject& sub = initial;
            
            for (int i = 0; i < nPoints; i++) {
                Particle& pi = sub[i];
                forset(pi.x, pi.w);
                forfill(pi.f, 0);
            }
            boundary.ConstrainPoint(sub);
            internalForce.ComputeForce(initial);
            boundary.ProjectForceAndVelocity(sub);
            
            for (int i = 0; i < nPoints; i++) {
                Particle& p = sub[i];
                LabelImage::IndexType pIdx;
                fordim (k) {
                    p.f[k] -= p.v[k];
                    p.v[k] += dt*p.f[k];
                    p.x[k] += dt*p.v[k];
                    pIdx[k] = p.x[k] + 0.5;
                }
                
                
                if (!boundary.IsBufferInside(pIdx)) {
                    cout << "Stop system: out of region" << endl;
                    goto quit;
                }
            }
            
            if (traceOn) {
                trace.Add(t, sub.m_Particles, 0);
            }
            
            double elapsedTime = timer.getElapsedTimeInSec();

            if (count % 1000 == 0) {
                cout << "t: " << t << " " << flush;
                cout << "; elapsed time: " << elapsedTime << " sec" << endl;
                timer.start();
            }
        }
        
        if (traceOn) {
            ofstream out(traceFile.c_str());
            trace.Write(out);
            out.close();
        }
    quit:
        cout << "Preprocessing done ..." << endl;
        return;
    }
    

    void ParticleSystemSolver::Setup() {
        ParticleSystem& system = m_System;
        int nSubz = system.size();

        m_Options.GetStringTo("RunTrace:", traceFile);
        traceOn = traceFile != "";

        m_Options.GetStringTo("SystemSnapshot:", systemSnapshot);

        trace.Clear();
        trace.Resize(nSubz);

        if (!m_Options.GetBool("use_previous_position")) {
            if (system.GetInitialSubject().size() > 0) {
                ImageProcessing proc;
                for (int n = 0; n < nSubz; n++) {
                    ParticleSubject& s = m_System[n];
                    s.Initialize(system.GetInitialSubject().m_Particles);
                    s.friendImage = m_ImageContext.GetLabel(n);
                    s.friendSampler = NNLabelInterpolatorType::New();
                    s.friendSampler->SetInputImage(s.friendImage);
                    s.realImage = m_ImageContext.GetRealImage(n);
                    s.realSampler = LinearImageInterpolatorType::New();
                    s.realSampler->SetInputImage(s.realImage);
                    s.gradImage = proc.ComputeGaussianGradient(m_System[n].realImage, 1);
                    s.gradSampler = GradientInterpolatorType::New();
                    s.gradSampler->SetInputImage(s.gradImage);
                }
            }
        }
        m_System.InitializeMean();


        useEnsemble = m_Options.GetBool("ensemble");
        useIntensity = m_Options.GetBool("intensity");
        noInternal = m_Options.GetBool("no_internal");
        noBoundary = m_Options.GetBool("no_boundary");
        useVelocity = m_Options.GetBool("use_velocity");
        useAlignment = m_Options.GetBool("use_alignment");

        m_Options.SetBool("use_velocity", useVelocity);
        m_Options.SetBool("intensity", useIntensity);
        m_Options.SetBool("no_internal", noInternal);
        m_Options.SetBool("no_boundary", noBoundary);
        m_Options.SetBool("ensemble", useEnsemble);
        m_Options.SetBool("use_alignment", useAlignment);
        
        // check label images are loaded
        if (!noBoundary && nSubz != m_ImageContext.GetLabelVector().size()) {
            cout << "the same number of boundary mask files are required" << endl;
            return;
        }

        m_Options.GetRealTo("InternalForceSigma:", internalForce.repulsionSigma);
        m_Options.GetRealTo("InternalForceCutoff:", internalForce.repulsionCutoff);
        m_Options.GetRealTo("InternalForceFriendSigma:", internalForce.friendSigma);
        m_Options.GetRealTo("InternalForceFriendCutoff:", internalForce.friendCutoff);
        m_Options.GetBoolTo("adaptive_sampling", internalForce.useAdaptiveSampling);
        m_Options.GetBoolTo("multiphase_force", internalForce.useMultiPhaseForce);


        if (internalForce.useAdaptiveSampling) {
            m_System.LoadKappaImages(m_Options, &m_ImageContext);
        }
        if (internalForce.repulsionCutoff < internalForce.repulsionSigma) {
            internalForce.repulsionCutoff = internalForce.repulsionSigma * 5;
            cout << "repulsion sigma: " << internalForce.repulsionSigma << endl;
            cout << "adjusted repulsion cutoff: " << internalForce.repulsionCutoff << endl;
        }



        ensembleForce.SetImageContext(&m_ImageContext);

        intensityForce.SetImageContext(&m_ImageContext);
        intensityForce.useAttributesAtWarpedSpace = true;

        internalForce.coeff = 0.3;
        ensembleForce.coeff = 0.3;
        intensityForce.coeff = 0.3;

        m_Options.GetRealVectorValueTo("ForceCoefficients:", 0, internalForce.coeff);
        m_Options.GetRealVectorValueTo("ForceCoefficients:", 1, ensembleForce.coeff);
        m_Options.GetRealVectorValueTo("ForceCoefficients:", 2, intensityForce.coeff);

        RealVector forceCoefficients = m_Options.GetRealVector("ForceCoefficients:");
        forceCoefficients.clear();
        forceCoefficients.push_back(internalForce.coeff);
        forceCoefficients.push_back(ensembleForce.coeff);
        forceCoefficients.push_back(intensityForce.coeff);


#ifdef ATTR_SIZE
        m_Options.SetInt("AttributeDimension:", ATTR_SIZE);
#endif


        ///////////////////////////////////////////////
        // Collison Handler Setup
        //
        collisionHandlers.resize(nSubz);

        for (int n = 0; n < nSubz; n++) {
            collisionHandlers[n].subject = &m_System[n];
            collisionHandlers[n].applyMaskSmoothing = true;
            m_Options.GetStringVectorValueTo("BinaryMaskCache:", n, collisionHandlers[n].binaryMaskCache);
            m_Options.GetStringVectorValueTo("BinaryMaskDistanceMapCache:", n, collisionHandlers[n].distanceMapCache);
            collisionHandlers[n].SetLabelImage(m_ImageContext.GetLabel(n));
            collisionHandlers[n].UpdateImages();
        }


        ////////////////////////////////////////////////
        // Time Range Setup
        //
        RealVector& times = m_Options.GetRealVector("TimeRange:");
        if (times.size() > 0) {
            t0 = times[0];
            dt = times[1];
            t1 = times[2];
        }

        ////////////////////////////////////////////////
        // Active Option Listing
        //
        if (useEnsemble) cout << "ensemble term enabled" << endl; else cout << "ensemble term disabled" << endl;
        if (useIntensity) cout << "intensity term enabled" << endl; else cout << "intensity term disabled" << endl;
        if (noInternal) cout << "internal force disabled" << endl; else cout << "internal force enabled" << endl;
        if (noBoundary) cout << "boundary term disabled" << endl; else cout << "boundary term enabled" << endl;
        if (useVelocity) cout << "velocity system enabled" << endl; else cout << "velocity system disabled" << endl;
        if (internalForce.useAdaptiveSampling) cout << "adaptive sampling enabled" << endl; else cout << "adaptive sampling disabled" << endl;
        if (internalForce.useMultiPhaseForce) cout << "multi-phase force enabled" << endl; else cout << "multi-phase force disabled" << endl;

        cout << "ForceCoefficients: " << internalForce.coeff << ", " << ensembleForce.coeff << ", " << intensityForce.coeff << endl;
        
        cout << Attr::NATTRS << " attributes per particle" << endl;

        m_System.currentIteration = -1;
        //        ofstream err("error.txt");


    }
    
    void ParticleSystemSolver::Run() {
        ParticleSystem& system = m_System;
        const int nSubz = system.GetNumberOfSubjects();
        Setup();
        try {
            int stepResult = 0;
            /////////////////////////////////
            // t is automatically increases
            ////////////////////////////////
            for (t = t0; t < t1;) {
                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                // IMPORTANT STEP
                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                stepResult = RunStep();
            }
            if (stepResult == 0) {
                // last collision handling
                for (int i = 0 ; i < nSubz; i++) {
                    collisionHandlers[i].HandleCollision(system[i]);
                }
                cout << "Done ..." << endl;
            }
        } catch (itk::ExceptionObject& ex) {
            ex.Print(cout);
        }
//        err.close();
        if (traceOn) {
            ofstream traceOut(traceFile.c_str());
            trace.Write(traceOut);
            traceOut.close();
        }
        return;
    }



    //////////////////////////////////////////////////////
    //
    // Actual Computation: IMPORTANT!!
    //
    //
    int ParticleSystemSolver::RunStep() {
        ParticleSystem& system = m_System;

        const int nSubz = system.GetNumberOfSubjects();
        const int nPoints = system.GetNumberOfParticles();
        ParticleSubjectArray& subs = system.GetSubjects();

        timer.start();
        m_System.currentTime = t;
        m_System.currentIteration++;

        RunStepBegin();

        // compute internal force
        for (int n = 0; n < nSubz; n++) {
            for (int i = 0; i < nPoints; i++) {
                Particle& pi = m_System[n][i];
                forfill(pi.f, 0);
            }
            collisionHandlers[n].ConstrainPoint(m_System[n]);
            m_System[n].ComputeDensity();
        }

        // compute alignment and apply internal force
        if (!noInternal) {
            if (nSubz > 1 && useAlignment) {
                m_System.ComputeXMeanSubject();
                for (int n = 0; n < nSubz; n++) {
                    internalForce.ComputeForce(m_System[n]);
                    // alignment for later procedures
                    m_System[n].ComputeAlignment(m_System.GetMeanSubject());
                    m_System[n].AlignmentTransformX2Y();
                }
                m_System.ComputeYMeanSubject();
            } else {
                for (int n = 0; n < nSubz; n++) {
                    internalForce.ComputeForce(m_System[n]);
                    m_System[n].TransformX2Y();
                    m_System[n].TransformY2Z();
                }
            }
        } else {
            for (int n = 0; n < nSubz; n++) {
                m_System[n].TransformX2Y();
            }
        }

        // work on coordinate Y
        if (useEnsemble) {
            ensembleForce.ComputeEnsembleForce(m_System);
        }

        // work on coordinate Z
        if (useIntensity) {
            intensityForce.ComputeIntensityForce(&m_System);
        }

        // system update
        for (int n = 0; n < nSubz; n++) {
            ParticleSubject& sub = subs[n];
            collisionHandlers[n].ProjectForceAndVelocity(sub);
            for (int i = 0; i < nPoints; i++) {
                Particle& p = sub[i];
                LabelImage::IndexType pIdx;
                if (useVelocity) {
                    fordim (k) {
                        p.f[k] -= p.v[k];
                        p.v[k] += dt*p.f[k];
                        p.x[k] += dt*p.v[k];
                        pIdx[k] = p.x[k] + 0.5;
                    }
                } else {
                    fordim (k) {
                        p.x[k] += dt*p.f[k];
                        pIdx[k] = p.x[k] + 0.5;
                    }
                }
                if (!collisionHandlers[n].IsBufferInside(pIdx)) {
                    cout << "\nStop system: out of region" << endl;
                    return 1;
                }
            }
        }
        
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // IMPORTANT
        t += dt;
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        RunStepEnd();
        return 0;
    }

    void ParticleSystemSolver::RunStepBegin() {
        if (verbose) {
            cout << "t: " << t << " " << flush;
        }
        if (systemSnapshot != "") {
            SaveConfig(systemSnapshot.c_str());
        }
    }

    
    void ParticleSystemSolver::RunStepEnd() {
        ParticleSubjectArray& subs = m_System.GetSubjects();
        const int nSubz = subs.size();

        double elapsedTime = timer.getElapsedTimeInSec();
        if (verbose) {
            cout << "; elapsed time: " << elapsedTime << " secs; estimated remaining time: about " << int((elapsedTime*(t1 - t)/dt)/60+1) << " mins                   \r" << flush;
        }

        if (traceOn) {
            for (int n = 0; n < nSubz; n++) {
                ParticleSubject& sub = subs[n];
                if (std::abs(int(t + dt * 0.1)-t) < dt * 0.1) {
                    trace.Add(t, sub);
                }
            }
        }
    }

    RealImage::Pointer ParticleSystemSolver::WarpImage(int i, int j) {
        if (j < 0) {
            // warp to the mean image
            ParticleSubject& meanSubj = m_System.ComputeXMeanSubject();
//            for (int i = 0; i < meanSubj.size(); i++) {
//                fordim (k) {
//                    cout << meanSubj[i].x[k] << " ";
//                }
//                cout << endl;
//            }
            ParticleBSpline bspline;
            bspline.SetReferenceImage(m_ImageContext.GetLabel(i));
            bspline.EstimateTransform(meanSubj, m_System[i]);
            return bspline.WarpImage(m_ImageContext.GetRealImage(i));
        }
        return RealImage::Pointer();
    }

    LabelImage::Pointer ParticleSystemSolver::WarpLabel(int i, int j) {
        return LabelImage::Pointer();
    }

    void ParticleSystemSolver::PrintPoints() {
        for (int j = 0; j < m_System.GetNumberOfParticles(); j++) {
            for (int i = 0; i < m_System.GetNumberOfSubjects(); i++) {
                fordim (k) {
                    cout << m_System[i][j].x[k] << ",";
                }
                cout << "\t";
            }
            cout << endl;
        }
    }
}