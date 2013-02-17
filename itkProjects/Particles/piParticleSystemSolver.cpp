//
//  piParticleSystemSolver.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/15/13.
//
//

#include "piParticleSystemSolver.h"
#include "piParticleCollision.h"
#include "piParticleForces.h"
#include "piOptions.h"
#include "piParticleTrace.h"

#include "fstream"
#include "sstream"

namespace pi {
    using namespace std;
    
    bool ParticleSystemSolver::LoadConfig(const char* name) {
        ifstream in(name);
        if (!in.is_open()) {
            cout << "Can't open " << name << endl;
            return false;
        }
        in >> m_Options;
        
        cout << m_Options << endl;
        
        // this requires 'NumberOfParticles:', 'Subjects:'
        m_System.InitializeSystem(m_Options);

        StringVector& labelImages = m_Options.GetStringVector("LabelImages:");
        for (int i = 0; i < labelImages.size(); i++) {
            m_ImageContext.LoadLabel(labelImages[i]);
        }
        
        StringVector& realImages = m_Options.GetStringVector("RealImages:");
        for (int i = 0; i < realImages.size(); i++) {
            m_ImageContext.LoadDoubleImage(realImages[i]);
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
                        m_System[subjId].ReadParticles(in, nPoints);
                    }
                } else if (name == "InitialParticles:") {
                    int nPoints;
                    ss >> nPoints;
                    m_System.GetInitialSubject().ReadParticlePositions(in, nPoints);
                    cout << "Read " << nPoints << " initial particles" << endl;
                }
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
        const int nInitialPoints = m_System.GetInitialSubject().GetNumberOfPoints();
        if (nInitialPoints > 0) {
            out << "InitialParticles: " << nInitialPoints << endl;
            m_System.GetInitialSubject().WriteParticlePositions(out);
        }
        const int nSubj = m_System.GetNumberOfSubjects();
        const int nParticles = m_System.GetNumberOfParticles();
        if (nSubj > 0) {
            for (int i = 0; i < nSubj; i++) {
                out << "Particles: " << i << " " << nParticles << endl;
                m_System[i].WriteParticles(out);
            }
        }
        return true;
    }
    
    void ParticleSystemSolver::Preprocessing() {
        string traceFile = m_Options.GetString("PreprocessingTrace:", "");
        const bool traceOn = traceFile != "";
        ParticleTrace trace;
        
        ParticleSubject& initial = m_System.GetInitialSubject();
        if (initial.GetNumberOfPoints() == 0) {
            itkcmds::itkImageIO<LabelImage> io;
            string intersectionCache = m_Options.GetString("InitialIntersectionMaskCache:", "");
            if (intersectionCache == "" || !io.FileExists(intersectionCache.c_str())) {
                m_ImageContext.ComputeIntersection();
                if (intersectionCache != "") {
                    io.WriteImageT(intersectionCache.c_str(), m_ImageContext.GetIntersection());
                }
            } else {
                m_ImageContext.SetIntersection(io.ReadImageT(intersectionCache.c_str()));
            }
            
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

        double t0 = m_Options.GetDoubleVectorValue("PreprocessingTimeRange:", 0);
        double dt = m_Options.GetDoubleVectorValue("PreprocessingTimeRange:", 1);
        double t1 = m_Options.GetDoubleVectorValue("PreprocessingTimeRange:", 2);
        
        
        ParticleCollision boundary;
        boundary.SetBinaryMask(m_ImageContext.GetIntersection());
        boundary.UseBinaryMaskSmoothing();
        boundary.UseBinaryMaskSmoothingCache(m_Options.GetString("InitialIntersectionMaskCache:", "").c_str());
        boundary.UseDistanceMapCache(m_Options.GetString("InitialIntersectionDistanceMapCache:", "").c_str());
        boundary.UpdateImages();
        
        //            SaveConfig("/tmpfs/temp.txt");
        
        EntropyInternalForce internalForce;
        
        const int nPoints = initial.GetNumberOfPoints();
        for (double t = t0; t < t1; t += dt) {
            cout << "t: " << t << endl;
            ParticleSubject& sub = initial;

            for (int i = 0; i < nPoints; i++) {
                Particle& pi = sub[i];
                forset(pi.x, pi.w);
                forfill(pi.f, 0);
            }
            
            internalForce.ComputeForce(initial);
            boundary.HandleCollision(initial);
            
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

    void Run2() {
        /*
         const double t0 = m_times[0];
         const double dt = m_times[1];
         const double t1 = m_times[2];
         
         
         cout << "Distance map generation ..." << endl;
         ParticleConstraint constraint;
         constraint.SetImageList(m_ImageContext.GetLabelVector());
         cout << "Distance map generation ... done" << endl;
         
         
         char trackName[128];
         int k = 0;
         boost::timer timer;
         for (double t = t0; t < t1; t += dt) {
         timer.restart();
         ComputeMeanSubject();
         if (m_InternalForceFlag) {
         InternalForce internalForce;
         internalForce.ComputeForce(m_Subjects);
         }
         
         EnsembleForce ensembleForce(m_EnsembleCoeff);
         ensembleForce.SetImageContext(&m_ImageContext);
         if (m_EnsembleForceFlag) {
         ensembleForce.ComputeEnsembleForce(m_Subjects);
         }
         
         IntensityForce IntensityForce(m_IntensityCoeff);
         IntensityForce.SetImageContext(&m_ImageContext);
         if (m_IntensityForceFlag) {
         IntensityForce.ComputeIntensityForce(this);
         }
         
         constraint.ApplyConstraint(m_Subjects);
         UpdateSystem(m_Subjects, dt);
         if (m_TrackingOutputPattern != "") {
         sprintf(trackName, m_TrackingOutputPattern.c_str(), ++k);
         SaveSystem(trackName);
         }
         cout << "Processing Time: " << t << "; Elapsed " << timer.elapsed() << " secs" << endl;
         }
         */
        
    }
    
    void ParticleSystemSolver::Run() {
        ParticleSystem& system = m_System;
        const int nSubz = system.GetNumberOfSubjects();
        const int nPoints = system.GetNumberOfParticles();
        ParticleSubjectArray& subs = system.GetSubjects();

        string traceFile = m_Options.GetString("RunTrace:", "");
        const bool traceOn = traceFile != "";

        ParticleTrace trace;
        trace.Resize(nSubz);

        for (int n = 0; n < nSubz; n++) {
            m_System[n].Initialize(system.GetInitialSubject().m_Particles);
        }
        
        EntropyInternalForce internalForce;

        EnsembleForce ensembleForce(1);
        ensembleForce.SetImageContext(&m_ImageContext);

        IntensityForce intensityForce(1);
        intensityForce.SetImageContext(&m_ImageContext);

        std::vector<ParticleCollision> collisionHandlers;
        collisionHandlers.resize(nSubz);
        
        for (int n = 0; n < nSubz; n++) {
            collisionHandlers[n].SetBinaryMask(m_ImageContext.GetLabel(n));
            collisionHandlers[n].UseBinaryMaskSmoothing();
            collisionHandlers[n].UseBinaryMaskSmoothingCache(m_Options.GetStringVectorValue("BinaryMaskSmoothingCache:", n).c_str());
            collisionHandlers[n].UseDistanceMapCache(m_Options.GetStringVectorValue("BinaryMaskDistanceMapCache:", n).c_str());
            collisionHandlers[n].UpdateImages();
        }
        
        const double t0 = system.GetSystemOptions().GetDoubleVectorValue("TimeRange:", 0);
        const double dt = system.GetSystemOptions().GetDoubleVectorValue("TimeRange:", 1);
        const double t1 = system.GetSystemOptions().GetDoubleVectorValue("TimeRange:", 2);
        const bool useEnsemble = m_Options.GetBool("+ensemble", false);
        const bool useIntensity = m_Options.GetBool("+intensity", false);
        const bool useNoInternal = m_Options.GetBool("-internal", false);
        
        for (double t = t0; t < t1; t += dt) {
            cout << "t: " << t << endl;
            m_System.ComputeMeanSubject();

            // compute internal force
            if (!useNoInternal) {
                for (int n = 0; n < nSubz; n++) {
                    ParticleSubject& sub = subs[n];
                    for (int i = 0; i < nPoints; i++) {
                        Particle& pi = sub[i];
                        forset(pi.x, pi.w);
                        forfill(pi.f, 0);
                    }
                    
                    internalForce.ComputeForce(sub);
                    collisionHandlers[n].HandleCollision(sub);
                }
            }

            if (useEnsemble) {
                ensembleForce.ComputeEnsembleForce(m_System);
            }

            if (useIntensity) {
                intensityForce.ComputeIntensityForce(&m_System);
            }
            
            // system update
            for (int n = 0; n < nSubz; n++) {
                ParticleSubject& sub = subs[n];
                for (int i = 0; i < nPoints; i++) {
                    Particle& p = sub[i];
                    LabelImage::IndexType pIdx;
                    fordim (k) {
                        p.f[k] -= p.v[k];
                        p.v[k] += dt*p.f[k];
                        p.x[k] += dt*p.v[k];
                        pIdx[k] = p.x[k] + 0.5;
                    }

                    if (!collisionHandlers[n].IsBufferInside(pIdx)) {
                        cout << "Stop system: out of region" << endl;
                        goto quit;
                    }
                }
                if (traceOn) {
                    trace.Add(t, sub);
                }
            }
        }

    quit:
        if (traceOn) {
            ofstream traceOut(traceFile.c_str());
            trace.Write(traceOut);
            traceOut.close();
        }

        return;
    }
}