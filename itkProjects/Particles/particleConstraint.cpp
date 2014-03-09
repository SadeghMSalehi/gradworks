//
//  particleConstraint.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/12/13.
//
//

#include "particleConstraint.h"
#include "piParticleCore.h"
#include "piOptions.h"
#include "piImageProcessing.h"
#include "piParticleCollision.h"
#include "piParticleForces.h"
#include "piVectorGrid.h"

#include "itkImageIO.h"


using namespace pi;

int main(int argc, char* argv[]) {

    typedef std::vector<Particle> ParticleVector;

    ParticleVector particleTrace;

    Options parser;
//    StringVector& args = parser.ParseOptions(argc, argv, NULL);

    itkcmds::itkImageIO<LabelImage> io;
//    LabelImage::Pointer field = io.ReadImageT("/NIRAL/work/joohwi/data/ellipse/ellipse3.nrrd");
    LabelImage::Pointer field = io.ReadImageT("/tmpfs/circle_ellipse.nrrd");
    ParticleCollision collision;
    collision.SetBinaryMask(field);
    collision.UseDistanceMapCache("/tmpfs/circle_ellipse_distancemap.nrrd");
    collision.UpdateImages();

    itkcmds::itkImageIO<DoubleImage> iox;
    DoubleImage::Pointer outField = iox.NewImageT<LabelImage>(field);


    double t0 = 0;
    double t1 = 100;
    double dt = 0.05;

    VectorGrid grid;
    grid.CreateScalar("time");
    grid.CreateVector("normal");
    grid.CreateVector("velocity");
    grid.CreateVector("tangent");
    
    LabelImage::IndexType pIdx;
    LabelImage::RegionType fieldRegion = field->GetBufferedRegion();

    VNLVector normal(__Dim, 0);
    VNLVector force(__Dim, 0);
    VNLVector tangent(__Dim, 0);

    const int Nx = 8, Ny = 8, Nz = 8;
    const double S = 0.5;

    ParticleSubjectArray subs;
    subs.resize(1);
    subs[0].NewParticles(Nx*Ny*Nz);

    const bool useRandomSampling = true;
    if (useRandomSampling) {
        subs[0].InitializeRandomPoints(field);
    } else {
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j ++) {
                for (int k = 0; k < Nz; k++) {
                    int m = j%2;
                    int n = 0;
                    arrayset3(subs[0][Ny*Nz*i+Nz*j+k].x, 40+(i-Nx/2.0+m/2.0)*S, 40+(j-Ny/2.0+n/2.0)*S, 40);
                    arrayset3(subs[0][Ny*Nz*i+Nz*j+k].v, 0, 0, 0);
                }
            }
        }
    }
    

    EntropyInternalForce internalForce;

    const int nPoints = subs[0].GetNumberOfPoints();
    for (double t = t0; t <= t1; t += dt) {
        cout << "t: " << t << endl;

        for (int n = 0; n < 1; n++) {
            for (int i = 0; i < subs[n].GetNumberOfPoints(); i++) {
                // system setup
                particleTrace.push_back(subs[n][i]);
                forset(subs[n][i].x, subs[n][i].w);
                forfill(subs[n][i].f, 0);
            }
        }

        internalForce.ComputeForce(subs);

        // collision handling
        collision.HandleCollision(subs);

        // system update
        for (int n = 0; n < 1; n++) {
            for (int i = 0; i < nPoints; i++) {
                fordim(k) {
                    subs[n][i].f[k] -= 1 * subs[n][i].v[k];
                    subs[n][i].v[k] += dt * subs[n][i].f[k];
                    subs[n][i].x[k] += dt * subs[n][i].v[k];
                }

                Particle &p = subs[n][i];
//                cout << t << ": x = (" << x2string(p.x) << "); v = (" << x2string(p.v) << "); f = (" << x2string(p.f) << "); " << endl;
                if (!fieldRegion.IsInside(pIdx)) {
                    cout << "Stop system: out of region" << endl;
                    goto quit;
                }
                fordim(k) {
                    pIdx[k] = p.x[k] + 0.5;
                }
                outField->SetPixel(pIdx, t);
                if (dimequal(p.v,0,0,0)) {
                    cout << "Stop system: " << t << endl;
                    goto quit;
                }
                subs[n][i].t = t;
            }
        }
    }

quit:
    for (int n = 0; n < 1; n++) {
        for (int i = 0; i < subs[n].GetNumberOfPoints(); i++) {
            // system setup
            particleTrace.push_back(subs[n][i]);
        }
    }

    ofstream of("/tmpfs/trace.txt");
    for (int i = 0; i < particleTrace.size(); i++) {
        of << particleTrace[i].t << " " << particleTrace[i].idx << " " << particleTrace[i] << endl;
    }
    of.close();
    const char* particleOutput = "/tmpfs/particle.nrrd";
    iox.WriteImageT(particleOutput, outField);
}