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
#include "itkImageIO.h"

using namespace pi;

int main(int argc, char* argv[]) {
    Options parser;
    StringVector& args = parser.ParseOptions(argc, argv, NULL);
    
    itkcmds::itkImageIO<LabelImage> io;
    LabelImage::Pointer field = io.ReadImageT("/NIRAL/work/joohwi/data/ellipse/circle3.nrrd");
    LabelImage::Pointer outField = io.NewImageT(field);

    ParticleCollision collision;
    collision.SetBinaryMask(field);
    collision.UpdateImages();
    

    Particle p;
    arrayset3(p.x, 50, 40, 40);
    arrayset3(p.v, 0, 0, 3);
    arrayset3(p.f, 0, 0.98, 0);

    double t0 = 0;
    double t1 = 300;
    double dt = 0.05;

    LabelImage::IndexType pIdx;
    LabelImage::RegionType fieldRegion = field->GetBufferedRegion();
    ContactPoint cp = { { 0 }, 0 };
    double center[__Dim] = { 40, 40, 40 };
    VNLVector normal(__Dim);
    VNLVector force(__Dim);
    VNLVector tangent(__Dim);
    int cnt = 0;
    for (double t = t0; t <= t1; t += dt) {
        // system setup
        forset(p.x, p.w);
        fordim(k) {
            force[k] = 0;//center[k] - p.x[k];
        }
        force.normalize();
        fordim(k) {
            p.f[k] = -force[k];
            p.v[k] += dt * p.f[k];
            p.x[k] += dt * p.v[k];
        }

        LabelImage::IndexType idx;
        fordim (k) {
            idx[k] = p.x[k];
        }
        if (!collision.IsRegionInside(idx)) {
            collision.ComputeClosestBoundary(p.x, p.x);
        }

        // system validity check
        forroundset(p.x, pIdx);
        if (!fieldRegion.IsInside(pIdx)) {
            cout << "Stop system: out of region" << endl;
            break;
        }

        // collision detection
        if (collision.ComputeContactPoint(p.w, p.x, cp)) {
            // collision response
            collision.ComputeNormal(cp.cp, normal.data_block());
            normal.normalize();
            double nv = dimdot(p.v, normal);

            const bool rebound = false;
            if (rebound) {
                fordim (k) {
                    p.v[k] -= 2 * nv * normal[k];
                    p.x[k] = cp.cp[k];
                }
            } else {
                double vv = sqrt(dimnorm2(p.v));
                fordim(k) {
                    tangent[k] = p.v[k] - nv * normal[k];
                }
                fordim(k) {
                    p.v[k] = tangent[k];
                    p.x[k] = cp.cp[k];
                }
                if (cp.status == STARTING_CONTACT) {
                    fordim (k) {
                        p.x[k] += dt * p.v[k];
                    }
                }
            }
        }


        cout << t << ": " << x2string(p.x) << "; " << x2string(p.v) << "; " << endl;
        if (!fieldRegion.IsInside(pIdx)) {
            cout << "Stop system: out of region" << endl;
            break;
        }
        outField->SetPixel(pIdx, ++cnt);
        if (dimequal(p.v,0,0,0)) {
            cout << "Stop system: " << t << endl;
            break;
        }
    }

    io.WriteImageT("/tmpfs/particle.nrrd", outField);
}