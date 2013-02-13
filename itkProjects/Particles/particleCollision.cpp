//
//  particleTest.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/11/13.
//
//

#include "particleTest.h"
#include "piOptions.h"
#include "piImageDef.h"
#include "itkImageIO.h"
#include "piParticleCollision.h"
#include "piImageProcessing.h"

using namespace pi;

int main(int argc, char* argv[]) {
    pi::Options argOpts;
    CSimpleOpt::SOption optSpecs[] = {
        { 0, "-h", SO_NONE },
        SO_END_OF_OPTIONS
    };
    argOpts.ParseOptions(argc, argv, optSpecs);
    pi::StringVector& args = argOpts.GetStringVector("args");

    pi::ImageProcessing imageProc;

    itkcmds::itkImageIO<LabelImage> io;

    int size[] = { 80, 80, 80 };
    double center[] = { 40, 40, 40 };
    double radius[] = { 20, 20, 20 };
    LabelImage::Pointer labelImage = imageProc.Ellipse(size, center, radius);

    io.WriteImageT("/NIRAL/work/joohwi/data/ellipse/circle3.nrrd", labelImage);

    pi::ParticleCollision dyn;
    dyn.SetBinaryMask(labelImage);
    dyn.UpdateImages();
    if (args.size() > 2) {
        dyn.Write(args[1].c_str(), args[2].c_str());
    }

    int count = 0;
    LabelImage::Pointer newImage = io.NewImageT(labelImage);
    for (double theta = 0; theta <= M_2_PI; theta += (M_2_PI/20.0)) {
        for (double rho = 0; rho <= M_PI; rho += (M_PI/10.0)) {
            const double r = 35;
            double x0[] = { 40, 40, 40 };
            double x1[] = { r*cos(theta)*sin(rho)+40, r*sin(theta)*sin(rho)+40, r*cos(rho)+40 };
            VNLVector normal(__Dim);
            VNLVector regularNormal(__Dim);
            fordim (k) {
                regularNormal[k] = x1[k] - x0[k];
            }
            regularNormal.normalize();
            ContactPoint cp = { { 0, } , 0 };
//            cout << "x0: " << x2string(x0) << endl;
//            cout << "x1: " << x2string(x1) << endl;
            if (dyn.ComputeContactPoint(x0, x1, cp)) {
                LabelImage::IndexType cpIdx;
                fordim (k) {
                    cpIdx[k] = cp.cp[k]+0.5;
                }
                newImage->SetPixel(cpIdx, ++count);
                dyn.ComputeNormal(cp.cp, normal.data_block());
                normal.normalize();
                double ang = acos(dot_product(regularNormal, normal)) * 180 / M_PI;
                cout << "cp: " << x2string(cp.cp) << ";" << x2string(normal) << ";" << x2string(regularNormal) << ";" << "angle = " << ang << ";" << endl;
            }
//            dyn.Write(args[1].c_str(), args[2].c_str());
        }
    }
    io.WriteImageT("/tmpfs/cp.nrrd", newImage);
    return 0;
}