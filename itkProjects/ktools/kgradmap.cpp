//
//  gradmaps.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/5/13.
//
//

#include <stdio.h>
#include "piImageDef.h"
#include "piImageProcessing.h"
#include "piOptions.h"
#include "piVTK.h"
#include "piTimer.h"
#include "piImageIO.h"
#include <vtkKdTreePointLocator.h>
#include <vtkPointData.h>

using namespace std;
using namespace pi;
using namespace pivtk;

ImageIO<RealImage> io;

int main(int argc, char* argv[]) {
    CSimpleOpt::SOption specs[] = {
        { 0, "-t", SO_REQ_SEP },
        { 1, "-s", SO_REQ_SEP },
        { 2, "--gradmag", SO_REQ_SEP },
        { 3, "--phi", SO_REQ_SEP },
        { 4, "--theta", SO_REQ_SEP },
        SO_END_OF_OPTIONS
    };

    Options parser;
    StringVector args = parser.ParseOptions(argc, argv, specs);

    cout << "option: " << parser << endl;
    double thresh = parser.GetStringAsReal("-t", 1000);
    cout << "threshold: " << thresh << endl;
    double sigma = parser.GetStringAsReal("-s", 1);
    cout << "gaussian sigma: " << sigma << endl;

    RealImage::Pointer img = io.ReadCastedImage(args[0].c_str());
    ImageProcessing proc;

    Timer t;
    t.start();
    GradientImage::Pointer gradImg = proc.ComputeGaussianGradient(img, sigma);
    RealImage::Pointer gradMag = proc.ComputeMagnitudeMap(gradImg);
    string gradMagMap = parser.GetString("--gradmag");
    if (gradMagMap != "") {
        io.WriteImage("gradmag.nrrd", gradMag);
    }
    cout << "Gradient Computation Done.. " << t.getElapsedTimeInSec() << endl;

    PolyDataPointer sphere = CreateSphere(parser.GetStringAsReal("--phi", 100), parser.GetStringAsReal("--theta", 100));

    __vtk(KdTreePointLocator);
    KdTreePointLocatorPointer pointsFinder  = KdTreePointLocatorPointer::New();
    pointsFinder->SetDataSet(sphere);
    pointsFinder->SetTolerance(0.0001);
    pointsFinder->BuildLocator();

    __vtk(IntArray);
    IntArrayPointer histogram = IntArrayPointer::New();
    histogram->SetName("GradientDensity");
    histogram->SetNumberOfValues(sphere->GetNumberOfPoints());
    for (int i = 0; i < sphere->GetNumberOfPoints(); i++) {
        histogram->SetValue(i, 0);
    }

    itk::ImageRegionIteratorWithIndex<GradientImage> gradIter(gradImg,gradImg->GetBufferedRegion());
    RealImage::SizeType sz = gradImg->GetBufferedRegion().GetSize();
    int numPixels = sz[0]*sz[1]*sz[2];

    gradIter.GoToBegin();


    int n = 0;
    VNLDoubleVector p(3);
    while (!gradIter.IsAtEnd()) {
        GradientPixel g = gradIter.Get();
        fordim(k) {
            p[k] = g[k];
        }
        double m = p.magnitude();
        if (m < thresh) {
            ++n;
            ++gradIter;
            continue;
        }
        p.normalize();
        int id = pointsFinder->FindClosestPoint(p.data_block());
        if (id >= 0 && id < sphere->GetNumberOfPoints()) {
            histogram->SetValue(id, histogram->GetValue(id) + 1);
        }
        ++gradIter;
        ++n;
        if (n % 100000 == 0) {
            printf("%5.3f\n", (n*100.0/numPixels));
        }
    }

    sphere->GetPointData()->AddArray(histogram);
    vtk_write_polydata(args[1].c_str(), sphere);
}