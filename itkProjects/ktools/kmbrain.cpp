//
//  ktest.cpp
//  ktools
//
//  Created by Joohwi Lee on 3/8/13.
//
//

#include "kmbrain.h"

#include <stdio.h>
#include <vtkKdTreePointLocator.h>
#include <vtkPointData.h>
#include "vtkUnstructuredGrid.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkXMLUnstructuredGridWriter.h"

#include "piImageDef.h"
#include "piImageProcessing.h"
#include "piOptions.h"
#include "piTimer.h"
#include "piImageIO.h"

#include "itkImageRegionMultidimensionalSplitter.h"

using namespace std;
using namespace pi;

int main(int argc, char* argv[]) {
    CSimpleOpt::SOption specs[] = {
        SO_END_OF_OPTIONS
    };

    Options parser;
    StringVector args = parser.ParseOptions(argc, argv, specs);

    if (args.size() < 1) {
        cout << argv[0] << " input-image ouptut-vtu" << endl;
        return 1;
    }

    DataReal sigma = parser.GetStringAsReal("-s", 0.3);
    ImageIO<RealImage> io;
    RealImage::Pointer img = io.ReadCastedImage(args[0]);

    ImageProcessing proc;
    GradientImage::Pointer gradImg = proc.ComputeGaussianGradient(img, sigma);
    GradientImage::RegionType region = gradImg->GetBufferedRegion();

    typedef itk::ImageRegionMultidimensionalSplitter<DIMENSIONS> SplitType;
    SplitType::Pointer split = SplitType::New();

    int nRegions = split->GetNumberOfSplits(region, 1000);
    cout << "n Regions = " << nRegions << endl;

    vtkPoints* points = vtkPoints::New();
    points->SetNumberOfPoints(nRegions);

    vtkDoubleArray* array = vtkDoubleArray::New();
    array->SetNumberOfComponents(3);
    array->SetNumberOfTuples(nRegions);
    array->SetName("MeanGradients");

    for (int i = 0; i < nRegions; i++) {
        GradientImage::RegionType subRegion = split->GetSplit(i, nRegions, region);
        itk::ImageRegionConstIterator<GradientImage> iter(gradImg, subRegion);
        GradientPixel m;
        m.Fill(0);
        int n = 0;
        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter, n++) {
            GradientPixel g = iter.Get();
            m += g;
        }
        m /= n;

        IntIndex uIdx = subRegion.GetUpperIndex();
        IntIndex lIdx = subRegion.GetIndex();
        IntIndex cIdx;
        fordim (k) {
            cIdx[k] = (uIdx[k] + lIdx[k]) / 2.0;
        }

        points->SetPoint(i, cIdx[0], cIdx[1], cIdx[2]);
        array->SetTuple3(i, m[0], m[1], m[2]);
    }

    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    grid->SetPoints(points);
    grid->GetPointData()->AddArray(array);

    vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetInput(grid);
    writer->SetCompressorTypeToZLib();
    writer->EncodeAppendedDataOn();
    writer->SetFileName(args[1].c_str());
    writer->Write();
}