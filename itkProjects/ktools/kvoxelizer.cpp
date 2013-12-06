//
//  kvoxelizer.cpp
//  ktools
//
//  Created by Joohwi Lee on 12/5/13.
//
//

#include "kvoxelizer.h"
#include "piImageDef.h"
#include "piImageIO.h"
#include "piOptions.h"

#include "vtkio.h"

#include <vtkPolyData.h>
#include <vtkVoxelModeller.h>
#include <vtkPolyDataReader.h>
#include <vtkMetaImageWriter.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkSphereSource.h>
#include <vtkMetaImageWriter.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>

using namespace std;
using namespace pi;

int main1(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "kvoxelizer input-vtk output-image" << endl;
        return 0;
    }

    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    reader->SetFileName(argv[1]);
    vtkVoxelModeller* modeler = vtkVoxelModeller::New();
    modeler->SetInputConnection(reader->GetOutputPort());
    modeler->SetSampleDimensions(20, 20, 20);
    modeler->SetModelBounds(reader->GetOutput()->GetBounds());
    modeler->SetScalarTypeToChar();
    modeler->SetForegroundValue(1);
    modeler->SetBackgroundValue(0);

    modeler->Update();

    vtkMetaImageWriter* imageWriter = vtkMetaImageWriter::New();
    imageWriter->SetFileName(argv[2]);
    imageWriter->SetInputConnection(modeler->GetOutputPort());
    imageWriter->Write();

}



/// perform scan conversion
/// [input-vtk] [reference-image] [output-image]
///
int runScanConversion(pi::Options& opts, pi::StringVector& args) {
    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    reader->SetFileName(args[0].c_str());
    reader->Update();

    vtkPolyData*pd = reader->GetOutput();

    // point flipping
    for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
        double p[3];
        pd->GetPoint(i, p);
        p[0] = -p[0];
        p[1] = -p[1];
        pd->GetPoints()->SetPoint(i, p);
    }

    vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();

    ImageIO<RealImage> imageIO;
    RealImage::Pointer refImage = imageIO.ReadImage(args[1].c_str());


    // compute bounding box
    RealImage::RegionType region = refImage->GetBufferedRegion();
    RealImage::IndexType lowerIndex = region.GetIndex();
    RealImage::IndexType upperIndex = region.GetUpperIndex();

    RealImage::PointType lowerPoint, upperPoint;
    refImage->TransformIndexToPhysicalPoint(lowerIndex, lowerPoint);
    refImage->TransformIndexToPhysicalPoint(upperIndex, upperPoint);

    // mesh bounds
    double bounds[6];

    // image bounds
    bounds[0] = lowerPoint[0];
    bounds[1] = upperPoint[0];
    bounds[2] = lowerPoint[1];
    bounds[3] = upperPoint[1];
    bounds[4] = lowerPoint[2];
    bounds[5] = upperPoint[2];


    // print bounds
    for (int i = 0; i < 6; i++) {
        cout << bounds[i] << ", ";
    }
    cout << endl;


    // make the same spacing as refImage
    double spacing[3]; // desired volume spacing
    for (int i = 0; i < 3; i++) {
        spacing[i] = refImage->GetSpacing()[i];
    }
    whiteImage->SetSpacing(spacing);

    // compute dimensions
    int dim[3];
    for (int i = 0; i < 3; i++)
    {
        dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i])) + 1;
    }
    whiteImage->SetDimensions(dim);
    whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

    double origin[3];
    origin[0] = bounds[0] + spacing[0] / 2;
    origin[1] = bounds[2] + spacing[1] / 2;
    origin[2] = bounds[4] + spacing[2] / 2;
    whiteImage->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
    whiteImage->SetScalarTypeToUnsignedChar();
    whiteImage->AllocateScalars();
#else
    whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif
    // fill the image with foreground voxels:
    unsigned char inval = 255;
    unsigned char outval = 0;
    vtkIdType count = whiteImage->GetNumberOfPoints();
    for (vtkIdType i = 0; i < count; ++i)
    {
        whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
    }

    // polygonal data --> image stencil:
    vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc =
    vtkSmartPointer<vtkPolyDataToImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
    pol2stenc->SetInput(pd);
#else
    pol2stenc->SetInputData(pd);
#endif
    pol2stenc->SetOutputOrigin(origin);
    pol2stenc->SetOutputSpacing(spacing);
    pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
    pol2stenc->Update();

    // cut the corresponding white image and set the background:
    vtkSmartPointer<vtkImageStencil> imgstenc =
    vtkSmartPointer<vtkImageStencil>::New();
#if VTK_MAJOR_VERSION <= 5
    imgstenc->SetInput(whiteImage);
    imgstenc->SetStencil(pol2stenc->GetOutput());
#else
    imgstenc->SetInputData(whiteImage);
    imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
#endif
    imgstenc->ReverseStencilOff();
    imgstenc->SetBackgroundValue(outval);
    imgstenc->Update();

    vtkSmartPointer<vtkMetaImageWriter> writer =
    vtkSmartPointer<vtkMetaImageWriter>::New();
    writer->SetFileName("SphereVolume.mhd");
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(imgstenc->GetOutput());
#else
    writer->SetInputData(imgstenc->GetOutput());
#endif
    writer->Write();
    
    return EXIT_SUCCESS;
}


/// create a new image which has the same attributes with the reference image
///
int runPointMarks(pi::Options& opts, pi::StringVector& args) {
    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    reader->SetFileName(args[0].c_str());
    reader->Update();

    vtkPolyData*pd = reader->GetOutput();

    // load an image
    ImageIO<RealImage> imageIO;
    RealImage::Pointer refImage = imageIO.ReadImage(args[1].c_str());


    // create an empty image
    ImageIO<LabelImage> labelIO;
    LabelImage::Pointer markImage = labelIO.NewImageS<RealImage>(refImage);
    markImage->FillBuffer(0);

    // point flipping
    for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
        RealImage::PointType p;
        pd->GetPoint(i, p.GetDataPointer());
        p[0] = -p[0];
        p[1] = -p[1];
        pd->GetPoints()->SetPoint(i, p.GetDataPointer());

        // point to index
        RealImage::IndexType idx;
        refImage->TransformPhysicalPointToIndex(p, idx);

        // mark image
        markImage->SetPixel(idx, 255);
    }

    labelIO.WriteImage(args[2], markImage);
    return EXIT_SUCCESS;
}

/**
 * This program generates a sphere (closed surface, vtkPolyData) and converts it into volume
 * representation (vtkImageData) where the foreground voxels are 1 and the background voxels are
 * 0. Internally vtkPolyDataToImageStencil is utilized. The resultant image is saved to disk
 * in metaimage file format (SphereVolume.mhd).
 */
int main(int argc, char * argv[])
{
    /*
    vtkSmartPointer<vtkSphereSource> sphereSource =
    vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetRadius(20);
    sphereSource->SetPhiResolution(30);
    sphereSource->SetThetaResolution(30);
    vtkSmartPointer<vtkPolyData> pd = sphereSource->GetOutput();
    sphereSource->Update();
     */

    Options opts;
    opts.addOption("--v", "voxelize a given mesh into a binary image [mesh] [refimage] [outimage]", SO_NONE);
    opts.addOption("--m", "create a binary image with a mark per point [mesh] [refimage] [outimage]", SO_NONE);
    StringVector args = opts.ParseOptions(argc, argv, NULL);


    if (opts.GetBool("--v")) {
        runScanConversion(opts, args);
    } else if (opts.GetBool("--m")) {
        runPointMarks(opts, args);
    }



}

