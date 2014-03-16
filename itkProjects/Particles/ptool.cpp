//
//  ptool.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 6/20/13.
//
//

#include "ptool.h"
#include "piOptions.h"
#include "piImageDef.h"
#include "piImageIO.h"

#include <itkZeroCrossingBasedEdgeDetectionImageFilter.h>
#include <itkPointSet.h>
#include <itkManifoldParzenWindowsPointSetFunction.h>
#include <itkCannyEdgeDetectionImageFilter.h>

using namespace pi;


void SVMInput(StringVector& args) {
    // read label and intensity images
    // using sliding window compute intensity patches, then label 1 if it is boundary of the label
    // how to detect boundaries?
    static ImageIO<RealImage2> realIO;
    static ImageIO<LabelImage2> labelIO;
    
    RealImage2::Pointer realImage = realIO.ReadImage(args[0]);
    LabelImage2::Pointer labelImage;
    if (args.size() > 1) {
        labelImage = labelIO.ReadImage(args[1]);
        itk::ZeroCrossingImageFilter<LabelImage2, LabelImage2>::Pointer filter =
            itk::ZeroCrossingImageFilter<LabelImage2, LabelImage2>::New();
        filter->SetInput(labelImage);
        filter->SetBackgroundValue(0);
        filter->SetForegroundValue(1);
        filter->Update();
        labelImage = filter->GetOutput();

        labelIO.WriteImage(args[2], labelImage);
    }

    RealImage2::SizeType szIm = realImage->GetBufferedRegion().GetSize();
    std::vector<float> values;
    values.resize(25);

    for (int j = 0; j < szIm[1] - 5; j++) {
        for (int i = 0; i < szIm[0] - 5; i++) {
            RealImage2::RegionType region;
            region.SetIndex(0, i);
            region.SetIndex(1, j);
            region.SetSize(0, 5);
            region.SetSize(1, 5);

            itk::ImageRegionConstIteratorWithIndex<RealImage2> iter(realImage, region);
            int label = -1;
            iter.GoToBegin();
            for (int k = 1; !iter.IsAtEnd(); k++) {
                float v = iter.Get();
                values[k-1] = v;
                if (k == 13) {
                    LabelImage2::IndexType idx = iter.GetIndex();
                    if (labelImage->GetPixel(idx) > 0) {
                        label = 1;
                    };
                }
                ++iter;
            }

            cout << label << " ";
            for (int k = 0; k < values.size(); k++) {
                cout << (k+1) << ":" << values[k] << " ";
            }
            cout << endl;
        }
    }
}

void SVMOutput(StringVector& args) {
    ImageIO<LabelImage2> io;
    LabelImage2::Pointer label = io.ReadImage(args[0]);
    LabelImage2::Pointer newImage = io.NewImage(label);
    LabelImage2::SizeType szIm = label->GetBufferedRegion().GetSize();

    string line;
    ifstream fin(args[1].c_str());
    for (int k = 0; fin.good(); k++) {
        int i = k % (szIm[0] - 5);
        int j = k / (szIm[0] - 5);


        getline(fin, line);
        if (line == "-1") {
        } else {
            cout << i << "," << j << endl;

            LabelImage2::IndexType idx;
            idx[0] = i + 2;
            idx[1] = j + 2;
            newImage->SetPixel(idx, 10);
        }
    }
    io.WriteImage(args[2], newImage);
}

void IdentityDirection(StringVector& args) {
    ImageIO<RealImage> io;
    ImageInfo inputImageInfo;
    RealImage::Pointer inputImage = io.ReadCastedImage(args[0], inputImageInfo);
    RealImage::DirectionType dir;
    dir.SetIdentity();
    inputImage->SetDirection(dir);
    io.WriteCastedImage(args[1], inputImage, inputImageInfo.componenttype);
}

void EnumPoints(StringVector& args) {
    for (int i = 0; i < args.size(); i++) {
        ImageIO<RealImage> io;
        RealImage::Pointer inputImage = io.ReadImage(args[i]);

        inputImage->Print(cout);
        
        RealImageIteratorType iter(inputImage, inputImage->GetBufferedRegion());
        iter.GoToBegin();
        while (!iter.IsAtEnd()) {
            RealImage::IndexType idx = iter.GetIndex();
            RealImage::PointType point;
            inputImage->TransformIndexToPhysicalPoint(idx, point);
            cout << idx << " => " << point << endl;
            ++iter;
        }
    }
}

void ComputeUnion(StringVector& args) {
    std::vector<float> unionBox;
    unionBox.resize(6);
    for (int i = 0; i < args.size(); i++) {
        ImageIO<LabelImage> io;
        LabelImage::Pointer labelImage = io.ReadImage(args[i]);
        std::vector<float> boundingBox = __labelImageTools.computeBoundingBox(labelImage, 1);
        if (i == 0) {
            unionBox = boundingBox;
        } else {
            for (int j = 0; j < 3; j++) {
                unionBox[j] = std::min(boundingBox[j], unionBox[j]);
            }
            for (int j = 3; j < 6; j++) {
                unionBox[j] = std::max(boundingBox[j], unionBox[j]);
            }
        }
        for (int j = 0; j < unionBox.size(); j++) {
            cout << unionBox[j] << ", ";
        }
        cout << endl;
    }
}

void DoParzenTest(StringVector& args) {
    typedef itk::PointSet<float,2> PointSetType;
    typedef PointSetType::PointType PointType;
    typedef PointSetType::PointsContainerPointer PointsContainerPointer;

    PointSetType::Pointer pointSet = PointSetType::New();
    PointsContainerPointer points = pointSet->GetPoints();

    for (int i = 0; i < 100; i++) {
        PointType p;
        p[0] = std::cos(i / 100.0 * M_PI * 2);
        p[1] = std::sin(i / 100.0 * M_PI * 2);
        points->InsertElement(i, p);

    }

    typedef itk::ManifoldParzenWindowsPointSetFunction<PointSetType> ParzenFuncType;
    ParzenFuncType::Pointer parzenFunc = ParzenFuncType::New();

    float kernelSigma = atof(args[0].c_str());
    float regularizationSigma = atof(args[1].c_str());

    cout << "Kernel Sigma: " << kernelSigma << endl;
    cout << "Regularization Sigma: " << regularizationSigma << endl;

    parzenFunc->SetKernelSigma(kernelSigma);
    parzenFunc->SetRegularizationSigma(regularizationSigma);
    parzenFunc->UseAnisotropicCovariancesOff();
    parzenFunc->SetInputPointSet(pointSet);

    typename PointSetType::PointIdentifier numPoints = parzenFunc->GetInputPointSet()->GetNumberOfPoints();
    cout << "# of points: " << numPoints << endl;

    static ImageIO<RealImage2> realIO;
    RealImage2::Pointer realImage = realIO.NewImageT(301, 301, 1);

    RealImage2::SpacingType spacing;
    spacing[0] = 0.01;
    spacing[1] = 0.01;
    realImage->SetSpacing(spacing);

    RealImage2::PointType origin;
    origin[0] = -150*spacing[0];
    origin[1] = -150*spacing[1];
    realImage->SetOrigin(origin);

    itk::ImageRegionIteratorWithIndex<RealImage2> regionIter(realImage, realImage->GetBufferedRegion());
    regionIter.GoToBegin();
    for (; !regionIter.IsAtEnd(); ++regionIter) {
        RealImage2::IndexType idx = regionIter.GetIndex();
        RealImage2::PointType idxPoint;
        realImage->TransformIndexToPhysicalPoint(idx, idxPoint);

        realImage->SetPixel(idx, parzenFunc->Evaluate(idxPoint));
    }
    realIO.WriteImage(args[2], realImage);

    return;
}


void CannyEdge(StringVector& args) {
    static ImageIO<RealImage2> realIO;
    RealImage2::Pointer realImage = realIO.ReadImage(args[0]);

    typedef itk::CannyEdgeDetectionImageFilter<RealImage2, RealImage2> FilterType;
    FilterType::Pointer filter = FilterType::New();
    FilterType::ArrayType var;
    var[0] = 0.15;
    var[1] = 0.15;
    filter->SetInput(realImage);
    filter->SetVariance(var);
    filter->Update();
    RealImage2::Pointer edgeImage = filter->GetOutput();
    realIO.WriteImage(args[1], edgeImage);
}

void DoParzenSample(StringVector& args) {
    static ImageIO<RealImage2> realIO;
    RealImage2::Pointer realImage = realIO.ReadImage(args[0]);
    float kernelSigma = atof(args[1].c_str());
    float regularizationSigma = atof(args[2].c_str());



    typedef itk::PointSet<float,2> PointSetType;
    typedef PointSetType::PointType PointType;
    typedef PointSetType::PointsContainerPointer PointsContainerPointer;

    PointSetType::Pointer pointSet = PointSetType::New();
    PointsContainerPointer points = pointSet->GetPoints();

    itk::ImageRegionIteratorWithIndex<RealImage2> regionIter(realImage, realImage->GetBufferedRegion());
    regionIter.GoToBegin();
    RealImage2::PointType point;

    int i = 0;
    for (; !regionIter.IsAtEnd(); ++regionIter) {
        RealImage2::IndexType idx = regionIter.GetIndex();
        if (regionIter.Get() > 0) {
            realImage->TransformIndexToPhysicalPoint(idx, point);
            points->InsertElement(i++, point);
            cout << point << endl;
        }
    }

    typedef itk::ManifoldParzenWindowsPointSetFunction<PointSetType> ParzenFuncType;
    ParzenFuncType::Pointer parzenFunc = ParzenFuncType::New();


    cout << "Kernel Sigma: " << kernelSigma << endl;
    cout << "Regularization Sigma: " << regularizationSigma << endl;

    parzenFunc->SetKernelSigma(kernelSigma);
    parzenFunc->SetRegularizationSigma(regularizationSigma);
    parzenFunc->SetInputPointSet(pointSet);

    typename PointSetType::PointIdentifier numPoints = parzenFunc->GetInputPointSet()->GetNumberOfPoints();
    cout << "# of points: " << numPoints << endl;

    RealImage2::Pointer outputImage = realIO.NewImage(realImage);
    itk::ImageRegionIteratorWithIndex<RealImage2> outputIter(outputImage, outputImage->GetBufferedRegion());
    outputIter.GoToBegin();
    for (; !outputIter.IsAtEnd(); ++outputIter) {
        RealImage2::IndexType idx = outputIter.GetIndex();
        RealImage2::PointType idxPoint;
        outputImage->TransformIndexToPhysicalPoint(idx, idxPoint);
        outputImage->SetPixel(idx, parzenFunc->Evaluate(idxPoint));
    }
    

    realIO.WriteImage(args[3], outputImage);
    
    return;
}

int main(int argc, char* argv[]) {
    CSimpleOpt::SOption options[] ={
        { 1, "--union", SO_NONE },
        { 2, "--crop", SO_REQ_SEP },
        { 3, "-o", SO_REQ_SEP },
        { 4, "--enumPoint", SO_NONE },
        { 5, "--identity-dir", SO_NONE },
        { 6, "--svm-input", SO_NONE },
        { 7, "--svm-output", SO_NONE },
        { 8, "--parzen-test", SO_NONE },
        { 9, "--canny-edge", SO_NONE },
        { 10, "--parzen-sample", SO_NONE },
        SO_END_OF_OPTIONS
    };

    Options parser;
    StringVector& args = parser.ParseOptions(argc, argv, options);

    if (parser.GetBool("--union")) {
        ComputeUnion(args);
    } else if (parser.GetBool("--enumPoint")) {
        EnumPoints(args);
    } else if (parser.GetBool("--identity-dir")) {
        IdentityDirection(args);
    } else if (parser.GetBool("--svm-input")) {
        SVMInput(args);
    } else if (parser.GetBool("--svm-output")) {
        SVMOutput(args);
    } else if (parser.GetBool("--parzen-test")) {
        DoParzenTest(args);
    } else if (parser.GetBool("--canny-edge")) {
        CannyEdge(args);
    } else if (parser.GetBool("--parzen-sample")) {
        DoParzenSample(args);
    }
}