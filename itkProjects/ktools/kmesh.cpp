
//
//  kmesh.cpp
//  ktools
//
//  Created by Joohwi Lee on 12/5/13.
//
//

#include "kmesh.h"

#include <set>
#include <iostream>

#include "piOptions.h"
#include "vtkio.h"

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkIdList.h>
#include <vtkFloatArray.h>
#include <vtkMath.h>
#include <vtkAppendPolyData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkStreamTracer.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkPolyLine.h>
#include <vtkMath.h>
#include <vtkSmoothPolyDataFilter.h>


#include <itkImage.h>
#include <itkVectorNearestNeighborInterpolateImageFunction.h>
#include "piImageIO.h"
#include "kimage.h"
#include "kstreamtracer.h"



#include <vnl/vnl_vector.h>


using namespace std;
using namespace pi;


// append polydatas into one
void runAppendData(Options& opts, StringVector& args) {
    vtkAppendPolyData* appender = vtkAppendPolyData::New();
    vtkIO io;

    // read all files
    for (int i = 0; i < args.size(); i++) {
        vtkPolyData* data = io.readFile(args[i]);
        appender->AddInput(data);
    }

    appender->Update();
    io.writeFile(opts.GetString("-appendData"), appender->GetOutput());

    return;
}

// add scalar value to a mesh
void runImportScalars(Options& opts, StringVector& args) {
    vtkIO io;

    // read polydata
    vtkPolyData* poly = io.readFile(args[0]);
    vtkFloatArray* scalar = vtkFloatArray::New();
    scalar->SetNumberOfValues(poly->GetNumberOfPoints());
    scalar->SetName(opts.GetString("-scalarName").c_str());

    ifstream file(args[1].c_str());
    for (int i = 0; i < poly->GetNumberOfPoints() && !file.eof(); i++) {
        float value = 0;
        file >> value;
        scalar->SetValue(i, value);
    }

    poly->GetPointData()->AddArray(scalar);

    io.writeFile(args[2], poly);
}


/// @brief export scalar values to a text file
void runExportScalars(Options& opts, StringVector& args) {
    vtkIO io;

    /// - Read polydata
    vtkPolyData* poly = io.readFile(args[0]);

    /// - Check if file is loaded
    if (poly == NULL) {
        cout << "can't read file: " << args[0];
        return;
    }

    bool isPointData = true;
    /// - Find the scalar attribute given as '-scalarName'
    vtkDataArray* scalar = poly->GetPointData()->GetScalars(opts.GetString("-scalarName").c_str());
    /// - Check if the scalar exists
    if (scalar == NULL) {
        /// - Try with cell arrays
        scalar = poly->GetCellData()->GetScalars(opts.GetString("-scalarName").c_str());
        if (scalar == NULL) {
            cout << "can't find the scalar attribute: " << opts.GetString("-scalarName") << endl;
            return;
        }
        isPointData = false;
    }

    ofstream file(args[1].c_str());

    int nScalars = scalar->GetNumberOfTuples();
    for (int i = 0; i < nScalars; i++) {
        file << scalar->GetTuple1(i) << endl;
    }
    file.close();
}



/// @brief Perform smoothing on manifold by iterative averaging as used in FreeSurfer
void runScalarSmoothing(Options& opts, StringVector& args) {
    // number of iterations and sigma affects the smoothing results
    vtkIO io;
    vtkMath* math = vtkMath::New();

    // sigma for smoothing
    double sigma2 = opts.GetStringAsReal("-sigma", 1);
    sigma2 *= sigma2;

    // number of iterations
    int numIters = opts.GetStringAsInt("-iter", 1);


    // read polydata
    vtkPolyData* poly = io.readFile(args[0]);

    // access data
    string scalarName = opts.GetString("-scalarName");
    vtkDataArray* scalars = poly->GetPointData()->GetScalars(scalarName.c_str());
    if (scalars == NULL) {
        cout << "can't find scalars: " << scalarName << endl;
        return;
    }

    // copy the scalars to iteratively apply smoothing
    vtkFloatArray* data = vtkFloatArray::New();
    data->DeepCopy(scalars);

    // prepare new data array
    vtkFloatArray* newData = vtkFloatArray::New();

    string outputScalarName = opts.GetString("-outputScalarName", "smoothed_" + scalarName);
    newData->SetName(outputScalarName.c_str());
    newData->SetNumberOfTuples(data->GetNumberOfTuples());
    poly->GetPointData()->AddArray(newData);


    // check if the scalar array exists
    if (data == NULL) {
        cout << "can't access scalar array: " << scalarName << endl;
        return;
    }

    // iterate over all points
    vtkIdList* cellIds = vtkIdList::New();
    vtkIdList* ptIds = vtkIdList::New();
    std::set<int> ptSet;

    // build cells
    poly->BuildCells();
    poly->BuildLinks();

    for (int n = 0; n < numIters; n++) {
        cout << "Iter: " << n << endl;
        for (int i = 0; i < poly->GetNumberOfPoints(); i++) {
            double center[3];
            poly->GetPoint(i, center);

            // collect neighbor cells
            ptSet.clear();

            cellIds->Reset();
            poly->GetPointCells(i, cellIds);

            // iterate over neighbor cells
            for (int j = 0; j < cellIds->GetNumberOfIds(); j++) {
                int cellId = cellIds->GetId(j);
                ptIds->Reset();

                // collect cell points
                poly->GetCellPoints(cellId, ptIds);

                // iterate over all cell points
                for (int k = 0; k < ptIds->GetNumberOfIds(); k++) {
                    int ptId = ptIds->GetId(k);
                    ptSet.insert(ptId);
                }
            }

            // iterate over all neighbor points
            std::set<int>::iterator iter = ptSet.begin();

            // compute weight
            vnl_vector<float> weights;
            weights.set_size(ptSet.size());

            for (int j = 0; iter != ptSet.end(); iter++, j++) {
                int ptId = *iter;
                double neighbor[3];
                poly->GetPoint(ptId, neighbor);

                double dist2 = math->Distance2BetweenPoints(center, neighbor);

                // apply the heat kernel with the sigma
                weights[j] = exp(-dist2/sigma2);
            }

            // add one for the center
            double weightSum = weights.sum() + 1;
            weights /= weightSum;


            // iterate over neighbors and compute weighted sum
            double smoothedValue = data->GetTuple1(i) / weightSum;
            iter = ptSet.begin();

            // compute the weighted averge
            for (uint j = 0; j < ptSet.size(); j++, iter++) {
                int ptId = *iter;
                int value = data->GetTuple1(ptId);
                smoothedValue += (value * weights[j]);
            }
            newData->SetTuple1(i, smoothedValue);
        }


        // prepare next iteration by copying newdata to data
        data->DeepCopy(newData);
    }

    // write to file
    io.writeFile(args[1], poly);
}



/// @brief Convert an ITK image to a VTKImageData
void runConvertITK2VTI(Options& opts, StringVector& args) {
    if (args.size() < 2) {
        cout << "requires input-image-file and output-vti-file" << endl;
        return;
    }
    
    int attrDim = opts.GetStringAsInt("-attrDim", 1);
    string scalarName = opts.GetString("-scalarName", "Intensity");
    
    string input = args[0];
    string output = args[1];
    
    /// - Read an image data
    vtkImageData* outputData = vtkImageData::New();
    if (attrDim == 1) {
        ConvertImageT<ImageType>(input, outputData, scalarName.c_str(), 1);
    } else if (attrDim == 3) {
        ConvertImageT<VectorImageType>(input, outputData, scalarName.c_str(), attrDim);
    }
    
    vtkXMLImageDataWriter* w = vtkXMLImageDataWriter::New();
    w->SetFileName(output.c_str());
    w->SetDataModeToAppended();
    w->EncodeAppendedDataOff();
    w->SetCompressorTypeToZLib();
    w->SetDataModeToBinary();
    
    w->SetInput(outputData);
    w->Write();
}


/// @brief Convert an itk image file to vtkUnstructuredGrid
void runConvertITK2VTU(Options& opts, StringVector& args) {
    if (args.size() < 2) {
        cout << "requires input-image-file and output-vti-file" << endl;
        return;
    }

    int attrDim = opts.GetStringAsInt("-attrDim", 1);
    string scalarName = opts.GetString("-scalarName", "Intensity");

    string input = args[0];
    string output = args[1];

    string maskImageFile = opts.GetString("-maskImage");

    typedef itk::Image<ushort,3> MaskImageType;
    ImageIO<MaskImageType> maskIO;
    MaskImageType::Pointer maskImage = maskIO.ReadCastedImage(maskImageFile);

    /// - Read an image data
    vtkUnstructuredGrid* outputData = vtkUnstructuredGrid::New();
    if (attrDim == 1) {
        ConvertImageT<ImageType, MaskImageType>(input, outputData, maskImage, scalarName.c_str(), 1);
    } else if (attrDim == 3) {
        ConvertVectorImageT<VectorImageType, MaskImageType>(input, outputData, maskImage, scalarName.c_str(), attrDim);
    }

    vtkXMLUnstructuredGridWriter* w = vtkXMLUnstructuredGridWriter::New();
    w->SetFileName(output.c_str());
    w->SetDataModeToAppended();
    w->EncodeAppendedDataOff();
    w->SetCompressorTypeToZLib();
    w->SetDataModeToBinary();

    w->SetInput(outputData);
    w->Write();
}

bool endswith(std::string str, std::string substr) {
    size_t i = str.rfind(substr);
    return (i != string::npos) && (i == (str.length() - substr.length()));
}

/// @brief Execute the stream tracer
void runStreamTracer(Options& opts, StringVector& args) {
    string inputVTUFile = args[0];
    string inputPointsFile = args[1];
    string outputVTKFile = args[2];
    string outputPointFile = args[3];
    bool zRotate = opts.GetBool("-zrotate", false);


    vtkDataSet* inputData;

    // FIXME - create a dataset reader
    if (endswith(inputVTUFile, string(".vtu"))) {
        vtkXMLUnstructuredGridReader* reader = vtkXMLUnstructuredGridReader::New();
        reader->SetFileName(inputVTUFile.c_str());
        reader->Update();

        vtkUnstructuredGrid* inputVTU = reader->GetOutput();
        inputData = inputVTU;
    } else if (endswith(inputVTUFile, ".vti")) {
        vtkXMLImageDataReader* reader = vtkXMLImageDataReader::New();
        reader->SetFileName(inputVTUFile.c_str());
        reader->Update();

        vtkImageData* inputVTI = reader->GetOutput();
        inputData = inputVTI;
    }

    vtkIO vio;
    vtkPolyData* inputPoints = vio.readFile(inputPointsFile);
    vtkPoints* points = inputPoints->GetPoints();

    /// - Converting the input points to the image coordinate
    const int nInputPoints = inputPoints->GetNumberOfPoints();
    for (int i = 0; i < nInputPoints; i++) {
        double p[3];
        points->GetPoint(i, p);
        // FixMe: Do not use a specific scaling factor
        if (zRotate) {
            p[0] = -p[0];
            p[1] = -p[1];
            p[2] = p[2];
        }
        points->SetPoint(i, p);
    }
    inputPoints->SetPoints(points);

    /// - Set up tracer (Use RK45, both direction, initial step 0.05, maximum propagation 500
    StreamTracer* tracer = StreamTracer::New();
    tracer->SetInput(inputData);
//    tracer->SetSource(inputPoints);
    double seedPoint[3];
    inputPoints->GetPoint(24745, seedPoint);
    tracer->SetStartPosition(seedPoint);
    tracer->SetIntegratorTypeToRungeKutta45();

    if (opts.GetString("-traceDirection") == "both") {
        tracer->SetIntegrationDirectionToBoth();
    } else if (opts.GetString("-traceDirection") == "backward") {
        tracer->SetIntegrationDirectionToBackward();
    } else {
        tracer->SetIntegrationDirectionToForward();
    }

    tracer->SetInterpolatorTypeToDataSetPointLocator();
    tracer->SetMaximumPropagation(500);
    tracer->SetInitialIntegrationStep(0.05);
    tracer->Update();



    vtkPolyData* streamLines = tracer->GetOutput();
    // loop over the cell and compute the length
    int nCells = streamLines->GetNumberOfCells();
    cout << "# of cells: " << nCells << endl;


    /// - Prepare the output as a scalar array
    vtkDoubleArray* streamLineLength = vtkDoubleArray::New();
    streamLineLength->SetNumberOfValues(nCells);
    streamLineLength->SetName("Length");
    streamLineLength->SetNumberOfComponents(1);


    if (true) {
        cout << "Computing the length of each stream line ..." << endl;
        for (int i = 0; i < nCells; i++) {
            vtkPolyLine* polyLine = vtkPolyLine::SafeDownCast(streamLines->GetCell(i));
            double length = 0;
            if (polyLine) {
                length = sqrt(polyLine->GetLength2());
            } else {
                vtkLine* line = vtkLine::SafeDownCast(streamLines->GetCell(i));
                if (line) {
                    length = sqrt(line->GetLength2());
                } else {
                    cout << "Can't find polyline: " << i << "; Type = " << streamLines->GetCell(i)->GetCellType() << endl;
                }
//                exit(0);
            }
            streamLineLength->SetValue(i, length);
        }
        streamLines->GetCellData()->SetScalars(streamLineLength);
    }


    /// - Prepare the output for the input points
    vtkDoubleArray* streamLineLengthPerPoint = vtkDoubleArray::New();
    streamLineLengthPerPoint->SetNumberOfTuples(nInputPoints);
    streamLineLengthPerPoint->SetName("Length");
    streamLineLengthPerPoint->SetNumberOfComponents(1);
    inputPoints->GetPointData()->SetScalars(streamLineLengthPerPoint);

    cout << "Assigning a length to each source vertex ..." << endl;
    vtkDataArray* seedIds = streamLines->GetCellData()->GetScalars("SeedId");
    if (seedIds) {
        for (int i = 0; i < nInputPoints; i++) {
            streamLineLengthPerPoint->SetValue(i, -1);
        }
        for (int i = 0; i < seedIds->GetNumberOfTuples(); i++) {
            int id = seedIds->GetTuple1(i);
            if (id > nInputPoints) {
                id = (id - nInputPoints);
            }
            double length2 = streamLineLength->GetValue(i);
            double length1 = streamLineLengthPerPoint->GetValue(id);
            if (length1 > 0) {
                double maxlength = std::max(length1, length2);
                streamLineLengthPerPoint->SetValue(id, maxlength);
            } else {
                streamLineLengthPerPoint->SetValue(id, length2);
            }
        }
    } else {
        cout << "Can't find SeedId" << endl;
    }


    vio.writeFile(outputPointFile, inputPoints);
    vio.writeXMLFile(outputVTKFile, streamLines);
}

/// @brief Apply a filter to each stream line
void runFilterStream(Options& opts, StringVector& args) {
    string inputStream = args[0];
    string inputSeeds = args[1];
    string outputStreamFile = args[2];
    string scalarName = opts.GetString("-scalarName");

    double lowThreshold = opts.GetStringAsReal("-thresholdMin", itk::NumericTraits<float>::min());
    double highThreshold = opts.GetStringAsReal("-thresholdMax", itk::NumericTraits<float>::max());

    vtkIO vio;
    vtkPolyData* streamLines = vio.readFile(inputStream);
    vtkPolyData* streamSeeds = vio.readFile(inputSeeds);

    vtkPolyData* outputStream = vtkPolyData::New();
    outputStream->SetPoints(streamLines->GetPoints());

    /// - Lookup SeedId and a given scalar array
    vtkDataArray* seedList = streamLines->GetCellData()->GetScalars("SeedId");
    vtkDataArray* seedScalars = streamSeeds->GetPointData()->GetScalars(scalarName.c_str());
    vtkCellArray* lines = vtkCellArray::New();
    vtkDoubleArray* filteredScalars = vtkDoubleArray::New();
    filteredScalars->SetName(scalarName.c_str());
    filteredScalars->SetNumberOfComponents(1);

    /// - Lookup a corresponding point scalar, apply threashold filter, and add to the new object
    for (int i = 0; i < seedList->GetNumberOfTuples(); i++) {
        int seedId = seedList->GetTuple1(i);
        double value = seedScalars->GetTuple1(seedId);

        if (lowThreshold <= value && value <= highThreshold) {
            lines->InsertNextCell(streamLines->GetCell(i));
            filteredScalars->InsertNextValue(value);
        }
    }

    outputStream->SetLines(lines);
    outputStream->GetCellData()->AddArray(filteredScalars);
    outputStream->BuildCells();
    outputStream->BuildLinks();

    vio.writeFile(outputStreamFile, outputStream);
}


/// @brief Fit a model into a binary image
void runFittingModel(Options& opts, StringVector& args) {
    if (args.size() < 3) {
        cout << "requires input-model input-image output-model" << endl;
        return;
    }
    string inputModelFile = args[0];
    string inputImageFile = args[1];
    string outputModelFile = args[2];

    vtkIO vio;
    vtkPolyData* inputModel = vio.readFile(inputModelFile);
    const int nPoints = inputModel->GetNumberOfPoints();

    /// Apply z-rotation
    for (int i = 0; i < nPoints; i++) {
        double point[3];
        inputModel->GetPoint(i, point);
        if (opts.GetBool("-zrotate")) {
            point[0] = -point[0];
            point[1] = -point[1];
        }
        inputModel->GetPoints()->SetPoint(i, point);
    }

    ImageIO<MaskImageType> itkIO;
    MaskImageType::Pointer maskImage = itkIO.ReadCastedImage(inputImageFile);

    // for test, rescale 10 times
    VectorImageType::SpacingType spacing = maskImage->GetSpacing();
    for (int i = 0; i < 3; i++) {
//        spacing[i] *= 10;
    }
    maskImage->SetSpacing(spacing);


    // test
    VectorImageType::Pointer distImage = ComputeDistanceMap(maskImage);
    typedef itk::VectorNearestNeighborInterpolateImageFunction<VectorImageType> InterpolatorType;

    ImageIO<VectorImageType> distIO;
    distIO.WriteImage("dist.mha", distImage);

    InterpolatorType::Pointer distInterp = InterpolatorType::New();
    distInterp->SetInputImage(distImage);


    // iterate over the input model and project to the boundary
    const int nIters = 50;
    for (int i = 0; i < nIters; i++) {
        // Compute laplacian smoothing by taking iterative average

        vtkSmoothPolyDataFilter* filter = vtkSmoothPolyDataFilter::New();
        filter->SetInput(inputModel);
        filter->SetNumberOfIterations(1);
        filter->Update();
        inputModel = filter->GetOutput();

        // projection
        for (int j = 0; j < nPoints; j++) {
            VectorImageType::PointType point, nextPoint;
            inputModel->GetPoint(j, point.GetDataPointer());
            VectorType offset = distInterp->Evaluate(point);
            for (int k = 0; k < 3; k++) {
                nextPoint = point + offset[k] * spacing[k] * spacing[k];
            }
            cout << point << " => " << nextPoint << endl;
            inputModel->GetPoints()->SetPoint(j, nextPoint.GetDataPointer());
        }
    }

    vio.writeFile(outputModelFile, inputModel);
}

int main(int argc, char * argv[])
{
    Options opts;
    opts.addOption("-exportScalars", "Export scalar values to a text file", "-exportScalars [in-mesh] [scalar.txt]", SO_NONE);
    opts.addOption("-importScalars", "Add scalar values to a mesh [in-mesh] [scalar.txt] [out-mesh]", SO_NONE);
    opts.addOption("-smoothScalars", "Gaussian smoothing of scalar values of a mesh. [in-mesh] [out-mesh]", SO_NONE);
    opts.addOption("-appendData", "Append input meshes into a single data [output-mesh]", SO_REQ_SEP);

    opts.addOption("-sigma", "sigma value [double]", SO_REQ_SEP);
    opts.addOption("-scalarName", "scalar name [string]", SO_REQ_SEP);
    opts.addOption("-outputScalarName", "scalar name for output [string]", SO_REQ_SEP);
    opts.addOption("-iter", "number of iterations [int]", SO_REQ_SEP);
    opts.addOption("-attrDim", "The number of components of attribute", "-attrDim 3 (vector)", SO_REQ_SEP);
    opts.addOption("-vti", "Convert an ITK image to VTI format (VTKImageData)", "-vti imageFile outputFile [-attrDim 3]", SO_NONE);
    opts.addOption("-vtu", "Convert an ITK image to VTU format (vtkUnstructuredGrid). This is useful when masking is needed.", "-vtu imageFile outputFile -maskImage maskImage", SO_NONE);
    opts.addOption("-maskImage", "A mask image for the use of -vtu", "-maskImage mask.nrrd", SO_REQ_SEP);
    opts.addOption("-traceStream", "Trace a stream line from a given point set", "-traceStream input-vtu-field input-vtk output-lines output-points", SO_NONE);
    opts.addOption("-traceDirection", "Choose the direction of stream tracing (both, forward, backward)", "-traceStream ... -traceDirection (both|forward|backward)", SO_REQ_SEP);
    opts.addOption("-zrotate", "Rotate all the points along the z-axis. Change the sign of x and y coordinate.", "-traceStream ... -zrotate", SO_NONE);
    opts.addOption("-filterStream", "Filter out stream lines which are lower than a given threshold", "-filterStream stream-line-input stream-seed-input stream-line-output -scalarName scalar -threshold xx", SO_NONE);
    opts.addOption("-thresholdMin", "Give a minimum threshold value for -filterStream", "-threshold 10 (select a cell whose attriubte is greater than 10)", SO_REQ_SEP);
    opts.addOption("-thresholdMax", "Give a maximum threshold value for -filterStream", "-threshold 10 (select a cell whose attriubte is lower than 10)", SO_REQ_SEP);
    opts.addOption("-fitting", "Fit a model into a binary image", "-fitting input-model binary-image output-model", SO_NONE);
    opts.addOption("-h", "print help message", SO_NONE);
    StringVector args = opts.ParseOptions(argc, argv, NULL);


    if (opts.GetBool("-h")) {
        cout << "## *kmesh* Usage" << endl;
        opts.PrintUsage();
        return 0;
    } else if (opts.GetBool("-smoothScalars")) {
        runScalarSmoothing(opts, args);
    } else if (opts.GetBool("-importScalars")) {
        runImportScalars(opts, args);
    } else if (opts.GetBool("-exportScalars")) {
        runExportScalars(opts, args);
    } else if (opts.GetString("-appendData", "") != "") {
        runAppendData(opts, args);
    } else if (opts.GetBool("-vti")) {
        runConvertITK2VTI(opts, args);
    } else if (opts.GetBool("-vtu")) {
        runConvertITK2VTU(opts, args);
    } else if (opts.GetBool("-traceStream")) {
        runStreamTracer(opts, args);
    } else if (opts.GetBool("-filterStream")) {
        runFilterStream(opts, args);
    } else if (opts.GetBool("-fitting")) {
        runFittingModel(opts, args);
    }
    return 0;
}
