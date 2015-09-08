//
//  vtkUtils.cpp
//  ktools
//
//  Created by Joowhi Lee on 8/17/15.
//
//

#include "vtkUtils.h"

#include <vtkNew.h>
#include <vtkMath.h>

#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkCellTreeLocator.h>
#include <vtkCellDerivatives.h>
#include <vtkGradientFilter.h>
#include <vtkVectorNorm.h>

#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPolyData.h>
#include <vtkThresholdPoints.h>
#include <vtkCleanPolyData.h>
#include <vtkModifiedBSPTree.h>


#include <vector>
#include <unordered_set>
#include <unordered_map>
//#include <vtkStreamTracer.h>
#include "kstreamtracer.h"
#include "kgeometry.h"
#include "vtkio.h"

using namespace std;
using namespace pi;

static vtkIO vio;

static bool endswith(std::string str, std::string substr) {
	size_t i = str.rfind(substr);
	return (i != string::npos) && (i == (str.length() - substr.length()));
}


// using the scalar values, find edges in the grid that intersects the given surface
void selectSurfaceIntersectingEdges(vtkDataSet* grid, string scalarName, vtkPolyData* surf) {
    
}


void sampleSurfaceScalar(vtkDataSet* grid, string markerScalar, vtkPolyData* surf, string surfScalarName) {
    
    if (surf == NULL || grid == NULL) {
        cout << "sampleSurfaceScalar(): An input is NULL: " << grid << ", " << surf << endl;
        return;
    }
    vtkNew<vtkCellLocator> surfCellLoc;
    surfCellLoc->SetDataSet(surf);
    surfCellLoc->AutomaticOn();
    surfCellLoc->BuildLocator();

    vtkDataArray* markers = grid->GetPointData()->GetArray(markerScalar.c_str());
    vtkDataArray* surfScalars = surf->GetPointData()->GetArray(surfScalarName.c_str());
    unordered_set<int> cellScalars(10);
    
    const size_t nPts = markers->GetNumberOfTuples();
    
    vtkIntArray* outputScalar = vtkIntArray::New();
    outputScalar->SetNumberOfComponents(1);
    outputScalar->SetNumberOfValues(nPts);
    outputScalar->SetName("SampledSurfaceScalars");
    
    vtkIntArray* surfCheckedCells = vtkIntArray::New();
    surfCheckedCells->SetNumberOfComponents(1);
    surfCheckedCells->SetNumberOfValues(surf->GetNumberOfCells());
    surfCheckedCells->SetName("SampledSurfaceCells");
    
    for (size_t j = 0; j < nPts; j++) {
        // exterior border is 3
        int pointScalar = 0;
        
        // if the grid point j is an exterior boundary point
        if (markers->GetTuple1(j) == 3) {
            double pts[3];
            grid->GetPoint(j, pts);
            
            // find the closest cell from the grid point j (pts)
            double closestPoint[3] = { 0, };
            vtkIdType closestCellId = 0;
            int subId = -1;
            double dist2 = 0;
            
            // what if the closest point is a vertex?
            surfCellLoc->FindClosestPoint(pts, closestPoint, closestCellId, subId, dist2);
            
            // report this cell was checekd
            surfCheckedCells->SetValue(closestCellId, 1);
            
            
            // find the closest cell of the surface from the grid point j
            vtkCell* surfCell = surf->GetCell(closestCellId);
            if (surfCell == NULL) {
                cout << "Can't find closest cell for: " << j << endl;
            } else {
                // inspect the point scalar values of the closest cell
                const size_t nCellPts = surfCell->GetNumberOfPoints();
                int scalars[4] = {0,};
                for (size_t k = 0; k < nCellPts; k++) {
                    vtkIdType cellPointId = surfCell->GetPointId(k);
                    int surfScalar = (int) surfScalars->GetTuple1(cellPointId);
                    scalars[surfScalar]++;
                }
                for (size_t k = 0; k < 4; k++) {
                    if (scalars[k] >= 2) {
                        pointScalar = k ;
                        break;
                    }
                }
            }
            outputScalar->SetValue(j, pointScalar);
        } else if (markers->GetTuple1(j) == 2 || markers->GetTuple1(j) == 1) {
            outputScalar->SetValue(j, -1);
        }
    }
    
    grid->GetPointData()->AddArray(outputScalar);
    surf->GetCellData()->AddArray(surfCheckedCells);
    
}


void findNeighborPoints(vtkDataSet* data, vtkIdType ptId, unordered_set<vtkIdType>& nbrs, vtkIdList* cellIds, vtkIdTypeArray* neighbors = NULL) {
	nbrs.clear();
    cellIds->Reset();
    data->GetPointCells(ptId, cellIds);
    for (size_t k = 0; k < cellIds->GetNumberOfIds(); k++) {
        vtkCell* cell = data->GetCell(cellIds->GetId(k));
        for (size_t l = 0; l < cell->GetNumberOfEdges(); l++) {
            vtkCell* edge = cell->GetEdge(l);
            vtkIdType s = edge->GetPointId(0);
            vtkIdType e = edge->GetPointId(1);
			vtkIdType n = -1;
            if (s == ptId) {
				n = e;
            } else if (e == ptId) {
				n = s;
            }
			if (n > -1) {
				if (nbrs.find(n) == nbrs.end()) {
					nbrs.insert(n);
					if (neighbors != NULL) {
						neighbors->SetComponent(ptId, nbrs.size() - 1, n);
					}
				}
			}
        }
    }
}

void buildAdjacencyList(vtkDataSet* surf, string roiLabel, vector<pair<vtkIdType,vector<vtkIdType> > >& graph) {
    
    // the roiArray extract the graph nodes within the ROI
    // Currently, I assume the roi is marked as -1 as I did in sampleSurfaceScalar
    vtkDataArray* roiArray = surf->GetPointData()->GetArray(roiLabel.c_str());
	if (roiArray == NULL) {
		cout << "Can't find an ROI array: " << roiLabel << endl;
		return;
	}

	const size_t nPts = roiArray->GetNumberOfTuples();
	
    graph.reserve(roiArray->GetNumberOfTuples());
    graph.clear();
	
	vtkIdTypeArray* neighborIds = vtkIdTypeArray::New();
	neighborIds->SetName("Neighbor6");
	neighborIds->SetNumberOfComponents(6);
	neighborIds->SetNumberOfTuples(nPts);
	
    unordered_set<vtkIdType> nbrs;
    vtkNew<vtkIdList> cellIds;
    for (size_t j = 0; j < nPts; j++) {
        if (roiArray->GetTuple1(j) == -1) {
            //  find the neighborhood of j
            findNeighborPoints(surf, j, nbrs, cellIds.GetPointer(), neighborIds);
            graph.push_back(make_pair(j, vector<vtkIdType>(nbrs.begin(), nbrs.end())));
        }
        if (j % 1000 == 1) {
            cout << "Processed " << (j+1) << " points ..." << endl;
        }
    }
    
	surf->GetPointData()->AddArray(neighborIds);
}









/// @brief Perform a line clipping task
void runTraceClipping(Options& opts, StringVector& args) {
	string inputStreamsFile = args[0];
	string inputObjectFile = args[1];
	string outputStreamsFile = args[2];
	
	vtkIO vio;
	vtkPolyData* inputStream = vio.readFile(inputStreamsFile);
	vtkPolyData* inputObject = vio.readFile(inputObjectFile);
	vtkPolyData* outputObject = vtkPolyData::New();
	
	
	vtkCellArray* lines = inputStream->GetLines();
	vtkModifiedBSPTree* tree = vtkModifiedBSPTree::New();
	tree->SetDataSet(inputObject);
	tree->BuildLocator();
	
	vtkPoints* outputPoints = vtkPoints::New();
	vtkCellArray* outputLines = vtkCellArray::New();
	
	for (int i = 0; i < lines->GetNumberOfCells(); i++) {
		vtkCell* line = inputStream->GetCell(i);
		double length = 0;
//		performLineClipping(inputStream, tree, i, line, inputObject, outputPoints, outputLines, length);
	}
	vio.writeFile("test.vtp", outputObject);
}






/// @brief Execute the stream tracer
//void runStreamTracer(Options& opts, StringVector& args) {
////    string inputVTUFile = args[0];
////    string inputSeedPointsFile = args[1];
////    string outputStreamFile = args[2];
////    string outputPointFile = args[3];
////    bool zRotate = opts.GetBool("-zrotate", false);
////    
////    vtkIO vio;
////    vtkDataSet* inputData = vio.readDataFile(inputVTUFile);
////    vtkPolyData* inputSeedPoints = vio.readFile(inputSeedPointsFile);
////	
////	vtkNew<vtkThresholdPoints> selector;
////	selector->SetInput(inputData);
////	selector->ThresholdBetween(1.5, 2.5);
////	selector->SetInputArrayToProcess(0, 0, 0, vtkDataSet::FIELD_ASSOCIATION_POINTS, "SampledSurfaceScalars");
////	selector->Update();
////	vtkPolyData* selectedSeeds = selector->GetOutput();
////	vio.writeFile("selectedSeeds.vtp", selectedSeeds);
////	
//////	selectedSeeds->Print(cout);
////	
//////    vtkPolyData* outputStream = performStreamTracer(opts, inputData, selectedSeeds, zRotate);
////	
////    cout << selectedSeeds->GetPointData()->GetArray("LineOK")->GetNumberOfTuples() << endl;
////    cout << selectedSeeds->GetPointData()->GetArray("Length")->GetNumberOfTuples() << endl;
////	
////	
////    vio.writeFile(outputPointFile, selectedSeeds);
////	if (outputStream) {
////	    vio.writeFile(outputStreamFile, outputStream);
////	}
//}


/// @brief Copy a scalar list from a seed object to a stream line object
void runTraceScalarCombine(Options& opts, StringVector& args) {
	if (args.size() < 3) {
		cout << "requires input-seed input-stream output-stream-file" << endl;
		return;
	}
	
	string inputSeedFile = args[0];
	string inputStreamFile = args[1];
	string outputStreamFile = args[2];
	string scalarName = opts.GetString("-scalarName");
	
	if (scalarName == "") {
		cout << "requires -scalarName scalarName" << endl;
		return;
	}
	
	vtkIO vio;
	vtkPolyData* inputSeed = vio.readFile(inputSeedFile);
	vtkPolyData* inputStream = vio.readFile(inputStreamFile);
	
	vtkDataArray* pointIds = inputStream->GetCellData()->GetScalars("PointIds");
	if (pointIds == NULL) {
		cout << "Can't find PointIds" << endl;
		return;
	}
	vtkDataArray* scalars = inputSeed->GetPointData()->GetScalars(scalarName.c_str());
	if (scalars == NULL) {
		cout << "Can't find scalars: " << scalarName << endl;
		return;
	}
	
	vtkDoubleArray* outputScalars = vtkDoubleArray::New();
	outputScalars->SetName(scalarName.c_str());
	for (int i = 0; i < pointIds->GetNumberOfTuples(); i++) {
		int ptId = pointIds->GetTuple1(i);
		double value = scalars->GetTuple1(ptId);
		outputScalars->InsertNextTuple1(value);
	}
	inputStream->GetCellData()->AddArray(outputScalars);
	
	if (opts.GetBool("-zrotate")) {
		cout << "The output is rotated!" << endl;
		vio.zrotate(inputStream);
	}
	vio.writeFile(outputStreamFile, inputStream);
}


/// @brief Apply a filter to each stream line
void runStreamLineThreshold(Options& opts, StringVector& args) {
	string inputStream = args[0];
	string inputSeeds = args[1];
	string outputStreamFile = args[2];
	string scalarName = opts.GetString("-scalarName");
	
	double lowThreshold = opts.GetStringAsReal("-thresholdMin", DBL_MIN);
	double highThreshold = opts.GetStringAsReal("-thresholdMax", DBL_MAX);
	
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

/// @brief Rescale the streamline with a given length
void runRescaleStream(Options& opts, StringVector& args) {

}





void runExtractSurfaceBorder(Options& opts, StringVector& args) {
    string inputFile = args[0];
    int scalarValue = atoi(args[1].c_str());
    string outputFile = args[2];
    
    vtkPolyData* ds = vio.readFile(inputFile);
    
    for (size_t j = 0; j < ds->GetNumberOfCells(); j++) {
        
    }
    
    
}


//
void processVTKUtilsOptions(pi::Options& opts) {
    opts.addOption("-sampleSurfaceScalars", "For each point marked as 3, sample the closest cell's majority scalar value.", "-sampleSurfaceScalars input-dataset input-surface output-dataset (-o output-surface)", SO_NONE);
    opts.addOption("-buildAdjacencyList", "Build an adjacency list for the input data and its point scalars", "-buildAdjacencyList input-data scalar-name", SO_NONE);
    
    opts.addOption("-extractSurfaceBorder", "Extract a borderline of a surface patch", "-extractSurfaceBorder input-surface scalar-value out-file -scalarName scalarName ", SO_NONE);
	
	opts.addOption("-traceScalarCombine", "Combine scalar values from a seed object to a stream line object. The stream line object must have PointIds for association. -zrotate option will produce the rotated output.", "-traceScalarCombine stream_seed.vtp stream_lines.vtp stream_lines_output.vtp -scalarName scalarToBeCopied", SO_NONE);
	opts.addOption("-rescaleStream", "Rescale streamlines to fit with given lengths", "-rescaleStream input-stream-lines length.txt or input.vtp -scalarName scalarname", SO_NONE);
	opts.addOption("-traceStream", "Trace a stream line from a given point set", "-traceStream input-vtu-field input-vtk output-lines output-points", SO_NONE);
	opts.addOption("-traceDirection", "Choose the direction of stream tracing (both, forward, backward)", "-traceStream ... -traceDirection (both|forward|backward)", SO_REQ_SEP);
	opts.addOption("-traceClipping", "Clip stream lines to fit with an object", "-traceClipping stream_lines.vtp stream_object.vtp stream_lines_output.vtp", SO_NONE);
	opts.addOption("-thresholdStream", "Remove stream lines which are lower than a given threshold", "-thresholdStream stream-line-input stream-seed-input stream-line-output -scalarName scalar -threshold xx", SO_NONE);

}


void processVTKUtils(pi::Options opts, pi::StringVector args) {
	vtkIO vio;
	string input1File, input2File, outputFile;
	if (opts.GetBool("-sampleSurfaceScalars")) {
        input1File = args[0];
        input2File = args[1];
        outputFile = args[2];
        vtkDataSet* data1 = vio.readDataFile(input1File);
        vtkPolyData* surf1 = vio.readFile(input2File);
        sampleSurfaceScalar(data1, opts.GetString("-scalarName", "BorderPoints").c_str(), surf1, "labels");
        vio.writeFile(outputFile, data1);
        if (opts.GetString("-o") != "") {
            vio.writeFile(opts.GetString("-o"), surf1);
        }
    } else if (opts.GetBool("-buildAdjacencyList")) {
        input1File = args[0];
		outputFile = args[1];
		
        string scalarName = opts.GetString("-scalarName", "SampledSurfaceScalars");
        vtkDataSet* data1 = vio.readDataFile(input1File);
        
        vector<pair<vtkIdType, vector<vtkIdType> > > graph;
        buildAdjacencyList(data1, scalarName, graph);
		
		vio.writeFile(outputFile, data1);
	} else  if (opts.GetBool("-thresholdStream")) {
		runStreamLineThreshold(opts, args);
	} else if (opts.GetBool("-rescaleStream")) {
		runRescaleStream(opts, args);
	} else if (opts.GetBool("-traceClipping")) {
		runTraceClipping(opts, args);
	} else if (opts.GetBool("-traceScalarCombine")) {
		runTraceScalarCombine(opts, args);
    } else if (opts.GetBool("-extractSurfaceBorder")) {
        runExtractSurfaceBorder(opts, args);
    }
}
