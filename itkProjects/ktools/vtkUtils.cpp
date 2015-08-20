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
#include "vtkio.h"

using namespace std;
using namespace pi;


static bool endswith(std::string str, std::string substr) {
	size_t i = str.rfind(substr);
	return (i != string::npos) && (i == (str.length() - substr.length()));
}

vtkIntArray* selectBoundaryCells(vtkDataSet* data, string scalarName) {
    const size_t nCells = data->GetNumberOfCells();
	vtkDataArray* array = data->GetPointData()->GetArray(scalarName.c_str());
	vector<int> scalarValue;
	vtkIntArray* intArray = vtkIntArray::New();
	intArray->SetNumberOfComponents(1);
	intArray->SetNumberOfTuples(nCells);
	
    for (size_t j = 0; j < nCells; j++) {
        vtkCell* cell = data->GetCell(j);
		scalarValue.clear();
        for (size_t k = 0; k < cell->GetNumberOfPoints(); k++) {
			vtkIdType cellPtId = cell->GetPointId(k);
			scalarValue.push_back(array->GetTuple1(cellPtId));
        }
		
		bool borderCell = false;
		if (scalarValue.size() > 0) {
			int value = scalarValue[0];
			for (size_t j = 1; j < scalarValue.size(); j++) {
				if (scalarValue[j] != value) {
					borderCell = true;
				}
			}
		}
		
		if (borderCell) {
			intArray->SetValue(j, 1);
		}
    }
	
	intArray->SetName("BorderCells");
	data->GetCellData()->AddArray(intArray);
	
	return intArray;
}


// mark the boundary points (2: interior boundary, 3: exterior boundary)
vtkDataArray* selectBoundaryPoints(vtkDataSet* data, std::string scalarName) {
    
    vtkDataArray* interiorMaker = data->GetPointData()->GetArray(scalarName.c_str());
    if (interiorMaker == NULL) {
        cout << "Can't find scalar values: " << scalarName << endl;
        return NULL;
    }
    
    const size_t npts = data->GetNumberOfPoints();
    
    vtkIntArray* boundaryMarker = vtkIntArray::New();
    boundaryMarker->SetNumberOfComponents(1);
    boundaryMarker->SetNumberOfTuples(npts);
	boundaryMarker->FillComponent(0, 0);
	

	unordered_set<vtkIdType> nbrs;
    vtkNew<vtkIdList> cellIds;

    for (size_t j = 0; j < npts; j++) {
        if (interiorMaker->GetTuple1(j) != 1) {
            continue;
        }
        
        // iterate over neighbor cells and find neighbor points
		nbrs.clear();
        cellIds->Reset();
        data->GetPointCells(j, cellIds.GetPointer());
        for (size_t k = 0; k < cellIds->GetNumberOfIds(); k++) {
            vtkCell* cell = data->GetCell(cellIds->GetId(k));
            for (size_t l = 0; l < cell->GetNumberOfEdges(); l++) {
                vtkCell* edge = cell->GetEdge(l);
                vtkIdType s = edge->GetPointId(0);
                vtkIdType e = edge->GetPointId(1);
                if (s == j) {
					nbrs.insert(e);
                } else if (e == j) {
					nbrs.insert(s);
                }
            }
        }
        
        // check neighbor points and find exterior points
        bool surfacePoint = false;
		unordered_set<vtkIdType>::const_iterator iter = nbrs.begin();
		for (; iter != nbrs.end(); iter++) {
			// if the neighbor is outside point
			vtkIdType nbrId = *iter;
            if (interiorMaker->GetTuple1(nbrId) == 0) {
                boundaryMarker->SetTuple1(nbrId, 3);
				// then j is a surface point
                surfacePoint = true;
            }
        }
        if (surfacePoint) {
            boundaryMarker->SetTuple1(j, 2);
        } else {
            boundaryMarker->SetTuple1(j, 1);
        }
    }
	
    boundaryMarker->SetName("BorderPoints");
    data->GetPointData()->AddArray(boundaryMarker);
    return boundaryMarker;
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




// Compute Laplace PDE based on the adjacency list and border
void computeLaplacePDE(vtkDataSet* data, string selectedPointsName, string neighborName, string boundaryCondName, const double low, const double high, const int nIters, const double dt) {
	
	if (data == NULL) {
		cout << "Data input is NULL" << endl;
		return;
	}
	
	// check interior points
	vtkDataArray* selectedPoints = data->GetPointData()->GetArray(selectedPointsName.c_str());
	if (selectedPoints == NULL) {
		cout << "No scalar values for SelectedPoints" << endl;
		return;
	}
	
	// check boundary points
	vtkDataArray* boundaryPoints = data->GetPointData()->GetArray(boundaryCondName.c_str());
	if (boundaryPoints == NULL) {
		cout << "No scalar values for BoundaryPoints" << endl;
		return;
	}
	
	// check boundary points
	vtkDataArray* neighborIds = data->GetPointData()->GetArray(neighborName.c_str());
	if (neighborIds == NULL) {
		cout << "No scalar values for Neighbor IDs" << endl;
		return;
	}
	

	struct LaplaceGrid {
		double low;
		double high;
		double dt;
		
		vector<vtkIdType> solutionDomain;
		vtkDataArray* neighborIds;
		vtkIntArray* boundaryCond;
		
		vtkDoubleArray* solution;
		vtkDoubleArray* tmpSolution;
		
		LaplaceGrid(double l, double h, double d, vtkDataArray* sp, vtkDataArray* nids, vtkDataArray* bcond): low(l), high(h), dt(d) {
			
			extractSolutionDomain(sp);
			neighborIds = nids;
			boundaryCond = vtkIntArray::SafeDownCast(bcond);
			
			initializeSolution();
		}
		
		void extractSolutionDomain(vtkDataArray* sp) {
			for (size_t j = 0; j < sp->GetNumberOfTuples(); j++) {
				const int domain = sp->GetTuple1(j);
				if (domain == 1 || domain == 2) {
					// keep an interior point id
					solutionDomain.push_back(j);
				}
			}
			cout << "# of solution grids: " << solutionDomain.size() << endl;
		}
		
		void initializeSolution() {
			solution = vtkDoubleArray::New();
			solution->SetName("LaplacianSolution");
			solution->SetNumberOfComponents(1);
			solution->SetNumberOfTuples(boundaryCond->GetNumberOfTuples());
			solution->FillComponent(0, 0);

			tmpSolution = vtkDoubleArray::New();
			tmpSolution->SetNumberOfComponents(1);
			tmpSolution->SetNumberOfTuples(boundaryCond->GetNumberOfTuples());

			const size_t nPts = boundaryCond->GetNumberOfTuples();
			for (size_t j = 0; j < nPts; j++) {
				if (boundaryCond->GetValue(j) == 1) {
					// high
					solution->SetValue(j, high);
				} else if (boundaryCond->GetValue(j) == 2){
					// low
					solution->SetValue(j, low);
				} else {
					solution->SetValue(j, 0);
				}
				tmpSolution->SetValue(j, solution->GetValue(j));
			}
		}
		
		void computeStep() {
			const size_t nPts = solutionDomain.size();
			for (size_t j = 0; j < nPts; j++) {
				vtkIdType centerId = solutionDomain[j];
				const double* nbrs = neighborIds->GetTuple(centerId);
				double u = 0;
				for (size_t k = 0; k < 6; k++) {
					const vtkIdType kId = nbrs[k];
					u += solution->GetValue(kId);
				}
				u = u/6.0;
				tmpSolution->SetValue(centerId, u);
			}
			solution->DeepCopy(tmpSolution);
		}
		
		
		void computeNormals(vtkDataSet* data) {
			/*
			vtkNew<vtkCellDerivatives> deriv;
			deriv->SetInput(data);
			deriv->SetVectorModeToComputeGradient();
			deriv->Update();
			vtkDataSet* derivOut = deriv->GetOutput();
			derivOut->GetCellData()->SetActiveVectors("ScalarGradient");
			vtkDataArray* scalarGradient = deriv->GetOutput()->GetCellData()->GetArray("ScalarGradient");
			scalarGradient->SetName("LaplacianGradient");
			*/
			
			vtkNew<vtkGradientFilter> gradFilter;
			gradFilter->SetInput(data);
			gradFilter->SetInputScalars(vtkDataSet::FIELD_ASSOCIATION_POINTS, "LaplacianSolution");
			gradFilter->SetResultArrayName("LaplacianGradient");
			gradFilter->Update();
			vtkDataArray* scalarGradient = gradFilter->GetOutput()->GetPointData()->GetArray("LaplacianGradient");
			
			vtkDoubleArray* norms = vtkDoubleArray::New();
			norms->SetName("LaplacianGradientNorm");
			norms->SetNumberOfComponents(3);
			norms->SetNumberOfTuples(scalarGradient->GetNumberOfTuples());
			
			const size_t nPts = norms->GetNumberOfTuples();
			for (size_t j = 0; j < nPts; j++) {
				double* vec = scalarGradient->GetTuple3(j);
				vtkMath::Normalize(vec);
				norms->SetTuple3(j, vec[0], vec[1], vec[2]);
			}
			
			data->GetPointData()->AddArray(scalarGradient);
			data->GetPointData()->SetVectors(norms);
		}
	};
	
	
	LaplaceGrid grid(low, high, dt, selectedPoints, neighborIds,boundaryPoints);

	// main iteration loop
	for (size_t i = 1; i <= nIters; i++) {
		cout << "iteration: " << i << endl;
		grid.computeStep();
	}
	
	
	// return the solution
	data->GetPointData()->AddArray(grid.solution);
	grid.computeNormals(data);
}



/// @brief perform a line clipping to fit within the object
bool performLineClipping(vtkPolyData* streamLines, vtkModifiedBSPTree* tree, int lineId, vtkCell* lineToClip, vtkPolyData* object, vtkPoints* outputPoints, vtkCellArray* outputLines, double &length) {

	/// - Iterate over all points in a line
	vtkIdList* ids = lineToClip->GetPointIds();
	/// - Identify a line segment included in the line
	
	int nIntersections = 0;
	bool foundEndpoint = false;
	std::vector<vtkIdType> idList;
	for (int j = 2; j < ids->GetNumberOfIds(); j++) {
		double p1[3], p2[3];
		streamLines->GetPoint(ids->GetId(j-1), p1);
		streamLines->GetPoint(ids->GetId(j), p2);
		
		// handle initial condition
		if (j == 2) {
			double p0[3];
			streamLines->GetPoint(ids->GetId(0), p0);
			idList.push_back(outputPoints->GetNumberOfPoints());
			outputPoints->InsertNextPoint(p0);
			
			idList.push_back(outputPoints->GetNumberOfPoints());
			outputPoints->InsertNextPoint(p1);
			
			length = sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
		}
		
		int subId;
		double x[3] = {-1,-1,-1};
		double t = 0;
		
		double pcoords[3] = { -1, };
		int testLine = tree->IntersectWithLine(p1, p2, 0.01, t, x, pcoords, subId);
		if (testLine) {
			nIntersections ++;
			if (nIntersections > 0) {
				idList.push_back(outputPoints->GetNumberOfPoints());
				outputPoints->InsertNextPoint(x);
				length += sqrt(vtkMath::Distance2BetweenPoints(p1, x));
				foundEndpoint = true;
				break;
			}
		}
		//        cout << testLine << "; " << x[0] << "," << x[1] << "," << x[2] << endl;
		
		
		idList.push_back(outputPoints->GetNumberOfPoints());
		outputPoints->InsertNextPoint(p2);
		length += sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
	}
	
	if (foundEndpoint) {
		outputLines->InsertNextCell(idList.size(), &idList[0]);
		return true;
	}
	return false;
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
		performLineClipping(inputStream, tree, i, line, inputObject, outputPoints, outputLines, length);
	}
	vio.writeFile("test.vtp", outputObject);
}





vtkPolyData* performStreamTracer(Options& opts, vtkDataSet* inputData, vtkPolyData* inputSeedPoints, bool zRotate = false) {

	// set active velocity field
	inputData->GetPointData()->SetActiveVectors("LaplacianGradientNorm");

	/// - Converting the input points to the image coordinate
	vtkPoints* points = inputSeedPoints->GetPoints();
	const int nInputPoints = inputSeedPoints->GetNumberOfPoints();
	if (zRotate) {
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
		inputSeedPoints->SetPoints(points);
	}

	inputSeedPoints->Print(cout);
	
	/// StreamTracer should have a point-wise gradient field
	/// - Set up tracer (Use RK45, both direction, initial step 0.05, maximum propagation 500
	StreamTracer* tracer = StreamTracer::New();
	tracer->SetInput(inputData);
	tracer->SetSource(inputSeedPoints);
	tracer->SetComputeVorticity(false);
	
	//    double seedPoint[3];
 //   inputPoints->GetPoint(24745, seedPoint);
 //   tracer->SetStartPosition(seedPoint);
	tracer->SetIntegratorTypeToRungeKutta45();

	
	bool isBothDirection = false;
	if (opts.GetString("-traceDirection") == "both") {
		tracer->SetIntegrationDirectionToBoth();
		isBothDirection = true;
		cout << "Forward/Backward Integration" << endl;
	} else if (opts.GetString("-traceDirection") == "backward") {
		tracer->SetIntegrationDirectionToBackward();
		cout << "Backward Integration" << endl;
	} else {
		tracer->SetIntegrationDirectionToForward();
		cout << "Forward Integration" << endl;
	}
	
	tracer->SetInterpolatorTypeToDataSetPointLocator();
//	tracer->SetInterpolatorTypeToCellLocator();
	tracer->SetMaximumPropagation(500);
	tracer->SetInitialIntegrationStep(0.05);
	tracer->Update();
	
	
	vtkPolyData* streamLines = tracer->GetOutput();
	streamLines->Print(cout);

	// remove useless pointdata information
	streamLines->GetPointData()->Reset();
	
	
	// loop over the cell and compute the length
	int nCells = streamLines->GetNumberOfCells();
	cout << "# of cells: " << nCells << endl;

	
	vtkIO vio;
	vio.writeFile("stream.vtp", streamLines);
	
	/// - Prepare the output as a scalar array
	//    vtkDataArray* streamLineLength = streamLines->GetCellData()->GetScalars("Length");
	
	/// - Prepare the output for the input points
	vtkDoubleArray* streamLineLengthPerPoint = vtkDoubleArray::New();
	streamLineLengthPerPoint->SetNumberOfTuples(nInputPoints);
	streamLineLengthPerPoint->SetName("Length");
	streamLineLengthPerPoint->SetNumberOfComponents(1);
	streamLineLengthPerPoint->FillComponent(0, 0);
	
	vtkIntArray* lineCorrect = vtkIntArray::New();
	lineCorrect->SetName("LineOK");
	lineCorrect->SetNumberOfValues(nInputPoints);
	lineCorrect->FillComponent(0, 0);
	
	inputSeedPoints->GetPointData()->SetScalars(streamLineLengthPerPoint);
	inputSeedPoints->GetPointData()->AddArray(lineCorrect);
	
	cout << "Assigning a length to each source vertex ..." << endl;
	vtkDataArray* seedIds = streamLines->GetCellData()->GetScalars("SeedId");
	if (seedIds) {
		// line clipping
		vtkPoints* outputPoints = vtkPoints::New();
		vtkCellArray* outputCells = vtkCellArray::New();
		
		/// construct a tree locator
		vtkModifiedBSPTree* tree = vtkModifiedBSPTree::New();
		tree->SetDataSet(inputSeedPoints);
		tree->BuildLocator();
		
		
		vtkDoubleArray* lengthArray = vtkDoubleArray::New();
		lengthArray->SetName("Length");
		
		vtkIntArray* pointIds = vtkIntArray::New();
		pointIds->SetName("PointIds");
		
		for (int i = 0; i < nCells; i++) {
			int pid = seedIds->GetTuple1(i);
			double length = 0;
			if (pid > -1) {
				vtkCell* line = streamLines->GetCell(i);
				/// - Assume that a line starts from a point on the input mesh and must meet at the opposite surface of the starting point.
				bool lineAdded = performLineClipping(streamLines, tree, i, line, inputSeedPoints, outputPoints, outputCells, length);
				
				if (lineAdded) {
					pointIds->InsertNextValue(pid);
					lengthArray->InsertNextValue(length);
					streamLineLengthPerPoint->SetValue(pid, length);
					lineCorrect->SetValue(pid, 1);
				} else {
					lineCorrect->SetValue(pid, 2);
				}
			}
		}
		
		vtkPolyData* outputStreamLines = vtkPolyData::New();
		outputStreamLines->SetPoints(outputPoints);
		outputStreamLines->SetLines(outputCells);
		outputStreamLines->GetCellData()->AddArray(pointIds);
		outputStreamLines->GetCellData()->AddArray(lengthArray);
		
		
		vtkCleanPolyData* cleaner = vtkCleanPolyData::New();
		cleaner->SetInput(outputStreamLines);
		cleaner->Update();
        return cleaner->GetOutput();
	} else {
		cout << "Can't find SeedId" << endl;
	}
	
	//    vio.writeXMLFile(outputVTKFile, streamLines);
    return NULL;
}


/// @brief Execute the stream tracer
void runStreamTracer(Options& opts, StringVector& args) {
    string inputVTUFile = args[0];
    string inputSeedPointsFile = args[1];
    string outputStreamFile = args[2];
    string outputPointFile = args[3];
    bool zRotate = opts.GetBool("-zrotate", false);
    
    vtkIO vio;
    vtkDataSet* inputData = vio.readDataFile(inputVTUFile);
    vtkPolyData* inputSeedPoints = vio.readFile(inputSeedPointsFile);
	
	vtkNew<vtkThresholdPoints> selector;
	selector->SetInput(inputData);
	selector->ThresholdBetween(1.5, 2.5);
	selector->SetInputArrayToProcess(0, 0, 0, vtkDataSet::FIELD_ASSOCIATION_POINTS, "SampledSurfaceScalars");
	selector->Update();
	vtkPolyData* selectedSeeds = selector->GetOutput();
	vio.writeFile("selectedSeeds.vtp", selectedSeeds);
	
//	selectedSeeds->Print(cout);
	
    vtkPolyData* outputStream = performStreamTracer(opts, inputData, selectedSeeds, zRotate);
    
    cout << selectedSeeds->GetPointData()->GetArray("LineOK")->GetNumberOfTuples() << endl;
    cout << selectedSeeds->GetPointData()->GetArray("Length")->GetNumberOfTuples() << endl;
	
	
    vio.writeFile(outputPointFile, selectedSeeds);
	if (outputStream) {
	    vio.writeFile(outputStreamFile, outputStream);
	}
}


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


void runMeasureThickness(Options& opts, StringVector& args) {
    vtkIO vio;
    string inputFile = args[0];
    string inputSeedFile = args[1];
    string outputStreamFile = args[2];
    
    vtkDataSet* data = vio.readDataFile(inputFile);
    vtkPolyData* inputSeedPoints = vio.readFile(inputSeedFile);

	vtkNew<vtkThresholdPoints> selector;
	selector->SetInput(inputSeedPoints);
	selector->ThresholdBetween(0.5, 1.5);
	selector->SetInputArrayToProcess(0, 0, 0, vtkDataSet::FIELD_ASSOCIATION_POINTS, "SampledSurfaceScalars");
	selector->Update();
	vtkPolyData* selectedSeeds = selector->GetOutput();

	
    computeLaplacePDE(data, "SelectedPoints", "Neighbor6", "SampledSurfaceScalars", 0, 10000, 2000, 0.065);
    
    vtkPolyData* outputStream = performStreamTracer(opts, data, selectedSeeds);
    vio.writeFile(outputStreamFile, outputStream);
   
}


//
void processVTKUtilsOptions(pi::Options& opts) {
    opts.addOption("-markBorderCells", "Mark border cells of an input dataset. The border cells have 1 in BorderCells data", "-markBorderCells input-data output-data", SO_NONE);
    opts.addOption("-markBorderPoints", "Mark border points of an input dataset. The border points will be marked as 2 and its exterior neighbors will be marked as 3.", "-markBorderPoints input-data output-data", SO_NONE);
    opts.addOption("-sampleSurfaceScalars", "For each point marked as 3, sample the closest cell's majority scalar value.", "-sampleSurfaceScalars input-dataset input-surface output-dataset (-o output-surface)", SO_NONE);
    opts.addOption("-buildAdjacencyList", "Build an adjacency list for the input data and its point scalars", "-buildAdjacencyList input-data scalar-name", SO_NONE);
	opts.addOption("-computeLaplacePDE", "Compute the Laplace PDE over the given domain", "-computeLaplacePDE input-data output-data selectedPointsScalarName, neighborIdScalarName, boundaryConditionScalarName", SO_NONE);
	
	opts.addOption("-traceScalarCombine", "Combine scalar values from a seed object to a stream line object. The stream line object must have PointIds for association. -zrotate option will produce the rotated output.", "-traceScalarCombine stream_seed.vtp stream_lines.vtp stream_lines_output.vtp -scalarName scalarToBeCopied", SO_NONE);
	opts.addOption("-rescaleStream", "Rescale streamlines to fit with given lengths", "-rescaleStream input-stream-lines length.txt or input.vtp -scalarName scalarname", SO_NONE);
	opts.addOption("-traceStream", "Trace a stream line from a given point set", "-traceStream input-vtu-field input-vtk output-lines output-points", SO_NONE);
	opts.addOption("-traceDirection", "Choose the direction of stream tracing (both, forward, backward)", "-traceStream ... -traceDirection (both|forward|backward)", SO_REQ_SEP);
	opts.addOption("-traceClipping", "Clip stream lines to fit with an object", "-traceClipping stream_lines.vtp stream_object.vtp stream_lines_output.vtp", SO_NONE);
	opts.addOption("-thresholdStream", "Remove stream lines which are lower than a given threshold", "-thresholdStream stream-line-input stream-seed-input stream-line-output -scalarName scalar -threshold xx", SO_NONE);

    opts.addOption("-measureThickness", "Measure the thickness of the solution domain via RK45 integration", "-measureThickness input-data", SO_NONE);
}


void processVTKUtils(pi::Options opts, pi::StringVector args) {
	vtkIO vio;
	string input1File, input2File, outputFile;
	if (opts.GetBool("-markBorderCells")) {
		input1File = args[0];
		outputFile = args[1];
		vtkDataSet* data1 = vio.readDataFile(input1File);
		string scalarName = opts.GetString("-scalarName", "SelectedPoints");
		selectBoundaryCells(data1, scalarName);
		vio.writeFile(outputFile, data1);
    } else if (opts.GetBool("-markBorderPoints")) {
		input1File = args[0];
		outputFile = args[1];
        vtkDataSet* data1 = vio.readDataFile(input1File);
        string scalarName = opts.GetString("-scalarName", "SelectedPoints");
        selectBoundaryPoints(data1, scalarName);
        vio.writeFile(outputFile, data1);
    } else if (opts.GetBool("-sampleSurfaceScalars")) {
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
	} else if (opts.GetBool("-computeLaplacePDE")) {
		input1File = args[0];
		outputFile = args[1];

		string selectedPointName = args[2];
		string neighborIdName = args[3];
		string boundaryCondName = args[4];
		
		vtkDataSet* data = vio.readDataFile(input1File);
		computeLaplacePDE(data, selectedPointName, neighborIdName, boundaryCondName, 0, 10000, 2000, 0.065);
		vio.writeFile(outputFile, data);
	} else if (opts.GetBool("-traceStream")) {
		runStreamTracer(opts, args);
	} else if (opts.GetBool("-thresholdStream")) {
		runStreamLineThreshold(opts, args);
	} else if (opts.GetBool("-rescaleStream")) {
		runRescaleStream(opts, args);
	} else if (opts.GetBool("-traceClipping")) {
		runTraceClipping(opts, args);
	} else if (opts.GetBool("-traceScalarCombine")) {
		runTraceScalarCombine(opts, args);
    } else if (opts.GetBool("-measureThickness")) {
        runMeasureThickness(opts, args);
    }
}
