//
//  kvolume.cpp
//  ktools
//
//  Created by Joowhi Lee on 9/3/15.
//
//

#include "kvolume.h"
#include "kstreamtracer.h"
#include "kgeometry.h"

#include "piOptions.h"

#include "vtkio.h"

#include <vtkPolyData.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkNew.h>
#include <vtkExtractGrid.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkGradientFilter.h>
#include <vtkMath.h>
#include <vtkThresholdPoints.h>
#include <vtkCleanPolyData.h>
#include <vtkModifiedBSPTree.h>
#include <vtkCellLocator.h>
#include <vtkPointLocator.h>
#include <vtkGenericCell.h>
#include <vtkPolyDataNormals.h>
#include <vtkInterpolatedVelocityField.h>

#include <unordered_map>
#include <unordered_set>

using namespace pi;
using namespace std;

static vtkIO vio;


void runExtractBorderline(Options& opts, StringVector& args) {
    string inputFile = args[0];
    string outputFile = args[1];
    string scalarName = opts.GetString("-scalarName", "labels");
    
    vtkIO vio;
    vtkPolyData* input = vio.readFile(inputFile);
    input->BuildCells();
    input->BuildLinks();
    
    cout << input->GetNumberOfPoints() << endl;
    cout << input->GetNumberOfLines() << endl;
    
    vtkPoints* points = input->GetPoints();
    vtkDataArray* scalar = input->GetPointData()->GetArray(scalarName.c_str());
    
    vector<pair<vtkIdType,vtkIdType> > edgeSet;
    
    for (size_t j = 0; j < input->GetNumberOfCells(); j++) {
        vtkCell* cell = input->GetCell(j);
        vtkIdType p[3]; int s[3];
        p[0] = cell->GetPointId(0);
        p[1] = cell->GetPointId(1);
        p[2] = cell->GetPointId(2);
        
        s[0] = scalar->GetTuple1(p[0]);
        s[1] = scalar->GetTuple1(p[1]);
        s[2] = scalar->GetTuple1(p[2]);
        
        if (s[0] == s[1] && s[1] == s[2] && s[2] == s[0]) {
            continue;
        }
        
        vtkIdType p1, p2;
        if (s[0] != s[1] && s[0] != s[2]) {
            p1 = p[1];
            p2 = p[2];
        } else if (s[1] != s[2] && s[1] != s[0]) {
            p1 = p[2];
            p2 = p[0];
        } else if (s[2] != s[0] && s[2] != s[1]) {
            p1 = p[0];
            p2 = p[1];
        } else {
            continue;
        }
        
        edgeSet.push_back(make_pair(p1, p2));
    }
    
    vtkPolyData* output = vtkPolyData::New();
    output->SetPoints(points);
    
    vtkCellArray* lines = vtkCellArray::New();
    for (size_t j = 0; j < edgeSet.size(); j++) {
        vtkIdList* ids = vtkIdList::New();
        ids->InsertNextId(edgeSet[j].first);
        ids->InsertNextId(edgeSet[j].second);
        lines->InsertNextCell(ids->GetNumberOfIds(), ids->GetPointer(0));
        ids->Delete();
    }
    
    output->SetLines(lines);
    output->BuildCells();
    
    vio.writeFile(outputFile, output);
    cout << "Length of Borderline: " << edgeSet.size() << endl;
}

vtkDataSet* createGridForSphereLikeObject(vtkPolyData* input, int& insideCount, int dims = 100, bool insideOutOn = true) {
    // x1-x2, y1-y2, z1-z2
    double* bounds = input->GetBounds();
    
    cout << bounds[0] << "," << bounds[1] << endl;
    cout << bounds[2] << "," << bounds[3] << endl;
    cout << bounds[4] << "," << bounds[5] << endl;
    cout << "Grid Dimension: " << dims << endl;
    
    vtkStructuredGrid* grid = vtkStructuredGrid::New();
    grid->SetDimensions(dims + 6, dims + 6, dims + 6);
    
    vtkPoints* gridPoints = vtkPoints::New();
    
    for (int k = 0; k < dims + 6; k++) {
        for (int j = 0; j < dims + 6; j++) {
            for (int i = 0; i < dims + 6; i++) {
                double x = bounds[0] + (i-3)*(bounds[1]-bounds[0])/dims;
                double y = bounds[2] + (j-3)*(bounds[3]-bounds[2])/dims;
                double z = bounds[4] + (k-3)*(bounds[5]-bounds[4])/dims;
                
                gridPoints->InsertNextPoint(x, y, z);
            }
        }
    }
    
    grid->SetPoints(gridPoints);
    
    vtkSelectEnclosedPoints* encloser = vtkSelectEnclosedPoints::New();
    encloser->SetInput(grid );
    encloser->SetSurface(input);
    encloser->CheckSurfaceOn();
    if (insideOutOn) {
        cout << "Inside Out Mode" << endl;
        encloser->InsideOutOn();
    }
    encloser->SetTolerance(0);
    encloser->Update();
    
    vtkDataArray* selectedPoints = encloser->GetOutput()->GetPointData()->GetArray("SelectedPoints");
    for (size_t j = 0; j < selectedPoints->GetNumberOfTuples(); j++) {
        if (selectedPoints->GetTuple1(j) == 1) {
            insideCount++;
        }
    }
    return encloser->GetOutput();
}


// create a structured grid with the size of input
// convert the grid to polydata
// create the intersection between the grid and the polydata
void runFillGrid(Options& opts, StringVector& args) {
    if (opts.GetBool("-twosided")) {
        string inputFileOut = args[0];
        string inputFileIn = args[1];
        string outputFile = args[2];
        
        vtkIO vio;
        vtkPolyData* inputOut = vio.readFile(inputFileOut);
        vtkPolyData* inputIn = vio.readFile(inputFileIn);
        
        // x1-x2, y1-y2, z1-z2
        double* bounds = inputOut->GetBounds();
        
        cout << bounds[0] << "," << bounds[1] << endl;
        cout << bounds[2] << "," << bounds[3] << endl;
        cout << bounds[4] << "," << bounds[5] << endl;
        
        int dims = opts.GetStringAsInt("-dims", 100);
        
        double maxbound = max(bounds[1]-bounds[0], max(bounds[3]-bounds[2], bounds[5]-bounds[4]));
        
        double gridSpacing = maxbound / dims;
        
        cout << "Grid Dimension: " << dims << "; Grid Spacing: " << gridSpacing << endl;
        
        
        size_t xdim = (bounds[1]-bounds[0])/gridSpacing;
        size_t ydim = (bounds[3]-bounds[2])/gridSpacing;
        size_t zdim = (bounds[5]-bounds[4])/gridSpacing;
        
        vtkStructuredGrid* grid = vtkStructuredGrid::New();
        grid->SetDimensions(xdim + 2, ydim + 2, zdim + 2);
        
        vtkPoints* gridPoints = vtkPoints::New();
        gridPoints->SetNumberOfPoints((xdim+2)*(ydim+2)*(zdim+2));
        //    gridPoints->SetNumberOfPoints(101*101*101);
        
        
        size_t u = 0;
        double x =bounds[0], y = bounds[2], z = bounds[4];
        for (int k = 0; k < zdim+2; k++) {
            for (int j = 0; j < ydim+2; j++) {
                for (int i = 0; i < xdim+2; i++) {
                    gridPoints->SetPoint(u, x, y, z);
                    x += gridSpacing;
                    u++;
                }
                y += gridSpacing;
                x = bounds[0];
            }
            z += gridSpacing;
            y = bounds[2];
        }
        
        grid->SetPoints(gridPoints);
        cout << "Grid construction done..." << endl;
        
        vtkSelectEnclosedPoints* encloserOut = vtkSelectEnclosedPoints::New();
        encloserOut->SetInput(grid);
        encloserOut->SetSurface(inputOut);
        encloserOut->CheckSurfaceOn();
        encloserOut->SetTolerance(0);
        cout << "Outside surface processing ..." << endl;
        encloserOut->Update();
        
        vtkDataArray* outLabel = encloserOut->GetOutput()->GetPointData()->GetArray("SelectedPoints");
        
        vtkSelectEnclosedPoints* encloserIn = vtkSelectEnclosedPoints::New();
        encloserIn->SetInput(grid);
        encloserIn->SetSurface(inputIn);
        encloserIn->CheckSurfaceOn();
        encloserIn->InsideOutOn();
        encloserIn->SetTolerance(0);
        cout << "Inside surface processing ..." << endl;
        encloserIn->Update();
        
        vtkDataArray* inLabel = encloserIn->GetOutput()->GetPointData()->GetArray("SelectedPoints");
        
        vtkIntArray* inOutLabel = vtkIntArray::New();
        inOutLabel->SetNumberOfComponents(1);
        inOutLabel->SetNumberOfValues(inLabel->GetNumberOfTuples());
        
        size_t insideCount = 0;
        cout << "Computing the intersection ..." << endl;
        for (size_t j = 0; j < inOutLabel->GetNumberOfTuples(); j++) {
            inOutLabel->SetValue(j, (outLabel->GetTuple1(j) == 1 && inLabel->GetTuple1(j) == 1) ? 1 : 0);
            insideCount++;
        }
        
        inOutLabel->SetName("InteriorPoints");
        grid->GetPointData()->SetScalars(inOutLabel);
        
        vio.writeFile(outputFile, grid);
        cout << "Inside Voxels: " << insideCount << endl;
    } else {
        string inputFile = args[0];
        string outputFile = args[1];
        
        vtkIO vio;
        vtkPolyData* input = vio.readFile(inputFile);
        int insideCount = 0;
        vtkDataSet* output = createGridForSphereLikeObject(input, insideCount, 100, true);
        vio.writeFile(outputFile, output);
        cout << "Inside Voxels: " << insideCount << endl;
    }
    
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
    
    vtkDataArray* interiorMarker = data->GetPointData()->GetArray(scalarName.c_str());
    if (interiorMarker == NULL) {
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
        int jInteriorMark = interiorMarker->GetTuple1(j);
        
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
        int surfaceInteriorExterior = jInteriorMark;
        unordered_set<vtkIdType>::const_iterator iter = nbrs.begin();
        for (; iter != nbrs.end(); iter++) {
            // if the neighbor is an exterior point
            vtkIdType nbrId = *iter;
            int nbrInteriorMark = interiorMarker->GetTuple1(nbrId);
            if (jInteriorMark == nbrInteriorMark) {
                continue;
            } else if (jInteriorMark == 0 && nbrInteriorMark == 1) {
                surfaceInteriorExterior = 10;
                break;
            } else if (jInteriorMark == 1 && nbrInteriorMark == 0) {
                // interior surface
                surfaceInteriorExterior = 11;
                break;
            } else {
                throw logic_error("invalid interior scalar makrer");
            }
        }
        
        boundaryMarker->SetTuple1(j, surfaceInteriorExterior);
    }
    
    boundaryMarker->SetName("BorderPoints");
    data->GetPointData()->AddArray(boundaryMarker);
    return boundaryMarker;
}

void runExtractStructuredGrid(Options& opts, StringVector& args) {
	string inputFile = args[0];
	string outputFile = args[1];
	
	vtkDataSet* inputPoly = vio.readDataFile(inputFile);
	double bounds[6];
	inputPoly->GetBounds(bounds);
	
	int extent[6];
	extent[0] = extent[2] = extent[4] = 0;
	extent[1] = (bounds[1]-bounds[0])/0.1;
	extent[3] = (bounds[3]-bounds[2])/0.1;
	extent[5] = (bounds[5]-bounds[4])/0.1;
	
	
	vtkNew<vtkExtractGrid> gridFilter;
	gridFilter->SetInput(inputPoly);
	gridFilter->IncludeBoundaryOn();
	gridFilter->SetVOI(extent);
	gridFilter->SetSampleRate(1, 1, 1);
	gridFilter->Update();
	
	vio.writeFile(outputFile, gridFilter->GetOutput());
	
}


// Compute Laplace PDE based on the adjacency list and border
void computeLaplacePDE(vtkDataSet* data, const double low, const double high, const int nIters, const double dt, vtkPolyData* surfaceData = NULL) {
	
	if (data == NULL) {
		cout << "Data input is NULL" << endl;
		return;
	}

	
	
	class LaplaceGrid {
	public:
		double low;
		double high;
		double dt;
		vtkDataSet* dataSet;
		vtkPolyData* samplePoints;
		
		vector<vtkIdType> solutionDomain;
		vtkIntArray* boundaryCond;
		vtkPolyData* boundarySurface;
		
		vtkDoubleArray* solution;
		vtkDoubleArray* tmpSolution;
		vtkDataArray* laplaceGradient;
		vtkDoubleArray* laplaceGradientNormals;

		
		Geometry geom;
		std::vector<Geometry::EdgeMap> edges;

		
        LaplaceGrid(double l, double h, double d, vtkDataSet* ds, vtkPolyData* pd= NULL): low(l), high(h), dt(d), dataSet(ds), boundarySurface(pd) {
			geom.extractEdges(ds, edges);
            
            // check boundary points
            boundaryCond = vtkIntArray::SafeDownCast(ds->GetPointData()->GetArray("SampledValue"));
            if (boundaryCond == NULL) {
                throw runtime_error("No scalar values for BoundaryPoints");
            }

            initializeSolution();
		}
        
		void initializeSolution() {
            // low-value 2
            // high-value 1
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
                int domain = boundaryCond->GetValue(j);
                double uValue = 0;
				if (domain == 700) {
					// high
                    uValue = high;
				} else if (domain == 300){
					// low
                    uValue = low;
				} else if (domain == 1) {
                    uValue = 0;
                    solutionDomain.push_back(j);
				}
                solution->SetValue(j, uValue);
				tmpSolution->SetValue(j, uValue);
			}
            cout << "# of solution grids: " << solutionDomain.size() << endl;
		}
		
		void computeStep() {
			const size_t nPts = solutionDomain.size();
			for (size_t j = 0; j < nPts; j++) {
				vtkIdType centerId = solutionDomain[j];
                Geometry::EdgeMap& edgeMap = edges[centerId];
                Geometry::EdgeMap::const_iterator iter = edgeMap.cbegin();
                
				double u = 0;
                size_t nNbrs = 0;
                for (; iter != edgeMap.cend(); iter++) {
                    if (iter->second.axisAligned == 2) {
                        const vtkIdType kId = iter->first;
                        const double du = solution->GetValue(kId);
                        u += du;
                        nNbrs ++;
                    }
//                    cout << iter->second.axisAligned << endl;
				}
				u = u/((double) nNbrs);
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
			laplaceGradient = gradFilter->GetOutput()->GetPointData()->GetArray("LaplacianGradient");
			
			laplaceGradientNormals = vtkDoubleArray::New();
			laplaceGradientNormals->SetName("LaplacianGradientNorm");
			laplaceGradientNormals->SetNumberOfComponents(3);
			laplaceGradientNormals->SetNumberOfTuples(laplaceGradient->GetNumberOfTuples());

			const size_t nPts = laplaceGradientNormals->GetNumberOfTuples();
			for (size_t j = 0; j < nPts; j++) {
				double* vec = laplaceGradient->GetTuple3(j);
				vtkMath::Normalize(vec);
				laplaceGradientNormals->SetTuple3(j, vec[0], vec[1], vec[2]);
			}
			
			data->GetPointData()->AddArray(laplaceGradient);
			data->GetPointData()->SetVectors(laplaceGradientNormals);
			laplaceGradientNormals->Delete();
		}
		
		void computeExteriorNormals(vtkPolyData* boundarySurface, const double radius = .1) {
			vtkNew<vtkPolyDataNormals> normalsFilter;
			normalsFilter->SetInput(boundarySurface);
			normalsFilter->ComputeCellNormalsOn();
			normalsFilter->ComputePointNormalsOn();
			normalsFilter->Update();
			vtkFloatArray* cellNormals = vtkFloatArray::SafeDownCast(normalsFilter->GetOutput()->GetCellData()->GetNormals());
			
			vtkNew<vtkCellLocator> cloc;
			cloc->SetDataSet(boundarySurface);
			cloc->AutomaticOn();
			cloc->BuildLocator();
			
			dataSet->GetPointData()->SetActiveScalars("SampledValue");
			
			vtkNew<vtkThresholdPoints> threshold;
			threshold->SetInput(dataSet);
			threshold->ThresholdByUpper(250);
			threshold->Update();
			vtkDataSet* inoutBoundary = threshold->GetOutput();
			vtkIntArray* inoutBoundaryCond = vtkIntArray::SafeDownCast(inoutBoundary->GetPointData()->GetArray("SampledValue"));

			
			vtkNew<vtkPointLocator> ploc;
			ploc->SetDataSet(inoutBoundary);
			ploc->AutomaticOn();
			ploc->BuildLocator();
			
			const size_t nPts = dataSet->GetNumberOfPoints();
			for (size_t j = 0; j < nPts; j++) {
				int domain = boundaryCond->GetValue(j);
				if (domain == 700 || domain == 300 || domain == 0) {
					double x[3] = { 0, }, closestPoint[3] = { 0, };
					vtkIdType cellId = -1;
					int subId = 0;
					double dist2 = -1;
					dataSet->GetPoint(j, x);
					vtkNew<vtkGenericCell> closestCell;
					cloc->FindClosestPointWithinRadius(x, radius, closestPoint, cellId, subId, dist2);
					
					float cellNormal[3];
					cellNormals->GetTupleValue(cellId, cellNormal);
					cellNormal[0] = 0;
					vtkMath::Normalize(cellNormal);

					if (domain == 0) {
						vtkIdType xId = ploc->FindClosestPoint(x);
						domain = inoutBoundaryCond->GetValue(xId);
						assert(domain == 300 || domain == 700);
					}
					

					if (domain == 300) {
						laplaceGradientNormals->SetTuple3(j, -cellNormal[0], -cellNormal[1], -cellNormal[2]);
					} else {
						laplaceGradientNormals->SetTuple3(j, cellNormal[0], cellNormal[1], cellNormal[2]);
					}
				}
			}
		}
	};
	
	
	LaplaceGrid grid(low, high, dt, data, surfaceData);
	
	// main iteration loop
	for (size_t i = 1; i <= nIters; i++) {
        if (i%500 == 0) {
            cout << "iteration: " << i << endl;
        }
		grid.computeStep();
	}
	
	
	// return the solution
	data->GetPointData()->AddArray(grid.solution);
	grid.computeNormals(data);
	grid.computeExteriorNormals(surfaceData);
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
	

	/// StreamTracer should have a point-wise gradient field
	/// - Set up tracer (Use RK45, both direction, initial step 0.05, maximum propagation 500
	StreamTracer* tracer = StreamTracer::New();
	tracer->SetInput(inputData);
	tracer->SetSource(inputSeedPoints);
	tracer->SetComputeVorticity(false);
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
	tracer->SetMaximumPropagation(500);
	tracer->SetInitialIntegrationStep(0.05);
	tracer->Update();
	
	
	vtkPolyData* streamLines = tracer->GetOutput();
//	streamLines->Print(cout);
	
	// remove useless pointdata information
	streamLines->GetPointData()->Reset();
	streamLines->BuildCells();
	streamLines->BuildLinks();
	
	
	// loop over the cell and compute the length
	int nCells = streamLines->GetNumberOfCells();
	
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
	
	cout << "Assigning length to each source vertex ..." << endl;
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







vtkDataSet* sampleSurfaceScalarsForGrid(vtkDataSet* ds, vtkDataSet* pd, string scalarName) {
    if (ds == NULL || pd == NULL) {
        throw runtime_error("Input is NULL");
    }
    vtkDataArray* scalars = NULL;
    if (scalarName == "") {
        scalars = pd->GetPointData()->GetScalars();
    } else {
        scalars = pd->GetPointData()->GetArray(scalarName.c_str());
    }
    if (scalars == NULL) {
        throw logic_error("No scalar available!");
    }
    
    
    vtkIntArray* sampledValue = vtkIntArray::New();
    sampledValue->SetName("SampledValue");
    sampledValue->SetNumberOfComponents(1);
    sampledValue->SetNumberOfTuples(ds->GetNumberOfPoints());
    
    vtkNew<vtkCellLocator> cloc;
    cloc->SetDataSet(pd);
    cloc->AutomaticOn();
    cloc->BuildLocator();
    
    
    
    vtkDataArray* insideOut = ds->GetPointData()->GetArray("BorderPoints");
	
	vtkNew<vtkIntArray> closestCell;
	closestCell->SetName("ClosestCell");
	closestCell->SetNumberOfComponents(1);
	closestCell->SetNumberOfValues(insideOut->GetNumberOfTuples());
	ds->GetPointData()->AddArray(closestCell.GetPointer());
	
    vtkNew<vtkGenericCell> genCell;
    for (size_t j = 0; j < insideOut->GetNumberOfTuples(); j++) {
        int jInsideOut = insideOut->GetTuple1(j);
        if (jInsideOut < 10) {
            closestCell->SetValue(j, -1);
            sampledValue->SetValue(j, jInsideOut);
            continue;
        } else if (jInsideOut == 11) {
            sampledValue->SetValue(j, 1);
            continue;
        }
        
        double jPt[3], x[3] = {0,};
        ds->GetPoint(j, jPt);
        
        vtkIdType cellId = -1;
        int subId = -1;
        double dist2 = -1;
        
        cloc->FindClosestPoint(jPt, x, cellId, subId, dist2);
		closestCell->SetValue(j, cellId);
        
        if (cellId == -1) {
            throw runtime_error("Can't find a closest cell");
        }
        vtkCell* cell = pd->GetCell(cellId);
        int values[3] = { 0, };
        for (size_t k = 0; k < cell->GetNumberOfPoints(); k++) {
            vtkIdType kId = cell->GetPointId(k);
            int scalarValue = scalars->GetTuple1(kId);
            values[scalarValue] ++;
        }
        if (j == 7970 || j == 8076 || j == 8182 || j == 8183 || j == 8289 || j == 8396) {
            cout << values[0] << "," << values[1] << "," << values[2] << endl;
            sampledValue->SetValue(j, 300);
        } else if (values[1] == 0 && values[2] > 0) {
            sampledValue->SetValue(j, 300);
        } else if (values[2] == 0 && values[1] > 0) {
            sampledValue->SetValue(j, 700);
        } else {
            sampledValue->SetValue(j, 300);
        }
//        cout << j << ": " <<  sampledValue->GetValue(j) << endl;
    }
    ds->GetPointData()->AddArray(sampledValue);
    
    Geometry geom;
    Geometry::EdgeList edges;
    geom.extractEdges(ds, edges);
	
    const size_t nPoints = edges.size();
    for (size_t j = 0; j < nPoints; j++) {
        Geometry::EdgeMap::const_iterator iter = edges[j].begin();
        for (; iter != edges[j].end(); iter++) {
            if (iter->second.u > iter->second.v) {
                continue;
            }
            vtkIdType uId = iter->second.u;
            vtkIdType vId = iter->second.v;
            int usv = insideOut->GetTuple1(uId);
            int vsv = insideOut->GetTuple1(vId);
            if ((usv == 0 && vsv == 1) || (usv == 1 && vsv == 0)) {
                throw logic_error("illegal boundary!");
            }
        }
    }
    return ds;
	

    
    vtkNew<vtkIdList> cells;
    
//    const size_t nPoints = edges.size();
    for (size_t j = 0; j < nPoints; j++) {
        Geometry::EdgeMap::const_iterator iter = edges[j].begin();
        for (; iter != edges[j].end(); iter++) {
            double uPt[3], vPt[3];
            if (iter->second.u > iter->second.v) {
                continue;
            }
            
            vtkIdType uId = iter->second.u;
            vtkIdType vId = iter->second.v;
            
            ds->GetPoint(uId, uPt);
            ds->GetPoint(vId, vPt);

            const int uIn = insideOut->GetTuple1(uId);
            const int vIn = insideOut->GetTuple1(vId);
            if (uIn == vIn) {
                continue;
            }
            
            vtkIdType insidePt = uIn > 0 ? iter->second.u : iter->second.v;
            vtkIdType outsidePt = uIn > 0 ? iter->second.v : iter->second.u;

            cells->Reset();
            cloc->FindCellsAlongLine(uPt, vPt, 0, cells.GetPointer());
            if (cells->GetNumberOfIds() > 1) {
                cout << uId << "," << vId << "," << cells->GetNumberOfIds() << " ";
				for (size_t k = 0; k < cells->GetNumberOfIds(); k++) {
					cout << cells->GetId(k) << " ";
				}
				cout << endl;
            }
            for (size_t k = 0; k < cells->GetNumberOfIds(); k++) {
                vtkIdType kId = cells->GetId(k);
                vtkCell* cell = pd->GetCell(kId);
                int values[3] = { 0, };
                for (size_t l = 0; l < cell->GetNumberOfPoints(); l++) {
                    int scalar = scalars->GetTuple1(cell->GetPointId(l));
                    values[scalar]++;
                }
                if (values[1] == 0 || values[2] == 0) {
                    int scalar = values[1] == 0 ? 2 : 1;
                    sampledValue->SetValue(outsidePt, scalar);
                } else {
                    sampledValue->SetValue(outsidePt, 3);
                }
            }
        }
    }
    return ds;
}


void runMeasureThickness(Options& opts, StringVector& args) {
    vtkIO vio;
    string inputFile = args[0];
    string outputStreamFile = args[1];
    
    vtkPolyData* input = vio.readFile(inputFile);
//    vtkPolyData* inputSeedPoints = vio.readFile(inputSeedFile);
    
//    vtkNew<vtkThresholdPoints> selector;
//    selector->SetInput(inputSeedPoints);
//    selector->ThresholdBetween(0.5, 1.5);
//    selector->SetInputArrayToProcess(0, 0, 0, vtkDataSet::FIELD_ASSOCIATION_POINTS, "meanLabels");
//    selector->Update();
//    vtkPolyData* selectedSeeds = selector->GetOutput();
    
//    cout << selectedSeeds->GetNumberOfPoints() << endl;
    
    int insideCount = 0;
    vtkDataSet* inOutGrid = createGridForSphereLikeObject(input, insideCount);
    cout << "Grid created..." << endl;
    selectBoundaryPoints(inOutGrid, "SelectedPoints");
    cout << "Boundary identified..." << endl;
    vtkDataSet* boundaryCondGrid = sampleSurfaceScalarsForGrid(inOutGrid, input, "meanLabels");
    cout << "Boundary condition assigned..." << endl;
    computeLaplacePDE(boundaryCondGrid, 0, 10000, 5000, 0.065, input);
    cout << "Laplace PDE computation done..." << endl;
    
    vtkPolyData* outputStream = performStreamTracer(opts, boundaryCondGrid, input);
    cout << "RK4 integration done..." << endl;
    vio.writeFile(outputStreamFile, outputStream);
    
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
	inputData->GetPointData()->SetActiveVectors("LaplacianGradientNorm");
	vtkPolyData* inputSeedPoints = vio.readFile(inputSeedPointsFile);
	
	vtkPolyData* outputStream = performStreamTracer(opts, inputData, inputSeedPoints, zRotate);
	
	vio.writeFile(outputPointFile, inputSeedPoints);
	if (outputStream) {
		vio.writeFile(outputStreamFile, outputStream);
	}
}


void processVolumeOptions(Options& opts) {
    opts.addOption("-markBorderCells", "Mark border cells of an input dataset. The border cells have 1 in BorderCells data", "-markBorderCells input-data output-data", SO_NONE);
    opts.addOption("-markBorderPoints", "Mark border points of an input dataset. The border points will be marked as 2 and its exterior neighbors will be marked as 3.", "-markBorderPoints input-data output-data", SO_NONE);
    
    opts.addOption("-extractBorderline", "Extract the borderlines between different labels", "-extractBorderline obj.vtp", SO_NONE);
	
    opts.addOption("-fillGrid", "Fill the inside of a polydata with a uniform grid (refer -twosided option)", "-fillGrid input.vtp output.vtp", SO_NONE);

	opts.addOption("-extractStructuredGrid", "Extract structured grid from a polydata ", "-extractStructuredGrid input.vtp output.vts", SO_NONE);
	
    opts.addOption("-dims", "x-y-z dimensions", "-dims 100", SO_REQ_SEP);
    opts.addOption("-twosided", "An option to generate the filled uniform grid", "-fillGrid CSF_GM_surface.vtk GM_WM_surface.vtk output.vts -twosided", SO_NONE);
	
    //
    opts.addOption("-sampleSurfaceScalarsForGrid", "Sample scalar values from a poly data at each grid point by finding a cell that intersects an edge of the grid", SO_NONE);

	// thickness measurement
	opts.addOption("-computeLaplacePDE", "Compute the Laplace PDE over the given domain", "-computeLaplacePDE input-data output-data ", SO_NONE);

	// RK45 stream tracer
	opts.addOption("-traceStream", "Trace a stream line from a given point set", "-traceStream input-vtu-field input-vtk output-lines output-points", SO_NONE);

    opts.addOption("-measureThickness", "Measure the thickness of the solution domain via RK45 integration", "-measureThickness input-polydata output-polydata", SO_NONE);

}

void processVolumeCommands(Options& opts, StringVector& args) {
    string input1File, outputFile;
    
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
    } else if (opts.GetBool("-extractBorderline")) {
        runExtractBorderline(opts, args);
    } else if (opts.GetBool("-fillGrid")) {
        runFillGrid(opts, args);
    } else if (opts.GetBool("-sampleSurfaceScalarsForGrid")) {
        input1File = args[0];
        string input2File = args[1];
        outputFile = args[2];
        
        string scalarName = opts.GetString("-scalarName", "");
        vtkDataSet* ds = vio.readDataFile(input1File);
        vtkPolyData* pd = vio.readFile(input2File);
        vtkDataSet* outDS = sampleSurfaceScalarsForGrid(ds, pd, scalarName);
        
        vio.writeFile(outputFile, outDS);
    } else if (opts.GetBool("-computeLaplacePDE")) {
		input1File = args[0];
		string input2File = args[1];
		outputFile = args[2];
		
		vtkDataSet* data = vio.readDataFile(input1File);
		vtkPolyData* surfaceData = vio.readFile(input2File);
		computeLaplacePDE(data, 0, 10000, 5000, 0.065, surfaceData);
		vio.writeFile(outputFile, data);
	} else if (opts.GetBool("-traceStream")) {
		// -traceStream 312.laplaceSol.vtp 312.sliceContour.vtk 312.thicknessSol.vtp 312.streams.vtp
		runStreamTracer(opts, args);
	} else if (opts.GetBool("-measureThickness")) {
		runMeasureThickness(opts, args);
	} else if (opts.GetBool("-extractStructuredGrid")) {
		runExtractStructuredGrid(opts, args);
	}
}
