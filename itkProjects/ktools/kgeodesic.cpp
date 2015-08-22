//
//  kgeodesic.cpp
//  ktools
//
//  Created by Joowhi Lee on 8/20/15.
//
//

#include "kgeodesic.h"
#include "vtkio.h"

#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyDataNormals.h>
#include <vtkMath.h>
#include <vtkTransform.h>

#include <vector>
#include <deque>

using namespace pi;
using namespace std;

static vtkIO vio;


struct ExpLogMap {
	vtkPoints* globalPoints;
	vtkPoints* centerPoints;
	vtkFloatArray* centerNormals;
	vtkIdTypeArray* ringId;
	
	vtkNew<vtkTransform> transform;
	
	ExpLogMap(vtkPoints* p, vtkPoints* cp, vtkFloatArray* cn, vtkIdTypeArray* ri): globalPoints(p), centerPoints(cp), centerNormals(cn), ringId(ri) {}
	
	double computeAngle(const double v1[3], const double v2[3]) {
		double cross[3] = { 0, };
		vtkMath::Cross(v1, v2, cross);
		double crossProdNorm = vtkMath::Norm(cross);
		double dotProd = vtkMath::Dot(v1, v2);
		return atan2(crossProdNorm, dotProd);
	}
	
	void transformExp(vtkPolygon* local, vtkIdType centerPtId, bool translate = false) {
		double centerPoint[3];
		double centerNormal[3];
		
		centerPoints->GetPoint(centerPtId, centerPoint);
		centerNormals->GetTuple(centerPtId, centerNormal);
		
		for (size_t j = 0; j < local->GetNumberOfPoints(); j++) {
			vtkIdType qId = local->GetPointId(j);
			vtkIdType gqId = ringId->GetValue(qId);
			
			// compute PQ vector
			double qPts[3], pqVec[3], Tpq[3];
			globalPoints->GetPoint(qId, qPts);
			vtkMath::Subtract(qPts, centerPoint, pqVec);
			
			double* qNormal = centerNormals->GetTuple3(gqId);
			double npDotNq = vtkMath::Dot(centerNormal, qNormal);
			
			// geodesic distance and angle
			double angle = computeAngle(centerNormal, pqVec) - M_PI_2;
//			angle = npDotNq < 0 ? M_PI - angle : angle;
			
			// map q to Tp
			double pqN[3] = { 0, };
			vtkMath::Cross(pqVec, centerNormal, pqN);
			vtkMath::Normalize(pqN);
			
			transform->Identity();
			transform->RotateWXYZ(vtkMath::DegreesFromRadians(angle), pqN);
			
			memcpy(Tpq, transform->TransformDoublePoint(pqVec), sizeof(Tpq));
			
			if (translate) {
				vtkMath::Add(Tpq, centerPoint, Tpq);
			}
			
			globalPoints->SetPoint(qId, Tpq);
			
		}
	}
};


void computeLocalTangentMap(vtkPolyData* g, vtkUnstructuredGrid* outputTangentMaps) {
	
//	vtkNew<vtkTransform> txfm;
//	
//	for (double angle = 90; angle <= 360; angle += 90) {
//		txfm->Identity();
//		txfm->RotateWXYZ(angle, 0, 0, 1);
//		
//		double x[3] = { 1, 0, 0 };
//		double *y = txfm->TransformDoublePoint(x);
//		
//		cout << angle << "; " << y[0] << "," << y[1] << "," << y[2] << endl;
//	}
//	
//	if (1) { return; }
	
	vtkNew<vtkPolyDataNormals> normalFilter;
	normalFilter->ComputePointNormalsOn();
	normalFilter->ConsistencyOn();
	normalFilter->SetInput(g);
	normalFilter->Update();
	g = normalFilter->GetOutput();
	
	const size_t nPts = g->GetNumberOfPoints();
	

    // the original points
	vtkPoints* globalPoints = vtkPoints::New();
	globalPoints->Allocate(nPts*5);
	
	vtkIdTypeArray* ringIds = vtkIdTypeArray::New();
    ringIds->SetName("RingID");
    ringIds->SetNumberOfComponents(1);
	ringIds->Allocate(nPts*5);
	
	
	// the cell array
    vtkCellArray* tangentMapArray = vtkCellArray::New();
    tangentMapArray->Allocate(nPts);
	
	vtkDoubleArray* centerPoints = vtkDoubleArray::New();
	centerPoints->SetName("CenterPoint");
	centerPoints->SetNumberOfComponents(3);
	centerPoints->SetNumberOfTuples(nPts);

	vtkFloatArray* centerNormals = vtkFloatArray::SafeDownCast(g->GetPointData()->GetNormals());
	centerNormals->SetName("CenterNormal");
	

	// Exponential Mapping
	ExpLogMap expMap(globalPoints, g->GetPoints(), centerNormals, ringIds);
	
    // temporary data types
    vtkNew<vtkIdList> cellIds, ptIds;
    vector<deque<vtkIdType> > neighbors;
    for (size_t p = 0; p < nPts; p++) {
		centerPoints->SetTupleValue(p, g->GetPoint(p));
		
        cellIds->Reset();
        ptIds->Reset();
		
        // inspect the neighbor cells of point p
        g->GetPointCells(p, cellIds.GetPointer());
		neighbors.clear();
        neighbors.resize(cellIds->GetNumberOfIds());
		
		if (p == 60) {
//			cout << cellIds->GetNumberOfIds() << endl;
			for (size_t j = 0; j < cellIds->GetNumberOfIds(); j++) {
				cout << cellIds->GetId(j) << endl;
			}
		}
        
        for (size_t j = 0; j < cellIds->GetNumberOfIds(); j++) {
            vtkIdType c = cellIds->GetId(j);
            g->GetCellPoints(c, ptIds.GetPointer());
            for (size_t k = 0; k < ptIds->GetNumberOfIds(); k++) {
                vtkIdType q = ptIds->GetId(k);
                neighbors[j].push_back(q);
            }
            while (neighbors[j].front() != p) {
                vtkIdType n = neighbors[j].front();
                neighbors[j].pop_front();
                neighbors[j].push_back(n);
            }
        }
        
        // construct a tangent plane with neighbor points
        vtkPolygon* tangentPlane = vtkPolygon::New();
		
        vtkIdType firstId = neighbors[0][1];
        vtkIdType lastId = -1;
        for (size_t j = 1; j < neighbors[0].size(); j++) {
            lastId = neighbors[0][j];
            ringIds->InsertNextValue(lastId);
            vtkIdType pId = globalPoints->InsertNextPoint(g->GetPoint(lastId));
            tangentPlane->GetPointIds()->InsertNextId(pId);
        }
        
        for (size_t j = 0; j < neighbors.size() && lastId != firstId; j++) {
            for (size_t k = 0; k < neighbors.size(); k++) {
                if (neighbors[k][1] == lastId) {
                    lastId = neighbors[k][2];
                    if (lastId == firstId) {
                        break;
                    }
                    ringIds->InsertNextValue(lastId);
                    vtkIdType pId = globalPoints->InsertNextPoint(g->GetPoint(lastId));
                    tangentPlane->GetPointIds()->InsertNextId(pId);
                }
            }
        }
		
		
		// compute the exponential mapping and add into the grid
		expMap.transformExp(tangentPlane, p, true);
		tangentMapArray->InsertNextCell(tangentPlane);

		
		if (p % 1000 == 0) {
			cout << "# of points processed: " << p << endl;
		}
    }
	
    outputTangentMaps->SetCells(VTK_POLYGON, tangentMapArray);
    outputTangentMaps->GetPointData()->SetScalars(ringIds);
    outputTangentMaps->SetPoints(globalPoints);
	outputTangentMaps->GetCellData()->SetVectors(centerPoints);
	outputTangentMaps->GetCellData()->AddArray(centerNormals);
}


// Compute geodesic distance to the set of destination points
void computeGeodesicDistance(vtkPolyData* data, vtkDataArray* dest) {
    vGraph vg(data);
    vg.buildAdjacencyList();
    
    vGraph::vtkIdVector source;
    for (size_t u = 0; u < vg.adjList.size(); u++) {
        int ul = dest->GetTuple1(u);
        for (size_t j = 0; j < vg.adjList[u].size(); j++) {
            size_t v = vg.adjList[u][j];
            int vl = dest->GetTuple1(v);
            if (vl != ul) {
                source.push_back(u);
                break;
            }
        }
    }
    
    vg.computeGeodesicDistance(source);
    
    vtkDoubleArray* dist = vtkDoubleArray::New();
    dist->SetName("distance");
    dist->SetNumberOfComponents(1);
    dist->SetArray(&vg.dist[0], vg.dist.size(), 1);
    
    data->GetPointData()->AddArray(dist);
}


void processGeodesicOptions(pi::Options& opts) {
    opts.addOption("-computeGeodesicDistance", "Compute geodesic distance to the set of points", "-computeGeodesicDistance input-data -scalarName destinationScalar", SO_NONE);
    
    opts.addOption("-computeLocalTangentMap", SO_NONE);
}

void processGeodesicCommands(pi::Options& opts, pi::StringVector& args) {
    if (opts.GetBool("-computeGeodesicDistance")) {
        string inputFile = args[0];
        string outputFile = args[1];
        
        string scalarName = opts.GetString("-scalarName", "BorderPoints");
        vtkPolyData* data = vio.readFile(inputFile);
        vtkDataArray* dest = data->GetPointData()->GetArray(scalarName.c_str());
        
        if (data == NULL || dest == NULL) {
            cout << "Input error!" << endl;
            return;
        }
        computeGeodesicDistance(data, dest);
        vio.writeFile(outputFile, data);
    } else if (opts.GetBool("-computeLocalTangentMap")) {
        string inputFile = args[0];
        string outputFile = args[1] == "" ? "localmap.vtu" : args[1];

        vtkPolyData* data = vio.readFile(inputFile.c_str());
        if (data == NULL) {
            cout << "Input error: " << inputFile << endl;
            return;
        }
        vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
        computeLocalTangentMap(data, ugrid);
        vio.writeFile(outputFile, ugrid);
    }
}