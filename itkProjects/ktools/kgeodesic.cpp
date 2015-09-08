//
//  kgeodesic.cpp
//  ktools
//
//  Created by Joowhi Lee on 8/20/15.
//
//

#include "kgeodesic.h"
#include "kgeometry.h"
#include "kdatastructure.h"

#include "vtkio.h"

#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyDataNormals.h>
#include <vtkMath.h>
#include <vtkTransform.h>
#include <vtkDijkstraGraphGeodesicPath.h>
#include <vtkLine.h>
#include <vtkPolyLine.h>
#include <vtkSphereSource.h>
#include <vtkGenericCell.h>
#include <vtkPointLocator.h>
#include <vtkSparseArray.h>

#include <vector>
#include <list>
#include <deque>
#include <queue>
#include <unordered_map>
#include <map>
#include <set>
#include <array>

using namespace pi;
using namespace std;

static vtkIO vio;
Geometry geom;

template <class T>
void print(T* p, size_t n = 3) {
	for (size_t j = 0; j < n; j++) {
		cout << p[j] << " ";
	}
}

class PolygonBuilder {
public:
	vtkNew<vtkUnstructuredGrid> grid;
	vtkNew<vtkPointLocator> pointLocator;
	
	PolygonBuilder() {
		grid->SetPoints(vtkPoints::New());
		grid->Allocate();
		pointLocator->SetDataSet(grid.GetPointer());
		//		pointLocator->BuildLocator();
	}
	
	vtkIdType addEdge(double* p1, double* p2) {
		vtkIdType ptIds[2];
		ptIds[0] = grid->GetPoints()->InsertNextPoint(p1);
		ptIds[1] = grid->GetPoints()->InsertNextPoint(p2);
		return grid->InsertNextCell(VTK_LINE, 2, ptIds);
	}
	
	vtkIdType addTriangle(double p1[3], double p2[3], double p3[3]) {
		vtkIdType ptIds[3];
		ptIds[0] = grid->GetPoints()->InsertNextPoint(p1);
		ptIds[1] = grid->GetPoints()->InsertNextPoint(p2);
		ptIds[2] = grid->GetPoints()->InsertNextPoint(p3);
		
		//        ptIds[0] = grid->FindPoint(p1);
		//        if (ptIds[0] < 0) {
		//            pointLocator->BuildLocator();
		//            pointLocator->InsertUniquePoint(p1, ptIds[0]);
		////            ptIds[0] = grid->GetPoints()->InsertNextPoint(p1);
		//        }
		//        ptIds[1] = grid->FindPoint(p2);
		//        if (ptIds[1] < 0) {
		//            pointLocator->InsertUniquePoint(p2, ptIds[1]);
		////            ptIds[1] = grid->GetPoints()->InsertNextPoint(p2);
		//        }
		//        ptIds[2] = grid->FindPoint(p2);
		//        if (ptIds[2] < 0) {
		//            pointLocator->InsertUniquePoint(p3, ptIds[2]);
		////            ptIds[2] = grid->GetPoints()->InsertNextPoint(p3);
		//        }
		return grid->InsertNextCell(VTK_TRIANGLE, 3, ptIds);
	}
	
	vtkIdType addCell(vtkCell* cell, vtkPoints* pts, vtkTransform* txfm) {
		vector<vtkIdType> newIds(cell->GetNumberOfPoints());
		for (size_t j = 0; j < newIds.size(); j++) {
			if (txfm == NULL) {
				newIds[j] = grid->GetPoints()->InsertNextPoint(pts->GetPoint(cell->GetPointId(j)));
			} else {
				newIds[j] = grid->GetPoints()->InsertNextPoint(txfm->TransformPoint(pts->GetPoint(cell->GetPointId(j))));
			}
		}
		return grid->InsertNextCell(cell->GetCellType(), newIds.size(), &newIds[0]);
	}
	
	void addFan(double *ct, vtkCell* neighbors) {
		for (size_t j = 0; j < neighbors->GetNumberOfPoints(); j++) {
			
		}
	}
	
	void writeToFile(string filename) {
		vtkIO vio;
		vio.writeFile(filename, grid.GetPointer());
	}
	
	vtkUnstructuredGrid* GetOutput() {
		grid->SetPoints(pointLocator->GetPoints());
		return grid.GetPointer();
	}
};

class LocalTangentPlane {
public:
	struct PointTag {
		vtkIdType uId;
		vtkIdType vId;
		vtkIdType uGlobalId;
		vtkIdType vGlobalId;
		double uvAngle;
		double uNormal[3];
		bool tag;
		inline bool isCenter() { return uId == vId; }
	};
	
	vtkPolyData* dataSet;
	const size_t nPoints;
	unordered_map<vtkIdType, LocalSurfaceTopology> topologyMap;
	vector<vector<vtkIdType> > neighborIds;
	vtkPolyData* geometry, *expGeom;
	vtkFloatArray* normals;
    
    typedef unordered_map<vtkIdType, PointTag> NeighborInfoType;
	vector<NeighborInfoType > globalIds;
	
	LocalTangentPlane(vtkPolyData* ds): dataSet(ds), nPoints(dataSet->GetNumberOfPoints()) {
        neighborIds.resize(nPoints);
        globalIds.resize(nPoints);
        
		buildTopology();
		computeNormals(dataSet);
		computeGeometry();
		normalizelGeometry();
		computeRotation();
	}
	
	vtkFloatArray* computeNormals(vtkPolyData* ds) {
		vtkNew<vtkPolyDataNormals> normalFilter;
		normalFilter->ComputePointNormalsOn();
		normalFilter->ConsistencyOn();
		normalFilter->SplittingOff();
		normalFilter->SetInput(dataSet);
		normalFilter->Update();
		vtkPolyData* g = normalFilter->GetOutput();
		normals = vtkFloatArray::New();
		normals->DeepCopy(g->GetPointData()->GetNormals());
		return normals;
	}
	
	vtkPolyData* getGeometry() {
		return geometry;
	}
	
	vtkPolyData* getLocalGeometry() {
		return expGeom;
	}
	
	void addEdge(vtkIdType j, vtkIdType u, vtkIdType v) {
		topologyMap[j].setCenterId(j);
		topologyMap[j].addEdge(u, v);
	}
	
	void buildTopology() {
		vtkNew<vtkIdList> cellIds;
		for (size_t j = 0; j < nPoints; j++) {
			dataSet->GetPointCells(j, cellIds.GetPointer());
			for (size_t k = 0; k < cellIds->GetNumberOfIds(); k++) {
				vtkCell* cell = dataSet->GetCell(cellIds->GetId(k));
				for (size_t h = 0; h < cell->GetNumberOfEdges(); h++) {
					vtkCell* edge = cell->GetEdge(h);
					vtkIdType u = edge->GetPointId(0);
					vtkIdType v = edge->GetPointId(1);
					if (u != j && v != j) {
						addEdge(j, u, v);
					}
				}
			}
		}
		
		for (size_t j = 0; j < nPoints; j++) {
			topologyMap[j].getOrderedNeighbors(neighborIds[j]);
		}
	}
	
	void computeGeometry() {
		vtkPoints* outputPoints = vtkPoints::New();
		vtkPolyData* output = vtkPolyData::New();
		output->Allocate();
		output->SetPoints(outputPoints);
		
		// source Cell Ids
		vtkIntArray* sourceIds = vtkIntArray::New();
		sourceIds->Allocate(nPoints*5);
		sourceIds->SetNumberOfComponents(1);
		sourceIds->SetName("SourceID");
		output->GetCellData()->SetScalars(sourceIds);
		
		
		vtkIntArray* cellIndex = vtkIntArray::New();
		cellIndex->SetName("CellIndex");
		cellIndex->SetNumberOfComponents(2);
		cellIndex->SetNumberOfTuples(nPoints);
		output->GetFieldData()->AddArray(cellIndex);
		
		// source Point Ids
		vtkIntArray* ringIds = vtkIntArray::New();
		ringIds->Allocate(dataSet->GetNumberOfPoints()*5);
		ringIds->SetNumberOfComponents(3);
		ringIds->SetName("RingID");
		output->GetPointData()->AddArray(ringIds);
		
		
		// exp transform
		vtkNew<vtkTransform> expTxfm;
		expTxfm->PostMultiply();
		
		// assume there is no boundary
		unordered_map<vtkIdType, vtkIdType> newIds;
		for (size_t j = 0; j < nPoints; j++) {
			assert(!topologyMap[j].boundary);
			newIds.clear();
			vtkIdType centerGlobalId = ringIds->GetMaxId() + 1;
			
			double centerPt[3];
			dataSet->GetPoint(j, centerPt);
			
			double* nu = normals->GetTuple3(j);
			
			newIds[j] = outputPoints->InsertNextPoint(dataSet->GetPoint(j));
			ringIds->InsertNextTuple3(j, j, centerGlobalId);
			
			globalIds[j][j].uId = globalIds[j][j].vId = j;
			globalIds[j][j].uGlobalId = globalIds[j][j].vGlobalId = newIds[j];
			
			for (size_t k = 0; k < neighborIds[j].size() - 1; k++) {
				vtkIdType nbrId = neighborIds[j][k];
				ringIds->InsertNextTuple3(nbrId, j, centerGlobalId);
				
				double uv[3] = { 0, }, expuv[3] = { 0, };
				dataSet->GetPoint(nbrId, uv);
				
				expTxfm->Identity();
				geom.tangentVector(centerPt, uv, nu, expuv, expTxfm.GetPointer());
				
				newIds[nbrId] = outputPoints->InsertNextPoint(expuv);
				
				globalIds[j][nbrId].uId = j;
				globalIds[j][nbrId].vId = nbrId;
				globalIds[j][nbrId].uGlobalId = newIds[j];
				globalIds[j][nbrId].vGlobalId = newIds[nbrId];
			}
			for (size_t k = 0; k < neighborIds[j].size() - 1; k++) {
				vtkIdType tri[3];
				tri[0] = newIds[j];
				tri[1] = newIds[neighborIds[j][k]];
				tri[2] = newIds[neighborIds[j][k+1]];
				vtkIdType newCellId = output->InsertNextCell(VTK_TRIANGLE, 3, tri);
				sourceIds->InsertNextValue(j);
			}
			cellIndex->SetTuple2(j, newIds[j], neighborIds[j].size());
		}
		
		geometry = output;
	}
	
	void normalizelGeometry() {
		expGeom = vtkPolyData::New();
		expGeom->DeepCopy(geometry);
		
		vtkIntArray* cellIndex = vtkIntArray::SafeDownCast(expGeom->GetFieldData()->GetArray("CellIndex"));
		
		vtkNew<vtkTransform> txfm;
		for (size_t j = 0; j < nPoints; j++) {
			int cellIdx[2];
			cellIndex->GetTupleValue(j, cellIdx);
			
			double u[3], n[3], cross[3] = { 0, };
			expGeom->GetPoint(cellIdx[0], u);
			normals->GetTuple(j, n);
			
			txfm->Identity();
			geom.normalizeToNorthPole(u, n, cross, txfm.GetPointer());
			expGeom->GetPoints()->SetPoint(cellIdx[0], txfm->TransformPoint(u));
			
			for (vtkIdType k = cellIdx[0] + 1; k < cellIdx[0] + cellIdx[1]; k++) {
				double kPt[3], tkPt[3];
				expGeom->GetPoint(k, kPt);
				txfm->TransformPoint(kPt, tkPt);
				expGeom->GetPoints()->SetPoint(k, tkPt);
			}
		}
	}
	
	void denormalizeGeometry() {
		vtkIntArray* cellIndex = vtkIntArray::SafeDownCast(expGeom->GetFieldData()->GetArray("CellIndex"));
		
		vtkNew<vtkTransform> txfm;
		for (size_t j = 0; j < nPoints; j++) {
			int cellIdx[2];
			cellIndex->GetTupleValue(j, cellIdx);
			
			double u[3], n[3], cross[3] = { 0, };
			dataSet->GetPoint(j, u);
			normals->GetTuple(j, n);
			
			txfm->Identity();
			geom.denormalizeFromNorthPole(u, n, cross, txfm.GetPointer());
			
			expGeom->GetPoint(cellIdx[0], u);
			expGeom->GetPoints()->SetPoint(cellIdx[0], txfm->TransformPoint(u));
			
			for (vtkIdType k = cellIdx[0] + 1; k < cellIdx[0] + cellIdx[1]; k++) {
				double kPt[3], tkPt[3];
				expGeom->GetPoint(k, kPt);
				txfm->TransformPoint(kPt, tkPt);
				expGeom->GetPoints()->SetPoint(k, tkPt);
			}
		}
	}
	
	void computeRotation() {
		vtkNew<vtkTransform> txf;
		for (size_t u = 0; u < nPoints; u++) {
			for (size_t v = 0; v < neighborIds[u].size() - 1; v++) {
				vtkIdType vId = neighborIds[u][v];
				PointTag& utag = globalIds[u][vId];
				PointTag& vtag = globalIds[vId][u];
				
				// compute the rotation to match along v-u
				double u1Pt[3], v1Pt[3], u2Pt[3], v2Pt[3];
				expGeom->GetPoint(utag.uGlobalId, u1Pt);
				expGeom->GetPoint(utag.vGlobalId, v1Pt);
				expGeom->GetPoint(vtag.uGlobalId, u2Pt);
				expGeom->GetPoint(vtag.vGlobalId, v2Pt);
				
				vtkMath::MultiplyScalar(v2Pt, -1);
				double cross[3];
				double angDeg = geom.rotateVector(v1Pt, v2Pt, txf.GetPointer(), cross);
				txf->TransformPoint(v1Pt, v1Pt);
				
				// the rotation angle must be counter clock-wise
				if (cross[2] < 0) {
					angDeg *= -1;
				}
				
				utag.uvAngle = angDeg;
				normals->GetTuple(u, utag.uNormal);
			}
		}
	}
	
	double findShortestPath(vtkIdType s, vtkIdType e, vtkPolyData* visLocalMap = NULL) {
		bool denormalizedVis = false;
		
		vtkNew<vtkDijkstraGraphGeodesicPath> dijkstraFilter;
		dijkstraFilter->SetInput(dataSet);
		dijkstraFilter->SetStartVertex(s);
		dijkstraFilter->SetEndVertex(e);
		dijkstraFilter->Update();
		vtkPolyData* shortestPath = dijkstraFilter->GetOutput();
		//        vio.writeFile("shortest_path.vtp", shortestPath);
		
		double dijkstraLength = 0;
		vtkIdList* pathElems = dijkstraFilter->GetIdList();
		
		
		// create a poly data with points from the expgeom
		
		vtkPoints* visPoints = vtkPoints::New();
		visPoints->DeepCopy(expGeom->GetPoints());
		
		double geodesicVector[3] = { 0, };
		vtkNew<vtkPoints> deltas;
		vtkNew<vtkIdList> localCells;
		vtkIntArray* sourceIds = vtkIntArray::SafeDownCast(expGeom->GetCellData()->GetArray("SourceID"));
		
		if (visLocalMap) {
			visLocalMap->SetPoints(visPoints);
			visLocalMap->Allocate();

			vtkIntArray* ringIds = vtkIntArray::SafeDownCast(expGeom->GetPointData()->GetArray("RingID"));
			
			vtkIntArray* visRingIds = vtkIntArray::New();
			visRingIds->DeepCopy(ringIds);
			visLocalMap->GetPointData()->AddArray(ringIds);
		}
		

		
		vtkNew<vtkTransform> txfm;
		vtkNew<vtkTransform> visTxfm;
		txfm->PostMultiply();
		visTxfm->PostMultiply();
		
		vtkNew<vtkIntArray> cellGroup;
		cellGroup->SetName("CellGroup");
		
		vtkNew<vtkIdList> cellList;
		
		double totAngle = 0;
		
		vtkIdType prevU = -1;
		for (size_t j = 0; j < pathElems->GetNumberOfIds() - 1; j++) {
			vtkIdType u = pathElems->GetId(j);
			vtkIdType v = pathElems->GetId(j+1);
			
			double uSurfPt[3], vSurfPt[3];
			dataSet->GetPoint(u, uSurfPt);
			dataSet->GetPoint(v, vSurfPt);
			double edgeDist = sqrt(vtkMath::Distance2BetweenPoints(uSurfPt, vSurfPt));
			dijkstraLength += edgeDist;
			
			
			double uvVec[3] = { 0, };
			expGeom->GetPoint(globalIds[u][v].vGlobalId, uvVec);
			txfm->TransformPoint(uvVec, uvVec);
			
			deltas->InsertNextPoint(uvVec);
			vtkMath::Add(geodesicVector, uvVec, geodesicVector);
			
			
			cellList->Reset();
			if (visLocalMap) {
				sourceIds->LookupValue(u, cellList.GetPointer());
				for (size_t k = 0; k < cellList->GetNumberOfIds(); k++) {
					vtkCell* cell = expGeom->GetCell(cellList->GetId(k));
					vtkIdType newId = visLocalMap->InsertNextCell(cell->GetCellType(), cell->GetPointIds());
					for (size_t l = 0; l < cell->GetNumberOfPoints(); l++) {
						vtkIdType cellPointId = cell->GetPointId(l);
						double cellPt[3];
						expGeom->GetPoint(cellPointId, cellPt);
						visTxfm->TransformPoint(cellPt, cellPt);
						visPoints->SetPoint(cellPointId, cellPt);
						
					}
					cellGroup->InsertValue(newId, j);
				}
			}
			
			
			double ang = globalIds[v][u].uvAngle;
			totAngle += ang;
			//            cout << "u: " << u << " v: " << v << "; " << ang << "; " << globalIds[u][v].uvAngle << endl;
			txfm->RotateZ(ang);
			
			
			
			if (visLocalMap) {
				visTxfm->Identity();
				visTxfm->RotateZ(totAngle);
				
				if (denormalizedVis) {
					geom.rotateVector(_northPole, normals->GetTuple3(v), visTxfm.GetPointer());
				}
				visTxfm->Translate(geodesicVector);
				
				if (denormalizedVis) {
					vtkMath::Subtract(vSurfPt, geodesicVector, vSurfPt);
					visTxfm->Translate(vSurfPt);
				}
			}
		}
		
		if (visLocalMap) {
			if (pathElems->GetNumberOfIds() > 1) {
				size_t j = pathElems->GetNumberOfIds() - 1;
				vtkIdType u = pathElems->GetId(pathElems->GetNumberOfIds() - 1);
				
				
				cellList->Reset();
				sourceIds->LookupValue(u, cellList.GetPointer());
				for (size_t k = 0; k < cellList->GetNumberOfIds(); k++) {
					vtkCell* cell = expGeom->GetCell(cellList->GetId(k));
					vtkIdType newId = visLocalMap->InsertNextCell(cell->GetCellType(), cell->GetPointIds());
					for (size_t l = 0; l < cell->GetNumberOfPoints(); l++) {
						vtkIdType cellPointId = cell->GetPointId(l);
						double cellPt[3];
						expGeom->GetPoint(cellPointId, cellPt);
						visTxfm->TransformPoint(cellPt, cellPt);
						visPoints->SetPoint(cellPointId, cellPt);
					}
					cellGroup->InsertValue(newId, j);
				}
			}
			
			if (visLocalMap) {
				double p[3], n[3], cross[3];
				dataSet->GetPoint(pathElems->GetId(0), p);
				normals->GetTuple(pathElems->GetId(0), n);
				
				visTxfm->Identity();
				visTxfm->PostMultiply();
				
				if (!denormalizedVis) {
					geom.denormalizeFromNorthPole(p, n, cross, visTxfm.GetPointer());
				}
				
				for (size_t j = 0; j < visPoints->GetNumberOfPoints(); j++) {
					visPoints->SetPoint(j, visTxfm->TransformPoint(visPoints->GetPoint(j)));
				}
			}
			//			cout << "Geodesic Distance: " << vtkMath::Norm(geodesicVector) << endl;
			//			cout << "Dijkstra Distance: " << dijkstraLength << endl;			return vtkMath::Norm(geodesicVector);
		}
		
		visPoints->Delete();
		return vtkMath::Norm(geodesicVector);
	}
	
	
	void getTangentVector(vtkIdType u, vtkIdType v, double uv[3]) {
		vtkIdType uvPtId = globalIds[u][v].vGlobalId;
		expGeom->GetPoint(uvPtId, uv);
	}
	
	PointTag& getPointTag(vtkIdType u, vtkIdType v) {
		return globalIds[u][v];
	}
	
	void print(ostream& os) {
		vector<vtkIdType> boundaries;
		for (size_t j = 0; j < nPoints; j++) {
			topologyMap[j].getOrderedNeighbors(boundaries);
			cout << j << ": ";
			for (size_t k = 0; k < boundaries.size(); k++) {
				cout << boundaries[k] << " ";
			}
			cout << endl;
		}
	}
};


void computeGeodesicDistance(vtkDataSet* dataSet, vtkDataArray* dest) {
	vtkNew<vtkIdList> ids;
	dest->LookupValue(1, ids.GetPointer());
	
	for (size_t j = 0; j < ids->GetNumberOfIds(); j++) {
		cout << ids->GetId(j) << endl;
	}
}

struct SurfaceNode {
	vtkIdType id;
	double dist;
	
	SurfaceNode(): id(0), dist(DBL_MAX) {}
	SurfaceNode(vtkIdType d): id(d), dist(DBL_MAX) {}
	SurfaceNode(vtkIdType d, double dst): id(d), dist(dst) {}
	bool operator<(const SurfaceNode& o) const {
		return this->dist > o.dist;
	}
	bool operator==(const SurfaceNode& o) const {
		return this->id == o.id;
	}
};


vtkPolyData* computeSingleSourceGeodesicDistance0(vtkPolyData* data, vtkIdType sourceId) {
	LocalTangentPlane tangentPlane(data);
	
	vtkDoubleArray* dist = vtkDoubleArray::New();
	dist->SetName("GeodesicDistance");
	dist->SetNumberOfComponents(1);
	dist->SetNumberOfValues(tangentPlane.nPoints);

	for (size_t k = 0; k < tangentPlane.nPoints; k++) {
		if ((k)%100 == 0) cout << "Processing " << (k+1) << " points ..." << endl;
		if (k == sourceId) {
			dist->SetValue(k, 0);
			continue;
		}
		double distSP = tangentPlane.findShortestPath(sourceId, k);
		dist->SetValue(k, distSP);
	}
	data->GetPointData()->SetScalars(dist);
	
	return data;
}


// perform BFS and compute geodesic distance to every other point
vtkPolyData* computeSingleSourceGeodesicDistance1(vtkPolyData* data, vtkIdType sourceId) {
	const int NEWNODE = 0;
	const int VISITED = 1;
	const int QUEUED = 2;
	
	LocalTangentPlane tangentPlane(data);
	
	unordered_set<vtkIdType> sourceIds;
	sourceIds.insert(sourceId);
	
	vector<vtkIdType> visited(tangentPlane.nPoints);
	fill(visited.begin(), visited.end(), 0);
	
	vtkDoubleArray* geodesicVectors = vtkDoubleArray::New();
	geodesicVectors->SetName("GeodesicVectors");
	geodesicVectors->SetNumberOfComponents(3);
	geodesicVectors->SetNumberOfTuples(tangentPlane.nPoints);
    for (size_t j = 0; j < geodesicVectors->GetNumberOfComponents(); j++) {
        geodesicVectors->FillComponent(j, 0);
    }
	
	
    vtkNew<vtkDoubleArray> geodesicDistance;
	geodesicDistance->SetName("GeodesicDistance");
	geodesicDistance->SetNumberOfComponents(1);
	geodesicDistance->SetNumberOfValues(tangentPlane.nPoints);
	geodesicDistance->FillComponent(0, DBL_MAX);
	data->GetPointData()->SetScalars(geodesicDistance.GetPointer());
	
	vector<double> geodesicAngle(tangentPlane.nPoints);
    
    
    struct DijkstraNodeQueue {
        typedef FibQueue<double> FibQueueD;
        typedef FibQueueD::FibNode FibNodeD;
        typedef pair<vtkIdType, FibNodeD*> NodeType;
        vector<NodeType> nodeMap;
        
        FibQueueD fq;
        
        DijkstraNodeQueue(size_t nNodes) {
            nodeMap.resize(nNodes);
            for (size_t j = 0; j < nodeMap.size(); j++) {
                nodeMap[j].first = -1;
                nodeMap[j].second = NULL;
            }
        }
        
        void push(vtkIdType n, double d) {
            if (hasKey(n)) {
                decreaseKey(n, d);
            } else {
                nodeMap[n] = NodeType(n, fq.push(d));
                nodeMap[n].second->payload = &nodeMap[n];
            }
        }
        
        vtkIdType pop() {
            FibNodeD* fn = fq.topNode();
            NodeType* tn = (NodeType*) fn->payload;
            vtkIdType retId = tn->first;
            nodeMap[retId].first = -1;
            nodeMap[retId].second = NULL;
            fq.pop();
            return retId;
        }
        
        bool hasKey(vtkIdType n) {
            return nodeMap[n].second != NULL;
        }
        
        void decreaseKey(vtkIdType k, double d) {
            FibNodeD* fn = nodeMap[k].second;
            fq.decrease_key(fn, d);
            return;
        }
        
        bool empty() {
            return fq.empty();
        }
    };
	
    
    
    DijkstraNodeQueue que(tangentPlane.nPoints);
    que.push(sourceId, 0);
    
	geodesicDistance->SetValue(sourceId, 0);

	vtkNew<vtkTransform> txfm;
	txfm->PostMultiply();
	txfm->Identity();
    
    
    // record for the previous path
    vector<vtkIdType> prevPath(tangentPlane.nPoints);
    fill(prevPath.begin(), prevPath.end(), 0);
    prevPath[sourceId] = sourceId;
	
	while (!que.empty()) {
        vtkIdType fId = que.pop();
		
		if (visited[fId] == VISITED) {
			continue;
		}
		visited[fId] = VISITED;
		
        double fVec[3];
        geodesicVectors->GetTupleValue(fId, fVec);
        
		double fPts[3];
		data->GetPoint(fId, fPts);
		txfm->RotateZ(geodesicAngle[fId]);
		
		vector<vtkIdType>& fNbrs = tangentPlane.neighborIds[fId];
		for (size_t j = 0; j < fNbrs.size() - 1; j++) {
			vtkIdType fnId = fNbrs[j];
			if (sourceIds.find(fNbrs[j]) != sourceIds.end()) {
				continue;
			}
			
			double fnPts[3];
			data->GetPoint(fnId, fnPts);
			
			// relax
			double delta[3];
			txfm->Push();
			tangentPlane.getTangentVector(fId, fnId, delta);
			
			double localAngle = tangentPlane.globalIds[fnId][fId].uvAngle;
			txfm->RotateZ(localAngle);
			txfm->TransformPoint(delta, delta);
			txfm->Pop();
            
            delta[2] = 0;
			double newVector[3];
			vtkMath::Add(fVec, delta, newVector);
			
			
			double distViaF = vtkMath::Norm(newVector);
			if (distViaF < geodesicDistance->GetValue(fnId)) {
				geodesicDistance->SetValue(fnId, distViaF);
				geodesicAngle[fnId] = geodesicAngle[fId] + localAngle;
				geodesicVectors->SetTupleValue(fnId, newVector);
                prevPath[fnId] = fId;
			}
			
			if (visited[fnId] != VISITED) {
                que.push(fnId, distViaF);
				visited[fnId] = QUEUED;
			}
            if (fnId == 39) {
                double fDist = geodesicDistance->GetValue(fId);
                double gF[3], gnF[3];
                geodesicVectors->GetTupleValue(fId, gF);
                geodesicVectors->GetTupleValue(fnId, gnF);
                cout << fId << " > " << fnId << "; " << fDist << "; " << distViaF << "; " << gF[0] << "," << gF[1] << "," << gF[2] << "; " << delta[0] << "," << delta[1] << "," << delta[2] << "; " << gnF[0] << "," << gnF[1] << "," << gnF[2] << endl;
            }
		}
//        break;
        cout << "dist: " << geodesicDistance->GetValue(39) << endl;
	}
	
    vtkIdType pp = 39;
    while (prevPath[pp] != pp) {
        cout << pp << " ";
        pp = prevPath[pp];
    }
    cout << endl;
	return data;
}

vtkPolyData* computeSingleSourceGeodesicDistance(vtkPolyData* data, vtkIdType sourceId) {
    LocalTangentPlane tangentPlane(data);
    
    vtkNew<vtkDoubleArray> dist;
    dist->SetName("Distance");
    dist->SetNumberOfValues(tangentPlane.nPoints);
    data->GetPointData()->SetScalars(dist.GetPointer());
    
    vtkNew<vtkDoubleArray> angle;
    angle->SetName("Angle");
    angle->SetNumberOfValues(tangentPlane.nPoints);
    data->GetPointData()->AddArray(angle.GetPointer());
    
    vtkNew<vtkDoubleArray> vect;
    vect->SetName("Vector");
    vect->SetNumberOfComponents(3);
    vect->SetNumberOfTuples(tangentPlane.nPoints);
    data->GetPointData()->SetVectors(vect.GetPointer());
    
    for (size_t j = 0; j < tangentPlane.nPoints; j++) {
        dist->SetValue(j, DBL_MAX);
        vect->SetTuple3(j, 0, 0, 0);
        angle->SetValue(j, 0);
    }
    
    dist->SetValue(sourceId, 0);
    angle->SetValue(sourceId, 0);
    
    
    vtkNew<vtkTransform> txfm;
    txfm->PostMultiply();
    
    for (size_t j = 0; j < tangentPlane.nPoints; j++) {
        if (j % 100 == 0) {
            cout << "Processing " << (j/tangentPlane.nPoints)*100 << "% ..." << endl;
        }
        for (size_t k = 0; k < tangentPlane.nPoints; k++) {
            vtkIdType e1 = k;
            
            // the accumulated angle up to e1
            float e1ang = angle->GetValue(e1);
            txfm->Identity();
            txfm->RotateZ(e1ang);

            // the geodesic vector up to e1
            double e1Vec[3];
            vect->GetTuple(e1, e1Vec);
            
            double e1dist = dist->GetValue(e1);
            if (e1dist == DBL_MAX) {
                continue;
            }
            
            vector<vtkIdType>& nbrIds = tangentPlane.neighborIds[k];
            const size_t nNbrIds = nbrIds.size();
            LocalTangentPlane::NeighborInfoType& nbrs = tangentPlane.globalIds[e1];
            
            for (size_t l = 0; l < nNbrIds - 1; l++) {
                vtkIdType e2 = nbrIds[l];
                
                double tv[3];
                tangentPlane.getTangentVector(e1, e2, tv);
                
                double e2Vec[3];
                vtkMath::Add(e1Vec, txfm->TransformPoint(tv), e2Vec);

                double e2dist = dist->GetValue(e2);
                double e2NewDist = vtkMath::Norm(e2Vec);
                
                if (e2NewDist < e2dist) {
                    dist->SetValue(e2, e2NewDist);
                    vect->SetTuple(e2, e2Vec);
                    angle->SetValue(e2, e1ang + nbrs[e1].uvAngle);
                }
            }
        }
    }
    return data;
}

void testGeodesic() {
	if (1) {
		vtkPolyData* plane = vio.readFile("cylinder_sample_clean.vtk");
		computeSingleSourceGeodesicDistance(plane, 0);
		vio.writeFile("cylinder_sample_dist.vtp", plane);
		return;
	}
	
	if (0) {
		LocalSurfaceTopology topo(0);
		topo.addEdge(1, 2);
		topo.addEdge(2, 3);
		topo.addEdge(4, 5);
		topo.addEdge(3, 4);
		topo.addEdge(5, 1);
		vector<vtkIdType> ids;
		topo.getOrderedNeighbors(ids);
		for (size_t j = 0; j < ids.size(); j++) {
			cout << ids[j] << " ";
		}
		cout << endl;
	}
	if (0) {
		vtkNew<vtkTransform> s;
		s->RotateZ(10);
		s->Push();
		s->RotateZ(80);
		double p1[3] = { 1, 0, 0 };
		double p2[3] = { 0, };
		s->TransformPoint(p1, p2);
		print(p1); print(p2);
		s->Pop();
		s->TransformPoint(p1, p2);
		print(p1); print(p2);
	}
	if (0) {
		vtkIO vio;
		vtkPolyData* sphere = vio.readFile("sphere_sample_8x8.vtk");
		
		LocalTangentPlane localMap(sphere);
		vtkPolyData* visLocalMap = vtkPolyData::New();
		cout << "distance: " << localMap.findShortestPath(0, 4, visLocalMap) << endl;
		vio.writeFile("sphere_sample_8x8_dist.vtp", visLocalMap);
		
		return;
		
		for (size_t j = 0; j < 10; j++) {
			for (size_t k = 0; k < 10; k++) {
				if (j == k) {
					continue;
				}
				char filename[256];
				sprintf(filename, "shortest_path_%03lu_to_%03lu.vtp", k, j);
				vtkPolyData* visLocalMap = vtkPolyData::New();
				localMap.findShortestPath(j, k, visLocalMap);
				vio.writeFile(string(filename), visLocalMap);
			}
		}
	}
	
	if (0) {
		vtkNew<vtkSparseArray<double> > sparseArr;
		sparseArr->Resize(10, 10);
		sparseArr->SetValue(2,3,9);
		sparseArr->SetValue(4,7,13);
		sparseArr->Print(cout);
		cout << sparseArr->GetValue(3,7) << endl;
		
	}
	
}


void processGeodesicOptions(pi::Options& opts) {
	opts.addOption("-computeGeodesicDistance", "Compute geodesic distance to the set of points", "-computeGeodesicDistance input-data -scalarName destinationScalar", SO_NONE);
	
	opts.addOption("-singleSourceGeodesicDistance", "Compute geodesic distance from a given point id", "-singleSourceGeodesicDistance input-vtp source-id output-vtp", SO_NONE);
	
	opts.addOption("-computeLocalTangentMap", SO_NONE);
	
	opts.addOption("-computeGeodesicPath", "Compute a geodesic path from a vertex to the other vertex", "-computeGeodesicPath inputPoly inputTangentMap", SO_NONE);
	
	opts.addOption("-sphereSource", SO_NONE);
	
	opts.addOption("-testGeodesic", SO_NONE);
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
		//        computeLocalTangentMap(data, ugrid);
		vio.writeFile(outputFile, ugrid);
	} else if (opts.GetBool("-computeGeodesicPath")) {
		string input1File = args[0];
		string input2File = args[1];
		//		string outputFile = args[2];
		
		vtkPolyData* data1 = vio.readFile(input1File);
		vtkUnstructuredGrid* data2 = vtkUnstructuredGrid::SafeDownCast(vio.readDataFile(input2File));
		
		//        computeGeodesicPath(data1, data2);
	} else if (opts.GetBool("-singleSourceGeodesicDistance")) {
		string inputFile = args[0];
		int sourceId = atoi(args[1].c_str());
		string outputFile = args[2];
		
		vtkPolyData* inputPoly = vio.readFile(inputFile);
		vtkPolyData* outputPoly = computeSingleSourceGeodesicDistance(inputPoly, sourceId);
		vio.writeFile(outputFile, outputPoly);
	}else if (opts.GetBool("-sphereSource")) {
		vtkNew<vtkSphereSource> sphere;
		sphere->SetPhiResolution(8);
		sphere->SetThetaResolution(8);
		sphere->LatLongTessellationOff();
		sphere->Update();
		vio.writeFile("spheretest.vtp", sphere->GetOutput());
	} else if (opts.GetBool("-testGeodesic")) {
		testGeodesic();
	}
}