//
//  ktest1.cpp
//  ktools
//
//  Created by Joowhi Lee on 8/28/15.
//
//

#include "ktest1.h"
#include "kdatastructure.h"
#include "kgeometry.h"

#include <vtkNew.h>
#include <vtkLine.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#include <vtkMath.h>
#include <vtkTransform.h>
#include <vtkCellLocator.h>
#include <vtkPointLocator.h>
#include <vtkTriangle.h>

#include <queue>

#include "vtkio.h"



using namespace std;

static vtkIO vio;

template <class T>
void print(T* p, size_t n = 3) {
    for (size_t j = 0; j < n; j++) {
        cout << p[j] << " ";
    }
	cout << endl;
}

void testLineIntersection() {
    double pt1[3] = { 0, 0, 0 };
    double pt2[3] = { 0, 1, 0 };
    double xt1[3] = { 2, 2, 0 };
    double xt2[3] = { 1, 0, 0 };
    
    vtkNew<vtkLine> line;
    line->GetPoints()->SetPoint(0, pt1);
    line->GetPoints()->SetPoint(1, xt1);

    double x[3] = { 0, };
    double t = 0;
    double pcoords = 0;
    int subId = 0;
    int ret = line->IntersectWithLine(pt2, xt2, 0, t, x, &pcoords, subId);
    cout << ret << " " << t << " " << pcoords << " " << x[0] << "," << x[1] << "," << x[2] << endl;
}


void testPolyLine() {
    vtkNew<vtkPolyLine> lineSet;

    
}

void testFaceEdgeEnum() {
    vtkIO vio;
    vtkPolyData* poly = vio.readFile("sphere_sample_8x8.vtk");
    poly->BuildLinks();
    poly->BuildCells();
    
    vtkCellArray* polys = poly->GetPolys();
    polys->InitTraversal();
    
    vtkNew<vtkIdList> pts;
    while (polys->GetNextCell(pts.GetPointer())) {
        
    }
}

void testRotation() {
	vtkNew<vtkTransform> txfm;
	Geometry g;
    
    double u[3] = { 0, 0, 3 };
    double v[3] = { 1, 1, 3-sqrt(2) };
    double n[3] = { 0, 0, 1 };
    double tv[3] = { 0, };
    
    cout << g.tangentVector(u, v, n, tv) << endl;
    print(tv);
    cout << vtkMath::Distance2BetweenPoints(u, v) << " ~= " << vtkMath::Distance2BetweenPoints(u, tv);
	
//	for (size_t j = 0; j < 100000; j++) {
//		double p[3] = { 0, }, q[3] = { 0, }, r[3] = { 0, };
//		for (size_t k = 0; k < 2; k++) {
//			p[k] = (rand()-RAND_MAX/2)/double(RAND_MAX);
//			q[k] = r[k] = (rand()-RAND_MAX/2)/double(RAND_MAX);
//		}
//
//		vtkMath::MultiplyScalar(r, -1);
//        r[2] = 0;
//
//		txfm->Identity();
//		double cross[3];
//		double angDeg = g.rotateVector(p, r, txfm.GetPointer(), cross);
//		
//        double pr[3];
//        txfm->TransformPoint(p, pr);
//
//        vtkMath::Normalize(pr);
////        cout << "prN: "; print(pr);
//        
//        vtkMath::Normalize(r);
////        cout << "rN: "; print(r);
//        
//        vtkMath::Subtract(pr, r, r);
//        if (vtkMath::Norm(r) > 1e-10) {
//            cout << "Diff: " << vtkMath::Norm(r) << endl;
//        }
//        
//	}
}

struct NodeX {
	int id;
	int val;
	
	NodeX(int i, int v): id(i), val(v) {}
	
	bool operator<(const NodeX& o) const {
		return val > o.val;
	}
	
	bool operator==(const NodeX& o) const {
		return id == o.id;
	}
};

void testPQ() {
	priority_queue<NodeX> pq;
	pq.push(NodeX(2,3));
	pq.push(NodeX(4,3));
	pq.push(NodeX(2,4));
	pq.push(NodeX(2,1));
	cout << pq.size() << endl;
	
	while (!pq.empty()) {
		cout << pq.top().id << "/" << pq.top().val << endl;
		pq.pop();
		cout << pq.size() << endl;
	}
	
	
    typedef FibQueue<double> FibQueueD;
    typedef FibQueueD::FibNode FibNodeD;
    
    FibQueueD fq;
    
    vector<FibNodeD*> idMap(10);
    unordered_map<FibNodeD*, vtkIdType> nodeMap;
    for (size_t j = 0; j < 10; j++) {
        idMap[j] = fq.push(10 - j*.1, &idMap[j]);
    }
    
//    
//    while (!fq.empty()) {
//        FibNodeD* tn = fq.topNode();
////        vtkIdType id = (vtkIdType) (-&idMap.back() + (long int) (tn->payload));
//        cout << id << endl;
//        cout << nodeMap[tn] << "; " << tn->key << endl;
//        fq.pop();
//    }
//	
	
}

void testPlane() {
    vtkDataSet* ds = vio.readDataFile("312.laplaceSol.vtp");
    vtkPolyData* pd = vtkPolyData::SafeDownCast(ds);
    pd->BuildCells();
    pd->BuildLinks();
    double x[3] = { -16.1316, -18.7324, 10.5392 }, closestPoint[3] = {0,}, pcoords[3] = {0,}, weights[3] = {0,};
    int subId;
    vtkIdType cellId = pd->FindCell(x, NULL, -1, 10, subId, pcoords, weights);
    double dist2 = -1;
    double dist2Min = DBL_MAX;
    vtkIdType minDistCellId = -1;
    double minDistPcoords[3] = {0,};
    for (size_t j = 0; j < pd->GetNumberOfCells(); j++) {
        vtkCell* cell = pd->GetCell(j);
//        vtkTriangle* tri = vtkTriangle::SafeDownCast(cell);
        int ret = cell->EvaluatePosition(x, closestPoint, subId, pcoords, dist2, weights);
        if (dist2 < dist2Min) {
            dist2Min = dist2;
            minDistCellId = j;
            memcpy(minDistPcoords, pcoords, sizeof(minDistPcoords));
        }
//        cout << ret << "; " << dist2 << endl;
    }
    cout << cellId << endl;
    cout << minDistCellId << "; " << dist2Min << endl;
    print(minDistPcoords);
    
}

void processTest1Options(pi::Options& opts) {
    opts.addOption("-testLineIntersection", SO_NONE);
    opts.addOption("-testPolyLine", SO_NONE);
	opts.addOption("-testFaceEdgeEnum", SO_NONE);
	opts.addOption("-testRotation", SO_NONE);
	opts.addOption("-testPQ", SO_NONE);
	opts.addOption("-testPlane", SO_NONE);

}

void processTest1Commands(pi::Options& opts, pi::StringVector& args) {
    if (opts.GetBool("-testLineIntersection")) {
        testLineIntersection();
    } else if (opts.GetBool("-testPolyLine")) {
        testPolyLine();
    } else if (opts.GetBool("-testFaceEdgeEnum")) {
        testFaceEdgeEnum();
	} else if (opts.GetBool("-testRotation")) {
		testRotation();
	} else if (opts.GetBool("-testPQ")) {
		testPQ();
    } else if (opts.GetBool("-testPlane")) {
        testPlane();
    }
}