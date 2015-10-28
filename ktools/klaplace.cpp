#include "klaplace.h"

#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>

#include "vtkio.h"

#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkDataSet.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkGenericCell.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkThresholdPoints.h>
#include <vtkNew.h>
#include <vtkExtractGrid.h>
#include <vtkGradientFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkCleanPolyData.h>
#include <vtkMath.h>
#include <vtkCellLocator.h>
#include <vtkPointLocator.h>
#include <vtkModifiedBSPTree.h>

#include "kgeometry.h"
#include "kstreamtracer.h"
#include "kvolume.h"

using namespace std;
using namespace pi;


static vtkIO vio;



//void processVTKUtils(pi::Options opts, pi::StringVector args) {
//    vtkIO vio;
//    string input1File, input2File, outputFile;
//    if (opts.GetBool("-sampleSurfaceScalars")) {
//        input1File = args[0];
//        input2File = args[1];
//        outputFile = args[2];
//        vtkDataSet* data1 = vio.readDataFile(input1File);
//        vtkPolyData* surf1 = vio.readFile(input2File);
//        sampleSurfaceScalar(data1, opts.GetString("-scalarName", "BorderPoints").c_str(), surf1, "labels");
//        vio.writeFile(outputFile, data1);
//        if (opts.GetString("-o") != "") {
//            vio.writeFile(opts.GetString("-o"), surf1);
//        }
//    } else if (opts.GetBool("-buildAdjacencyList")) {
//        input1File = args[0];
//        outputFile = args[1];
//
//        string scalarName = opts.GetString("-scalarName", "SampledSurfaceScalars");
//        vtkDataSet* data1 = vio.readDataFile(input1File);
//
//        vector<pair<vtkIdType, vector<vtkIdType> > > graph;
//        buildAdjacencyList(data1, scalarName, graph);
//
//        vio.writeFile(outputFile, data1);
//    } else  if (opts.GetBool("-thresholdStream")) {
//        runStreamLineThreshold(opts, args);
//    } else if (opts.GetBool("-rescaleStream")) {
//        runRescaleStream(opts, args);
//    } else if (opts.GetBool("-traceClipping")) {
//        runTraceClipping(opts, args);
//    } else if (opts.GetBool("-traceScalarCombine")) {
//        runTraceScalarCombine(opts, args);
//    } else if (opts.GetBool("-extractSurfaceBorder")) {
//        runExtractSurfaceBorder(opts, args);
//    }
//}
//

int main(int argc, char* argv[]) {
    Options opts;
    opts.addOption("-h", "print help message", SO_NONE);

    processVolumeOptions(opts);

    StringVector args = opts.ParseOptions(argc, argv, NULL);

    if (opts.GetBool("-h")) {
        cout << "## *kmesh* Usage" << endl;
        opts.PrintUsage();
        return 0;
    }

    processVolumeCommands(opts, args);

    return 0;
}