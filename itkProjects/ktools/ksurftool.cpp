#include "vtkio.h"
#include "piOptions.h"
#include "ksurftool.h"



using namespace std;
using namespace pi;


struct vtkProcess {
    void extractBoundaryPoints(vtkPoints* points, vtkDataArray* scalars, vtkPolyData* surf) {
        
    }
};


int main(int argc, char* argv[]) {
    Options opts;
    StringVector args = opts.ParseOptions(argc, argv, NULL);
    
}
