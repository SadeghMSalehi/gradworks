#include "piParticleCore.h"
#include "piParticleBSpline.h"
#include "piParticleSystemSolver.h"
#include "piOptions.h"

using namespace std;
using namespace pi;

int main(int argc, char* argv[]) {
    Options parser;
    parser.ParseOptions(argc, argv, NULL);
    
    StringVector& args = parser.GetStringVector("args");
    if (args.size() < 2) {
        return 0;
    }
    
    ParticleSystemSolver solver;
    solver.LoadConfig(args[0].c_str());
    solver.Preprocessing();
    solver.Run();
    solver.SaveConfig(args[1].c_str());
    
//    if (!system.LoadSystem(argv[1])) {
//        cout << "Can't open " << argv[1] << endl;
//        return 0;
//    }

//    {
//        pi::ParticleBSpline bspline;
//        bspline.SetReferenceImage(system.GetImageContext().GetLabel(1));
//        bspline.EstimateTransform(system[1], system[0]);
//        pi::FieldTransformType::Pointer transform = bspline.GetTransform();
//
//        pi::FieldTransformType::InputPointType point;
//        fordim(k) {
//            point[k] = system[1][0].x[k];
//        }
//        pi::FieldTransformType::OutputPointType outPoint = transform->TransformPoint(point);
//        cout << system[1][0].x[0] << "," << system[1][0].x[1] << "; " << system[0][0].x[0] << "," << system[0][0].x[1] << " ~= " << outPoint << endl;
//    }
//
//    {
//        pi::ParticleBSpline bspline;
//        bspline.SetReferenceImage(system.GetImageContext().GetLabel(1));
//        bspline.EstimateTransform(system[0], system[1]);
//        pi::FieldTransformType::Pointer transform = bspline.GetTransform();
//
//        pi::FieldTransformType::InputPointType point;
//        fordim(k) {
//            point[k] = system[0][0].x[k];
//        }
//        pi::FieldTransformType::OutputPointType outPoint = transform->TransformPoint(point);
//        cout << system[1][0].x[0] << "," << system[1][0].x[1] << "; " << system[0][0].x[0] << "," << system[0][0].x[1] << " ~= " << outPoint << endl;
//    }

//    if (argc > 2 && 0 == strcmp(argv[2], "1")) {
//        system.RunPreprocessing(true, 0.1);
//    } else if (argc > 2 && 0 == strcmp(argv[2], "0")) {
//        system.RunPreprocessing(false);
//    }
//    if (argc > 3) {
//        system.SaveSystem(argv[3]);
//    }
//
//    // if the second argument is equal to '0', don't run main procedure
//    if (argc > 2 && 0 != strcmp(argv[2], "0")) {
//        system.Run();
//        system.SaveSystem(argv[2]);
//    }

    return 0;
}
