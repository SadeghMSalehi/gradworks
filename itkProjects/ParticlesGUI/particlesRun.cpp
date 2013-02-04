#include "myParticleCore.h"
#include "myParticleBSpline.h"

using namespace std;

int main(int argc, char* argv[]) {
    pi::ParticleSystem system;
    if (argc < 2) {
        return 0;
    }
    system.LoadSystem(argv[1]);

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


    system.RunPreprocessing();
    if (argc > 3) {
        system.SaveSystem(argv[3]);
    }

    // if the second argument is equal to '0', don't run main procedure
    if (argc > 2 && 0 != strcmp(argv[2], "0")) {
        system.Run();
        system.SaveSystem(argv[2]);
    }

    return 0;
}
