#include "myParticleCore.h"

using namespace std;

int main(int argc, char* argv[]) {
    pi::ParticleSystem system;
//    my::StringVector labels;
//    labels.push_back("/data/Particles/image_square_label.nrrd");
//    labels.push_back("/data/Particles/image_circle_label.nrrd");
//    system.LoadShapes(labels);
//    system.SaveSystem("/tmpfs/output.txt");
    if (argc < 3) {
        return 0;
    }
    system.LoadSystem(argv[1], 1);
    system.RunPreprocessing("/data/Particles/CircleSquares/Preprocessing.txt");
    system.LoadPreprocessing("/data/Particles/CircleSquares/Preprocessing.txt");

	cout << "Done Load Preprocessing .." << endl;
    system.Run();
    system.SaveSystem(argv[2]);
    return 0;
}
