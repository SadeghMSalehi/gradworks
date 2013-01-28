#include "myParticleCore.h"

int main(int argc, char* argv[]) {
    pi::ParticleSystem system;
//    my::StringVector labels;
//    labels.push_back("/data/Particles/image_square_label.nrrd");
//    labels.push_back("/data/Particles/image_circle_label.nrrd");
//    system.LoadShapes(labels);
//    system.SaveStatus("/tmpfs/output.txt");
    system.LoadStatus("/tmpfs/output.txt", 1);
    system.RunPreprocessing("/tmpfs/preprocessing.txt");
    system.LoadPreprocessing("/tmpfs/preprocessing.txt");
		system.Run();
    system.SaveStatus("/tmpfs/output2.txt");
    return 0;
}
