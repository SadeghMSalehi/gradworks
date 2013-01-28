#include "myParticleCore.h"

using namespace std;

int main(int argc, char* argv[]) {
    pi::ParticleSystem system;
//    my::StringVector labels;
//    labels.push_back("/data/Particles/image_square_label.nrrd");
//    labels.push_back("/data/Particles/image_circle_label.nrrd");
//    system.LoadShapes(labels);
//    system.SaveSystem("/tmpfs/output.txt");
    system.LoadSystem("/data/Particles/Output/output_1_0843.txt", 1);
//    system.RunPreprocessing("/tmpfs/preprocessing_2.txt");
//    system.LoadPreprocessing("/tmpfs/preprocessing_2.txt");
//    system.SaveSystem("/tmpfs/2.txt");

	cout << "Done Load Preprocessing .." << endl;
    system.Run();
    system.SaveSystem("/data/Particles/Output/ensemble_rodent_0016_2.txt");
    return 0;
}
