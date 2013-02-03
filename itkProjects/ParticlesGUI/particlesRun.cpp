#include "myParticleCore.h"

using namespace std;

int main(int argc, char* argv[]) {
    pi::ParticleSystem system;
    if (argc < 3) {
        return 0;
    }
    system.LoadSystem(argv[1]);
    system.RunPreprocessing();
    if (argc > 3) {
        system.SaveSystem(argv[3]);
    }

    // if the second argument is equalt to '0', don't run
    if (0 != strcmp(argv[2], "0")) {
        system.Run();
        system.SaveSystem(argv[2]);
    }
    return 0;
}
