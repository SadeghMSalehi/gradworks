//
//  peval.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/22/13.
//
//

#include "piImageProcessing.h"
#include "iostream"
#include "itkImageIO.h"

using namespace std;
using namespace pi;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << argv[0] << " requires input-label1 input-label2" << endl;
        return 0;
    }
    itkcmds::__noverbose = 1;
    itkcmds::itkImageIO<LabelImage> io;
    AtlasSimilarityScore score;
    score.Compute(io.ReadImageT(argv[1]), io.ReadImageT(argv[2]));
    cout << score;
}