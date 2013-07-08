//
//  peval.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/22/13.
//
//

#include "piImageProcessing.h"
#include "iostream"
#include "piImageIO.h"

using namespace std;
using namespace pi;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << argv[0] << " requires input-label1 input-label2" << endl;
        return 0;
    }
    pi::ImageIO<LabelImage> io;
    AtlasSimilarityScore score;
    score.Compute(io.ReadCastedImage(argv[1]), io.ReadCastedImage(argv[2]));
    cout << score;
}