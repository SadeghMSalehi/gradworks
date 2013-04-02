//
//  kextract.cpp
//  ktools
//
//  Created by Joohwi Lee on 4/1/13.
//
//

#include "kextract.h"
#include "piImageIO.h"
#include "piOptions.h"

using namespace pi;

int main(int argc, char* argv[]) {
    // Option declaration
    CSimpleOpt::SOption specs[] = {
        { 0, "-i", SO_REQ_SEP },
        { 1, "-o", SO_REQ_SEP },
        SO_END_OF_OPTIONS
    };

    // Option parsing
    Options parser;
    StringVector args = parser.ParseOptions(argc, argv, specs);
    if (args.size() < 2) {
        cout << "kextract [input-image] [output-image-pattern: slice-%1.nii.gz]" << endl;
        return 0;
    }

}