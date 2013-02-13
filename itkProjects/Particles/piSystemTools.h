//
//  piSystemTools.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/13/13.
//
//

#ifndef __ParticlesGUI__piSystemTools__
#define __ParticlesGUI__piSystemTools__

#include <iostream>

namespace pi {
    class SystemTools {
        static inline bool file_exists(const char* filename) {
            ifstream ifile(filename);
            return ifile;
        }
    }
}
#endif /* defined(__ParticlesGUI__piSystemTools__) */
