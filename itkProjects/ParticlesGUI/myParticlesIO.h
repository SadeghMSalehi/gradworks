//
//  myTimeSeries.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 1/16/13.
//
//

#ifndef __ParticlesGUI__myTimeSeries__
#define __ParticlesGUI__myTimeSeries__

#include <iostream>
#include "vnlCommon.h"

namespace my {
    class ParticlesContext {

    };

    class ParticlesIO {
    public:
        void OpenWriter(std::string filename);
        void WriteVector(std::string name, VNLVector& vector);
        void WriteMatrix(std::string name, VNLMatrix& mat);
        void CloseWriter(std::string filename);

        void Read(std::string filename);
    };
}

#endif /* defined(__ParticlesGUI__myTimeSeries__) */
