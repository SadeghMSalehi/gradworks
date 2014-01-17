//
//  piTools.h
//  PxImageReg
//
//  Created by Joohwi Lee on 1/16/14.
//
//

#ifndef __PxImageReg__piTools__
#define __PxImageReg__piTools__

#include <iostream>
#include <fstream>
#include <vector>
#include "piImageIO.h"

#include "libconfig.h++"

namespace pi {
    class PxTools {
    public:
        bool checkFile(std::string filename) {
            std::ifstream i(filename.c_str());
            bool file = i.is_open();
            i.close();
            return file;
        }

        template<class T>
        void readImages(std::vector<typename T::Pointer>& data, libconfig::Setting& files) {
            for (int i = 0; i < files.getLength(); i++) {
                std::string in = files[i];
                if (!checkFile(in)) {
                    cout << "can't read file " << in << endl;
                    exit(0);
                }
                ImageIO<T> io;
                data.push_back(io.ReadCastedImage(in));
            }
        }

        template<class T>
        void writeImages(std::vector<typename T::Pointer>& data, libconfig::Setting& files) {
            for (int i = 0; i < files.getLength(); i++) {
                std::string in = files[i];
                ImageIO<T> io;
                io.WriteImage(in, data[i]);
            }
        }
    };
}


#endif /* defined(__PxImageReg__piTools__) */

