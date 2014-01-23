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
    /// A utiliity class to perform I/O
    class PxTools {
    public:
        bool checkFile(std::string filename) {
            std::ifstream i(filename.c_str());
            bool file = i.is_open();
            i.close();
            return file;
        }

        /// @brief Read images given in files
        /// @param data A vector will contain loaded images (Output)
        /// @param files An instance of libconfig::Setting that stores a list of files
        template<class T>
        void readImages(std::vector<typename T::Pointer>& data, libconfig::Setting& files) {
            for (int i = 0; i < files.getLength(); i++) {
                std::string in = files[i];
                if (!checkFile(in)) {
                    cout << "can't read file " << in << endl;
                    exit(0);
                }
                cout << "Reading ... " << in << endl;
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

