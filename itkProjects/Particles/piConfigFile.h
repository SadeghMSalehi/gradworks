//
//  ConfigFile.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/24/13.
//
//

#ifndef __ParticleGuidedRegistration__ConfigFile__
#define __ParticleGuidedRegistration__ConfigFile__

#include <iostream>
#include <vector>
#include <string>
#include <libconfig.h++>
#include "piImageIO.h"
#include "piImageDef.h"


namespace pi {
    class ConfigFile {
    public:
        ConfigFile();
        ConfigFile(std::string file);

        void load(std::string file);

        inline const int imageCount() { return _images.size(); }
        inline const int labelCount() { return _labels.size(); }

        inline RealImage::Pointer image(int i) {
            return _images[i];
        }

        inline LabelImage::Pointer label(int i) {
            return _labels[i];
        }

        bool checkFile(std::string filename) {
            ifstream i(filename.c_str());
            bool file = i.is_open();
            i.close();
            return file;
        }
        
        bool exists(std::string path) {
            return _config.exists(path);
        }

        inline libconfig::Setting& operator[](std::string path) {
            return _config.lookup(path);
        }

        bool getBool(std::string path);

        RealImage::RegionType offsetRegion(RealImage::RegionType region, std::string path);

        template<class T>
        void readImages(std::vector<typename T::Pointer>& data, std::string path) {
            libconfig::Setting& files = _config.lookup(path);
            for (int i = 0; i < files.getLength(); i++) {
                std::string in = files[i];
                ImageIO<T> io;
                data.push_back(io.ReadCastedImage(in));
            }
        }

        template<class T>
        void readImages(std::vector<typename T::Pointer>& data, libconfig::Setting& files) {
            for (int i = 0; i < files.getLength(); i++) {
                std::string in = files[i];
                ImageIO<T> io;
                data.push_back(io.ReadCastedImage(in));
            }
        }

        template<class T>
        void writeImages(std::vector<typename T::Pointer>& data, std::string path) {
            libconfig::Setting& files = _config.lookup(path);
            for (int i = 0; i < files.getLength(); i++) {
                std::string in = files[i];
                ImageIO<T> io;
                io.WriteImage(in, data[i]);
            }
        }

    private:
        std::vector<RealImage::Pointer> _images;
        std::vector<LabelImage::Pointer> _labels;

        libconfig::Config _config;
    };
}
#endif /* defined(__ParticleGuidedRegistration__ConfigFile__) */