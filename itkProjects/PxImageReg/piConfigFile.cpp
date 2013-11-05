//
//  ConfigFile.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/24/13.
//
//

#include "piConfigFile.h"

using namespace std;
using namespace libconfig;

namespace pi {
    ConfigFile::ConfigFile() {

    }

    ConfigFile::ConfigFile(std::string file) {
        load(file);
    }

    bool ConfigFile::getBool(std::string path) {
        if (_config.exists(path)) {
            return _config.lookup(path);
        }
        return false;
    }

    RealImage::RegionType ConfigFile::offsetRegion(RealImage::RegionType region, std::string path) {
        if (_config.exists(path)) {
            RealImage::IndexType lowerIdx = region.GetIndex();
            RealImage::IndexType upperIdx = region.GetUpperIndex();
            Setting& config = _config.lookup(path);
            for (int i = 0; i < lowerIdx.GetIndexDimension(); i++) {
                lowerIdx[i] += (int) config[i];
                upperIdx[i] += (int) config[i + lowerIdx.GetIndexDimension()];
            }
            region.SetIndex(lowerIdx);
            region.SetUpperIndex(upperIdx);
        }
        return region;
    }


    void ConfigFile::load(std::string file) {
        _config.readFile(file.c_str());

        string useImages = "images";
        string useLabels = "labels";

        if (_config.lookupValue("use-images", useImages)) {
            readImages<RealImage>(_images, useImages);
        }

        if (_config.lookupValue("use-labels", useLabels)) {
            readImages<LabelImage>(_labels, useLabels);
        }
    }
}