//
//  airCLI.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/1/13.
//
//

#ifndef __ParticleGuidedRegistration__airCLI__
#define __ParticleGuidedRegistration__airCLI__

#include <iostream>

#include "piOptions.h"
#include "piImageIO.h"

namespace air {
    typedef itk::Image<double,3> Image;
    typedef itk::Image<unsigned char,3> Label;
    typedef itk::Image<double,2> ImageSlice;
    typedef itk::Image<unsigned char,2> LabelSlice;

    extern pi::ImageIO<Image> __imageIO;
    extern pi::ImageIO<Label> __labelIO;

    class CommandLineTools {
    public:
        void ExtractSlice(Image::Pointer img, pi::SliceDirectionEnum dir, pi::IntVector range, std::string outputPattern);
        void PasteSlice(pi::StringVector input, std::string output);
        void PasteLabel(pi::StringVector input, std::string output);
        
        int Run(pi::Options* parser, pi::StringVector args);
    };
    
}

#endif /* defined(__ParticleGuidedRegistration__airCLI__) */
