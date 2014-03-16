//
//  piImageProfile.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/11/13.
//
//

#ifndef __ParticlesGUI__piImageProfile__
#define __ParticlesGUI__piImageProfile__

#include <iostream>

#include "myImageDef.h"
#include "itkLineIterator.h"


namespace pi {
    template<typename TInputImage>
    class ImageProfile {
    public:
        void SetInput(TInputImage::Pointer image);

    private:
        TInputImage::Pointer m_SrcImage;
    };
}
#endif /* defined(__ParticlesGUI__piImageProfile__) */
