//
//  piImageDef.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/2/13.
//
//

#include <itkExtractImageFilter.h>
#include "piImageDef.h"


namespace pi {
    typedef itk::ExtractImageFilter<RealImage3, RealImage2> FilterType;

    // external variable for easy access
    RealImageTools __realImageTools;

    // split a volume image into slices
    RealImage2Vector RealImageTools::sliceVolume(RealImage3::Pointer volume, int dim) {
        RealImage2Vector slices;
        RealImage3::RegionType region = volume->GetBufferedRegion();
        RealImage3::SizeType sz = region.GetSize();
        const RealImage3::IndexType idx = region.GetIndex();
        
        const int nSlices = sz[dim];
        slices.resize(nSlices);
        region.SetSize(dim, 0);
        
        for (int i = idx[dim]; i < idx[dim] + nSlices; i++) {
            region.SetIndex(dim, i);
            
            FilterType::Pointer filter = FilterType::New();
            filter->SetInput(volume);
            filter->SetExtractionRegion(region);
            filter->SetDirectionCollapseToGuess();
            filter->Update();
            slices[i] = filter->GetOutput();
            slices[i]->DisconnectPipeline();
        }

        return slices;
    }
}