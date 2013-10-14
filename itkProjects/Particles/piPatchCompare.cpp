//
//  piPatchCompare.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/13/13.
//
//

#include "piPatchCompare.h"
#include "piImageIO.h"

#include <vector>
#include <algorithm>

namespace pi {
    struct PatchSimilarity {
        int originAtlas;
        int originLabel;
        double metric;

        PatchSimilarity(int o, int l, double s) {
            originAtlas = o;
            originLabel = l;
            metric = s;
        }

        PatchSimilarity(const PatchSimilarity& o) {
            originAtlas = o.originAtlas;
            originLabel = o.originLabel;
            metric = o.metric;
        }

        PatchSimilarity& operator=(const PatchSimilarity& o) {
            originAtlas = o.originAtlas;
            originLabel = o.originLabel;
            metric = o.metric;
            return (*this);
        }
    };

    bool patch_compare_descending(const PatchSimilarity& a, const PatchSimilarity& b) {
        return (a.metric < b.metric);
    }

    void PatchCompare::setAtlasImages(PatchImageVector atlasImages) {
        _atlasImages = atlasImages;
    }

    void PatchCompare::setAtlasLabels(LabelImageVector atlasLabels) {
        _atlasLabels = atlasLabels;
    }

    void PatchCompare::setTargetImage(RealImage::Pointer target) {
        _targetImage = target;
    }

    LabelImage::Pointer PatchCompare::getTargetLabel() {
        return _targetLabel;
    }

    void PatchCompare::setTargetROI(LabelImage::Pointer target) {
        _targetROI = target;
    }

    void PatchCompare::setParticleSystem(pi::ParticleSystem *system) {
        _system = system;
    }

    void PatchCompare::setTargetRadius(int r) {
        _targetRadius = r;
    }

    // assume isotropic patch size
    void PatchCompare::setPatchRadius(int r) {
        _patchRadius = r;
    }

    // run estimation comparing patches
    void PatchCompare::estimateLabel() {
        RealImage::RegionType targetRegion = _targetImage->GetBufferedRegion();
        PatchImage::Pointer targetPatch = buildPatchImage((_targetImage));

        ImageIO<LabelImage> labelIO;
        _targetLabel = labelIO.CastImageFromS<RealImage>(_targetImage);
        _targetLabel->FillBuffer(0);

        itk::ImageRegionConstIteratorWithIndex<PatchImage> iter(targetPatch, targetRegion);
        // for each particle
        // for each voxel v1 inside the target radius
        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            RealImage::IndexType idx = iter.GetIndex();

            bool insideROI = _targetROI->GetPixel(idx) > 0;
            if (!insideROI) {
                continue;
            }

            // extract patch P1 at v1
            LocalPatch targetPatch = iter.Get();

            // for each voxel v2 inside the search radius of v1

            // define search area
            int searchRadius = 3;

            RealImage::RegionType searchRegion;
            RealImage::IndexType searchIndex = idx;
            RealImage::SizeType searchSize;
            searchSize.Fill(searchRadius);
            RealImage::OffsetType searchOffset;
            searchOffset.Fill(searchRadius/2);
            searchIndex -= searchOffset;
            searchRegion.SetIndex(idx);
            searchRegion.SetSize(searchSize);

            // for each atlas
            std::vector<PatchSimilarity> similarityRanking;
            for (unsigned int i = 0; i < _atlasImages.size(); i++) {
                itk::ImageRegionConstIteratorWithIndex<PatchImage> searchIter(_atlasImages[i], searchRegion);
                for (searchIter.GoToBegin(); !searchIter.IsAtEnd(); ++searchIter) {
                    // extract patch P2 at v2
                    RealImage::IndexType originIdx = searchIter.GetIndex();
                    LocalPatch atlasPatch = searchIter.Get();

                    // compute patch distance between P1 and P2
                    double ssd = 0;
                    for (unsigned int j = 0; j < targetPatch.Size(); j++) {
                        double diff = atlasPatch[j] - targetPatch[j];
                        ssd += (diff * diff);
                    }

                    int originLabel = _atlasLabels[i]->GetPixel(originIdx);
                    similarityRanking.push_back(PatchSimilarity(i, originLabel, ssd));
                }
                // end for
            }
            // end for

            // choose top k nearest patches
            std::sort(similarityRanking.begin(), similarityRanking.end(), patch_compare_descending);

//            cout << similarityRanking[0].metric << ", " << similarityRanking[1].metric << ", " << similarityRanking[2].metric << endl;

            // estimate label for v1
            std::vector<int> counter(10);
            std::fill(counter.begin(), counter.end(), 0);
            for (unsigned int k = 0; k < 5; k++) {
                counter[similarityRanking[k].originLabel] ++;
            }
            int estimatedLabel = std::distance(counter.begin(), std::max_element(counter.begin(), counter.end()));
            _targetLabel->SetPixel(idx, estimatedLabel);
        }
        // end for
    }

    PatchImage::Pointer PatchCompare::buildPatchImage(RealImage::Pointer image) {
        const int patchWidth = PATCH_SIZE;
        ImageIO<PatchImage> io;

        RealImage::RegionType region = image->GetBufferedRegion();
        RealImage::SizeType imageSize = region.GetSize();

        PatchImage::Pointer outputImage = io.NewImageT(imageSize);
        LocalPatch zeroVector;
        zeroVector.Fill(0);
        outputImage->FillBuffer(zeroVector);

        // set region to contain full patch
        RealImage::IndexType regionIdx;
        regionIdx[0] = patchWidth/2 + 1;
        regionIdx[1] = patchWidth/2 + 1;

        RealImage::SizeType regionSize;
        regionSize[0] = imageSize[0] - patchWidth;
        regionSize[1] = imageSize[1] - patchWidth;

        region.SetIndex(regionIdx);
        region.SetSize(regionSize);

        // iterate over the region
        RealImageIteratorType iter(image, region);
        RealImage::RegionType patchRegion;
        RealImage::OffsetType patchOffset;
        patchOffset.Fill(patchWidth/2);
        RealImage::SizeType patchSize;
        patchSize.Fill(patchWidth);

        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            LocalPatch patch;
            // extract patch
            RealImage::IndexType idx = iter.GetIndex();
            idx -= patchOffset;
            patchRegion.SetIndex(idx);
            patchRegion.SetSize(patchSize);
            RealImageIteratorType patchIter(image, patchRegion);
            int i = 0;
            for (patchIter.GoToBegin(); !patchIter.IsAtEnd(); ++patchIter) {
                patch[i] = patchIter.Get();
                i++;
            }
            outputImage->SetPixel(idx, patch);
        }
        return outputImage;
    }


    void transferLabelsWithPatch(StringVector& args, std::string output) {
        int nAtlas = (args.size() - 2) / 2;
        if (args.size() == 0 || args.size() != (nAtlas+1)*2 || nAtlas < 1) {
            cout << "LabelTransfer: target-image target-roi atlas-patch-#1 atlas-patch-#2 ... atlas-label-#1 atlas-label-#2 ..." << endl;
        }

        ImageIO<RealImage> realIO;
        ImageIO<LabelImage> labelIO;
        ImageIO<PatchImage> patchIO;

        RealImage::Pointer targetImage = realIO.ReadCastedImage(args[0]);
        LabelImage::Pointer targetROI = labelIO.ReadImage(args[1]);

        cout << "# of atlases: " << nAtlas << endl;
        PatchImageVector atlasPatches;
        for (int i = 2; i < nAtlas + 2; i++) {
            PatchImage::Pointer patchAtlas = patchIO.ReadImage(args[i]);
            cout << "Reading patch #" << (i-1) << endl;
            atlasPatches.push_back(patchAtlas);
        }

        LabelImageVector atlasLabels;
        for (int i = nAtlas + 2; i < nAtlas * 2 + 2; i++) {
            LabelImage::Pointer labelAtlas = labelIO.ReadImage(args[i]);
            cout << "Reading label #" << (i-nAtlas-1) << endl;
            atlasLabels.push_back(labelAtlas);
        }

        PatchCompare patchRun;
        patchRun.setAtlasImages(atlasPatches);
        patchRun.setAtlasLabels(atlasLabels);
        patchRun.setTargetImage(targetImage);
        patchRun.setTargetROI(targetROI);
        patchRun.estimateLabel();

        labelIO.WriteImage(output, patchRun.getTargetLabel());
    }
}