//
//  piParticleTools.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/16/13.
//
//

#include "piParticleTools.h"
#include "piParticleSystemSolver.h"
#include "piImageIO.h"

namespace pi {
    static void end() {
        exit(EXIT_SUCCESS);
    }
    static void die() {
        exit(EXIT_FAILURE);
    }

    static ImageIO<RealImage> __io;

    #pragma mark -
    ParticleTools::ParticleTools(Options& o, StringVector& a): opts(o), args(a) {
        if (opts.GetBool("--coverlabel")) {
            runCoverLabel();
        }
    }

    void ParticleTools::runCoverLabel() {
        if (args.size() < 3) {
            cout << "--coverlabel [config-input] [image-input] [label-output]" << endl;
            die();
        }

        ParticleSystemSolver solver;
        solver.LoadConfig(args[0].c_str());


        // prepare input & output images
        RealImage::Pointer inputImage = __io.ReadCastedImage(args[1]);

        RealImage::Pointer outputRealImage = __io.CopyImage(inputImage);
        outputRealImage->FillBuffer(0);

        LabelImage::Pointer outputImage = __io.CastImageToS<LabelImage>(inputImage);
        outputImage->FillBuffer(0);


        // set patch size
        RealImage::RegionType srcRegion;
        RealImage::RegionType dstRegion;
        RealImage::SizeType patchSize;
        fordim (k) {
            patchSize[k] = ATTR_SIZE;
        }

        // loop over to copy patches
        ParticleSubject& srcSubj = solver.m_System[0];
        ParticleSubject& dstSubj = solver.m_System[1];
        for (int i = 0; i < srcSubj.size(); i++) {
            RealImage::PointType px;
            RealImage::PointType qx;
            fordim (k) {
                px[k] = srcSubj[i].x[k];
                qx[k] = dstSubj[i].x[k];
            }
            RealImage::IndexType srcIdx, dstIdx;
            inputImage->TransformPhysicalPointToIndex(px, srcIdx);
            inputImage->TransformPhysicalPointToIndex(qx, dstIdx);
            fordim (k) {
                srcIdx[k] -= (patchSize[k] / 2);
                dstIdx[k] -= (patchSize[k] / 2);
            }
            srcRegion.SetIndex(srcIdx);
            srcRegion.SetSize(patchSize);
            dstRegion.SetIndex(dstIdx);
            dstRegion.SetSize(patchSize);

            RealImageIteratorType inputIter(inputImage, srcRegion);
            LabelImageIteratorType outputIter(outputImage, dstRegion);
            RealImageIteratorType outputRealIter(outputRealImage, dstRegion);
            outputRealIter.GoToBegin();
            for (inputIter.GoToBegin(), outputIter.GoToBegin(); !inputIter.IsAtEnd() && !outputIter.IsAtEnd(); ++inputIter, ++outputIter) {
                outputIter.Set(outputIter.Get() + 1);
                outputRealIter.Set(outputRealIter.Get() + inputIter.Get());
                ++outputRealIter;
            }
        }

        int nVoxels = outputRealImage->GetPixelContainer()->Size();
        RealImage::PixelType* outputRealBuff = outputRealImage->GetBufferPointer();
        LabelImage::PixelType* outputBuff = outputImage->GetBufferPointer();
        for (int i = 0; i < nVoxels; i++) {
            outputRealBuff[i] /= outputBuff[i];
        }
        __io.WriteImage(args[3], outputRealImage);
        __io.WriteImageS<LabelImage>(args[2], outputImage);
        end();
    }
}