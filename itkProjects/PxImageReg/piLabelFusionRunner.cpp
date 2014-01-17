//
//  piLabelFusionRunner.cpp
//  PxImageReg
//
//  Created by Joohwi Lee on 1/16/14.
//
//

#include "piLabelFusionRunner.h"

#include <iostream>

#include "piOptions.h"
#include "piTools.h"
#include "libconfig.h++"

using namespace std;
using namespace libconfig;


namespace pi {

    static PxTools __tools;
    static ImageIO<LabelImage> __labeIO;

    /**
     * SPIE 2014 experiments for label fusion
     *
     * config file format:
     *
     *  imageSet = {
     *      labels = ( image1, image2, image3, ... )
     *      images = ( image1, image2, image3, ... )
     *  }
     */
    void executeLabelFusionRunner(Options& opts, StringVector& args) {
        // check if the output is given
        if (args.size() == 0) {
            cout << "--fusion [config-file] output-file (target-image)" << endl;
            return;
        }

        RealImage::Pointer targetImage;
        if (args.size() > 1) {
            cout << "reading the target image " << args[1] << endl;
            ImageIO<RealImage> io;
            targetImage = io.ReadImage(args[1]);
        }

        // read config file
        // if the file isn't read, exceptions will be thrown
        string configFile = opts.GetString("--fusion");
        cout << "reading config file ..." << configFile << endl;

        // read the configuration file
        libconfig::Config config;
        try {
            config.readFile(configFile.c_str());
        } catch (libconfig::ParseException& e) {
            cout << e.getFile() << ": Line " << e.getLine() << "; " << e.getError()<< endl;
            return;
        }

        // load images for fusion
        Setting& imageSet = config.lookup("imageSet");

        // read all the labels and images into memory
        LabelImageVector labels;
        __tools.readImages<LabelImage>(labels, imageSet["labels"]);

        RealImageVector images;
        __tools.readImages<RealImage>(images, imageSet["images"]);

        // perform label fusion and write out the result
        LabelImage::Pointer finalOutput = performLabelFusion(labels, images, targetImage);
        __labeIO.WriteImage(args[0], finalOutput);
    }



    /**
     * Perform Label Fusion with Majority Voting
     *
     */
    LabelImage::Pointer performLabelFusion(LabelImageVector& labels, RealImageVector& images, RealImage::Pointer targetImage) {
        // prepare the output as an empty image with the same dimension and size
        LabelImage::Pointer outputImage = __labeIO.CopyImage(labels[0]);
        outputImage->FillBuffer(0);

        bool simpleMajorityVoting = true;
        // the simplest majority voting
        if (simpleMajorityVoting) {
            const int m = labels.size();
            const int n = labels[0]->GetPixelContainer()->Size();
            // for each pixel
            for (int i = 0; i < n; i++) {
                // for each image
                static const int UNKNOWN = 255;
                int major = UNKNOWN;
                int count = 0;
                for (int j = 0; j < m; j++) {
                    // read a label from an image
                    int label = labels[j]->GetBufferPointer()[i];
                    if (major == UNKNOWN) {
                        // when the major is unknown
                        major = label;
                        count = 1;
                    } else {
                        // increment the major count
                        if (label == major) {
                            count ++;
                        } else if (count > 0) {
                            count --;
                            if (count == 0) {
                                major = UNKNOWN;
                            }
                        }
                    }
                }
                outputImage->GetBufferPointer()[i] = major;
            }
        }

        return outputImage;
    }
}