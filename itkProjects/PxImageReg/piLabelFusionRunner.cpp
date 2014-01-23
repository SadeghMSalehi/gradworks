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
#include "piImageProc.h"

#include <itkLabelOverlapMeasuresImageFilter.h>
#include <itkLabelStatisticsImageFilter.h>
#include <itkVectorMagnitudeImageFilter.h>
#include <itkVectorImageToImageAdaptor.h>

using namespace std;
using namespace libconfig;


namespace pi {

    static PxTools __tools;
    static ImageIO<LabelImage> __labeIO;

    void executeLabelFusionRunner(Options& opts, StringVector& args) {
        /// ---
        /// Label Fusion Format
        /// <code><br/>
        ///     imageSet = {<br/>
        ///         labels = ( label1, label2, label3, ... )<br/>
        ///         images = ( image1, image2, image3, ... )<br/>
        ///     }<br/>
        /// </code>
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

        /// - Read config file
        /// - If the file isn't read, exceptions will be thrown
        string configFile = opts.GetString("--fusion");
        cout << "reading config file ..." << configFile << endl;

        libconfig::Config config;
        try {
            config.readFile(configFile.c_str());
        } catch (libconfig::ParseException& e) {
            cout << e.getFile() << ": Line " << e.getLine() << "; " << e.getError()<< endl;
            return;
        }

        /// - Load images for fusion
        Setting& imageSet = config.lookup("imageSet");

        /// - Read all the labels and images into memory
        LabelImageVector labels;
        __tools.readImages<LabelImage>(labels, imageSet["labels"]);

        RealImageVector images;
        __tools.readImages<RealImage>(images, imageSet["images"]);

        /// - Perform label fusion and write out the result
        LabelImage::Pointer finalOutput = performLabelFusion(labels, images, targetImage);
        __labeIO.WriteImage(args[0], finalOutput);
    }




    LabelImage::Pointer performLabelFusion(LabelImageVector& labels, RealImageVector& images, RealImage::Pointer targetImage) {
        // prepare the output as an empty image with the same dimension and size
        LabelImage::Pointer outputImage = __labeIO.CopyImage(labels[0]);
        outputImage->FillBuffer(0);

        bool simpleMajorityVoting = true;
        // the simplest majority voting
        if (simpleMajorityVoting) {
            const int m = labels.size();
            const int n = labels[0]->GetPixelContainer()->Size();

            cout << "Processing " << n << " voxels ..." << endl;
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


    void executeVolumeOverlaps(Options& opts, StringVector& args) {
        typedef itk::LabelOverlapMeasuresImageFilter<LabelImage> MeasureFilterType;
        typedef itk::LabelStatisticsImageFilter<LabelImage, LabelImage> LabelFilterType;

        /// This can measure the Dice or the Jaccard coefficient.
        bool usingDice = opts.GetString("--overlap") == "dice";


        /// Set the output text file from the 0-th argument
        string outputFile = args[0];
        ofstream output(outputFile);

        /// Compute the volume overlap ratio using itk::LabelOverlapMeasuresImageFilter
        ImageIO<LabelImage> io;
        string refFile = args[1];

        /// If the reference image is not loaded, exit.
        LabelImage::Pointer refImage = io.ReadImage(refFile);
        if (refImage.IsNull()) {
            cout << "Can't read " << refFile << endl;
            exit(0);
        }

        /// - Enumerate valid labels
        LabelFilterType::Pointer labelFilter = LabelFilterType::New();
        labelFilter->SetInput(refImage);
        labelFilter->SetLabelInput(refImage);
        labelFilter->Update();
        LabelFilterType::ValidLabelValuesContainerType labels = labelFilter->GetValidLabelValues();

        /// - First, print valid labels in the reference image separated by a tab character
        cout << "# of labels: " << labels.size() << endl;
        output << refFile << "\t";
        for (int j = 0; j < labels.size(); j++) {
            output << (int) labels[j] << "\t";
        }
        output << endl;

        MeasureFilterType::Pointer filter = MeasureFilterType::New();
        filter->SetInput(0, refImage);

        for (int i = 2; i < args.size(); i++) {
            /// - Load a compared label image
            string compFile = args[i];
            LabelImage::Pointer compImage = io.ReadImage(compFile);
            filter->SetInput(1, compImage);
            filter->Update();

            std::vector<double> measures;
            measures.resize(labels.size());

            for (uint j = 0; j < labels.size(); j++) {
                /// - Compute the measurement for each label(Dice or Jaccard)
                measures[j] = 0;
                if (usingDice) {
                    measures[j] = filter->GetDiceCoefficient(labels[j]);
                } else {
                    measures[j] = filter->GetJaccardCoefficient(labels[j]);
                }
            }
            /// - Print the measurement in the order of valid labels
            output << compFile << "\t";
            for (int j = 0; j < labels.size(); j++) {
                output << measures[j] << "\t";
            }
            output << endl;
        }
    }


    void executeComputeDistanceMap(Options& parser, StringVector &args) {
        ImageIO<LabelImage> io;
        LabelImage::Pointer image = io.ReadCastedImage(args[0]);
        VectorImage::Pointer distanceMap = ComputeDistanceMap(image, args[2]);

        ImageIO<VectorImage> vio;
        vio.WriteImage(args[1], distanceMap);

        typedef itk::VectorMagnitudeImageFilter<VectorImage, RealImage> MagnitudeFilterType;
        MagnitudeFilterType::Pointer magFilter = MagnitudeFilterType::New();
        magFilter->SetInput(distanceMap);
        magFilter->Update();
        ImageIO<RealImage> rio;
        RealImage::Pointer distanceMagnitude = magFilter->GetOutput();
        rio.WriteImage("distmag.nii.gz", distanceMagnitude);


        fordim (k) {
            RealImage::Pointer outputImage = rio.NewImage(distanceMagnitude);
            outputImage->FillBuffer(0);

            ImageReal* outputImagePointer = outputImage->GetBufferPointer();
            VectorImage::PixelType* distanceMapPointer = distanceMap->GetBufferPointer();
            const int nPixels = outputImage->GetPixelContainer()->Size();
            for (int i = 0; i < nPixels; i++) {
                outputImagePointer[i] = distanceMapPointer[i][k];
            }

            if (k == 0) {
                rio.WriteImage("x.nii.gz", outputImage);
            } else if (k == 1) {
                rio.WriteImage("y.nii.gz", outputImage);
            } else if (k == 2) {
                rio.WriteImage("z.nii.gz", outputImage);
            }
        }

    }
}