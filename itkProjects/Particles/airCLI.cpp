//
//  airCLI.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/1/13.
//
//

#include "airCLI.h"
#include "airImageAlgorithm.h"

#include <itkExtractImageFilter.h>
#include <itkGenerateImageSource.h>
#include <itkPasteImageFilter.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include <QString>

using namespace pi;
using namespace std;


namespace air {
    ImageIO<Image> __imageIO;
    ImageIO<Label> __labelIO;
    ImageIO<ImageSlice> __imageSliceIO;
    ImageIO<LabelSlice> __labelSliceIO;

    typedef itk::ImageRegionIteratorWithIndex<Image> ImageIterator;
    typedef itk::ImageRegionIteratorWithIndex<ImageSlice> ImageSliceIterator;
    typedef itk::ImageRegionIteratorWithIndex<Label> LabelIterator;
    typedef itk::ImageRegionIteratorWithIndex<LabelSlice> LabelSliceIterator;

    typedef air::ImageAlgorithm<Image,Label> ImageAlgo;

    void CommandLineTools::ExtractSlice(Image::Pointer img, pi::SliceDirectionEnum dir, pi::IntVector range, std::string outputPattern) {

        for (int i = range[0]; i <= range[1]; i++) {
            Image::RegionType extractRegion = img->GetBufferedRegion();
            if (i < 0 || i >= extractRegion.GetSize(dir)) {
                continue;
            }
            extractRegion.SetIndex(dir,i);
            extractRegion.SetSize(dir,0);
            typedef itk::ExtractImageFilter<Image, ImageSlice> Filter;
            Filter::Pointer filter = Filter::New();
            filter->SetInput(img);
            filter->SetExtractionRegion(extractRegion);
            filter->SetDirectionCollapseToGuess();
            filter->Update();
            QString filePattern = QString::fromStdString(outputPattern);
            QString fileOutput = filePattern.arg(i,3,10,QChar('0'));
            __imageSliceIO.WriteImage(fileOutput.toStdString(), filter->GetOutput());
        }

        return;
    }

    void CommandLineTools::PasteSlice(StringVector args, string output) {
        int nSlices = args.size() - 1;
        ImageInfo imageInfo;
        ImageSlice::Pointer firstImage = __imageSliceIO.ReadCastedImage(args[1], imageInfo);
        ImageSlice::RegionType firsrtRegion = firstImage->GetBufferedRegion();
        Image::Pointer pasteImage = __imageIO.NewImageT(firsrtRegion.GetSize(0), firsrtRegion.GetSize(1), nSlices);
        pasteImage->FillBuffer(0);

        ImageIterator pasteIter(pasteImage, pasteImage->GetBufferedRegion());
        pasteIter.GoToBegin();
        for (int i = 0; i < nSlices; i++) {
            cout << "Processing #" << i << ": " << args[i+1] << endl;

            ImageSlice::Pointer sliceImage = __imageSliceIO.ReadCastedImage(args[i+1]);
            if (firsrtRegion != sliceImage->GetBufferedRegion()) {
                cout << "Error Processing #" << i << ": " << args[i+1] << endl;
                return;
            }
            ImageSliceIterator sliceIter(sliceImage, sliceImage->GetBufferedRegion());
            sliceIter.GoToBegin();
            while (!sliceIter.IsAtEnd()) {
                pasteIter.Set(sliceIter.Get());
                ++pasteIter;
                ++sliceIter;
            }
        }

        __imageIO.WriteCastedImage(output, pasteImage, imageInfo.componenttype);
        return;
    }

    void CommandLineTools::PasteLabel(StringVector args, string output) {
        int nSlices = args.size() - 1;
        LabelSlice::Pointer firstImage = __labelSliceIO.ReadCastedImage(args[1]);
        LabelSlice::RegionType sliceRegion = firstImage->GetBufferedRegion();
        Label::Pointer pasteImage = __labelIO.NewImageT(sliceRegion.GetSize(0), sliceRegion.GetSize(1), nSlices);

        LabelIterator pasteIter(pasteImage, pasteImage->GetBufferedRegion());
        pasteIter.GoToBegin();
        for (int i = 0; i < nSlices; i++) {
            cout << "Processing #" << i << ": " << args[i+1] << endl;

            LabelSlice::Pointer sliceImage = __labelSliceIO.ReadCastedImage(args[i+1]);
            LabelSliceIterator sliceIter(sliceImage, sliceImage->GetBufferedRegion());
            sliceIter.GoToBegin();
            while (!sliceIter.IsAtEnd()) {
                pasteIter.Set(sliceIter.Get());
                ++pasteIter;
                ++sliceIter;
            }
        }
        __labelIO.WriteImage(output, pasteImage);
        return;
    }


    int CommandLineTools::Run(Options* parser, StringVector args) {
        if (parser->GetBool("--isoRG")) {
            int fg = parser->GetStringAsInt("-f", 2);
            int bg = parser->GetStringAsInt("-b", 3);
            if (args.size() < 3) {
                cout << "Isolated Connected Filter: -f [fgId=2] -b [bgId=3] [input-image] [input-label] [output-label]" << endl;
                return 0;
            }

            ImageAlgo algo;
            Image::Pointer gray = __imageIO.ReadCastedImage(args[0]);
            Label::Pointer label = __labelIO.ReadCastedImage(args[1]);
            Label::Pointer labelOut = algo.ExecuteIsolatedConnectedImageFilter(label, gray, fg, bg);

            __labelIO.WriteImage(args[2], labelOut);
        } else if (parser->GetBool("--extractSlice")) {
            if (args.size() < 2) {
                cout << "ExtractSlice --dir IJ|JK|KI --range 10,100 [input-image] [output-image-pattern: slice-%1.nii.gz]" << endl;
                return 0;
            }

            Image::Pointer inputImg = __imageIO.ReadCastedImage(args[0]);
            if (inputImg.IsNull()) {
                cout << "can't read " << args[0] << endl;
                return 0;
            }
            string dirString = parser->GetString("--dir");

            SliceDirectionEnum dir = Unknown;
            if (dirString == "IJ") {
                dir = IJ;
            } else if (dirString == "JK") {
                dir = JK;
            } else if (dirString == "KI") {
                dir = KI;
            }

            IntVector rangeValues = parser->GetStringAsIntVector("--range");
            ExtractSlice(inputImg, dir, rangeValues, args[1]);
        } else if (parser->GetBool("--pasteSlice")) {
            if (args.size() < 2) {
                cout << "PasteSlice [output-volume] [input1] [input2] ... [inputN]" << endl;
                return 0;
            }
            PasteSlice(args, args[0]);
        } else if (parser->GetBool("--pasteLabel")) {
            if (args.size() < 2) {
                cout << "PasteLabel [output-volume] [input1] [input2] ... [inputN]" << endl;
                return 0;
            }
            PasteLabel(args, args[0]);
        }
        return 0;
    }
}