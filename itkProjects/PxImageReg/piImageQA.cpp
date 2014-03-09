//
//  piImageQA.cpp
//  PxImageReg
//
//  Created by Joohwi Lee on 11/9/13.
//
//

#include "piImageQA.h"
#include "piOptions.h"
#include "piImageHistogram.h"
#include "libconfig.h++"

#include "itkARGBColorFunction.h"
#include "itkScalarToARGBColormapImageFilter.h"

#include <itkExtractImageFilter.h>
#include <itkFlipImageFilter.h>

namespace pi {
    typedef itk::RGBAPixel<unsigned char> RGBPixel;
    typedef itk::Image<RGBPixel, 2> RGBImage2;
    typedef itk::Image<RGBPixel, 3> RGBImage3;

    typedef itk::ScalarToARGBColormapImageFilter<LabelImage3, RGBImage3> LabelToRGBFilter;


    /// execute QA image processing
    void executeQARunner(Options& parser, StringVector& args) {
        if (args.size() < 5) {
            cout << "--qa input-image input-label slice axis label-alpha output-image [--config config-file]" << endl;
            return;
        }

        ImageQA qa;
        qa.loadConfig(parser.GetString("--config"));
        qa.slice = atoi(args[2].c_str());
        qa.axis = atoi(args[3].c_str());
        qa.alpha = atof(args[4].c_str()) / 255.0;
        qa.load(args[0], args[1]);
        qa.save(args[5]);
    }

    /// ABGR
    static uint __labelColors[254] = { 0,
        0xff0000ff, 0xff00ff00, 0xffff0000, 0xff00ffff, 0xffffff00, 0xffff00ff, 0xffd5efff, 0xffcd0000, 0xff3f85cd,
        0xff8cb4d2, 0xffaacd66, 0xff8b0000, 0xff8b8b00, 0xff578b2e, 0xffe1e4ff, 0xffcd5a6a, 0xffdda0dd, 0xff7a96e9,
        0xff2a2aa5, 0xfffafaff, 0xffdb7093, 0xffd670da, 0xff4b0082, 0xff3a9708, 0xffc6853c, 0xfff84707, 0xff73dda6,
        0xff8392b4, 0xff54f584
    };


    /// templatead slice extraction function
    /// @dir 0-(yz), 1-(zx), 2-(xy)
    template <class T3, class T2>
    static typename T2::Pointer ExtractSlice(typename T3::Pointer srcImg, int idx, int dir) {
        typename T2::Pointer emptyImg;
        if (srcImg.IsNull() || dir == Unknown) {
            cout << "source is null or direction unknown:" << __FILE__ << ":" << __LINE__ << endl;
            return emptyImg;
        }

        typename T3::RegionType sliceRegion = srcImg->GetBufferedRegion();
        typename T3::SizeType srcSize = sliceRegion.GetSize();
        if (srcSize[dir] <= idx || idx < 0) {
            cout << "slice index is wrong: " << idx << " >= " << srcSize[dir] << __FILE__ << ":" << __LINE__ << endl;
            return emptyImg;
        }
        sliceRegion.SetIndex(dir, idx);
        sliceRegion.SetSize(dir,0);

        typedef itk::ExtractImageFilter<T3, T2> ExtractFilterType;
        typename ExtractFilterType::Pointer filter = ExtractFilterType::New();
        filter->SetInput(srcImg);
        filter->SetExtractionRegion(sliceRegion);
        filter->SetDirectionCollapseToGuess();
        filter->Update();
        typename T2::Pointer sliceImg = filter->GetOutput();
        sliceImg->DisconnectPipeline();


        typedef itk::FlipImageFilter<T2> FlipFilter;
        typename FlipFilter::FlipAxesArrayType flipAxes;
        flipAxes.Fill(false);
        if (dir == 1) {
            flipAxes[1] = true;
        } else if (dir == 0) {
            flipAxes[1] = true;
        }
        typename FlipFilter::Pointer flipper = FlipFilter::New();
        flipper->SetInput(sliceImg);
        flipper->SetFlipAxes(flipAxes);
        flipper->Update();
        sliceImg = flipper->GetOutput();
        sliceImg->DisconnectPipeline();

        return sliceImg;
    }

    RGBImage2::Pointer convertImageToRGB(RealImage2::Pointer image, RealImage2::PixelType rangeMin, RealImage2::PixelType rangeMax) {
        typedef itk::ScalarToARGBColormapImageFilter<RealImage2, RGBImage2> ImageToRGBFilter;

        ImageToRGBFilter::Pointer rgbFilter = ImageToRGBFilter::New();
        rgbFilter->SetInput(image);
        rgbFilter->SetColormap(itk::Grey);


        rgbFilter->SetMinimumValue(rangeMin);
        rgbFilter->SetMaximumValue(rangeMax);
        rgbFilter->SetUseManualScaling(true);

        rgbFilter->Update();
        RGBImage2::Pointer output = rgbFilter->GetOutput();
        output->DisconnectPipeline();
        return output;
    }

    RGBImage2::Pointer convertLabelToRGB(LabelImage2::Pointer label) {
        ImageIO<RGBImage2> io;
        RGBImage2::Pointer newImage = io.NewImageS<LabelImage2>(label);

        LabelImage2::PixelType* lPtr = label->GetBufferPointer();
        RGBImage2::PixelType* rPtr = newImage->GetBufferPointer();

        for (int i = 0; i < label->GetPixelContainer()->Size(); i++) {
            unsigned int *rgbptr = reinterpret_cast<unsigned int*>(rPtr[i].GetDataPointer());
            *rgbptr = __labelColors[lPtr[i]];
        }

        return newImage;
    }

    /// Composite two rgb images with given alpha value
    /// @param input1
    /// @param input2
    /// @param alpha float applied to input2
    RGBImage2::Pointer compositeImages2(RGBImage2::Pointer input1, RGBImage2::Pointer input2, float alpha) {
        ImageIO<RGBImage2> io;
        RGBImage2::Pointer outputImage = io.NewImageS<RGBImage2>(input1);

        RGBImage2::PixelType* i1 = input1->GetBufferPointer();
        RGBImage2::PixelType* i2 = input2->GetBufferPointer();
        RGBImage2::PixelType* o = outputImage->GetBufferPointer();

        for (int i = 0; i < input1->GetPixelContainer()->Size(); i++) {
            // Alpha
            o[i][3] = 255;
            // rgb
            for (int k = 0; k < 3; k++) {
                float labelAlpha = alpha * float(i2[i][3])/255.0;
                o[i][k] = ::round((1 - alpha) * i1[i][k] + labelAlpha * i2[i][k]);
            }
        }

        return outputImage;
    }


    void ImageQA::loadConfig(std::string file) {
        using namespace libconfig;
        if (file == "") {
            return;
        }
        Config config;
        config.readFile(file.c_str());
    }

    void ImageQA::load(std::string image, std::string label) {
        ImageIO<RealImage3> io;
        this->image = io.ReadCastedImage(image);
        this->label = io.ReadImageS<LabelImage3>(label);
    }

    void ImageQA::save(std::string output) {
        /*
        RGBImage3::Pointer colorImage = convertImageToRGB(this->image);
        RGBImage3::Pointer colorLabel = convertLabelToRGB(this->label);

        RGBImage2::Pointer sliceImage = ExtractSlice<RGBImage3, RGBImage2>(colorImage, slice, axis);
        RGBImage2::Pointer sliceLabel = ExtractSlice<RGBImage3, RGBImage2>(colorLabel, slice, axis);
         */


        // first extract a given slice from the volume
        RealImage2::Pointer sliceImage = ExtractSlice<RealImage3, RealImage2>(this->image, slice, axis);

        LabelImage2::Pointer sliceLabel = ExtractSlice<LabelImage3, LabelImage2>(this->label, slice, axis);


        // compute 10% intensity range
        ImageHistogram<RealImage3> imageHistogram;
        imageHistogram.SetImage(this->image);
        imageHistogram.FitRange();

        RGBImage2::Pointer colorImage = convertImageToRGB(sliceImage, imageHistogram.rangeMin, imageHistogram.rangeMax);

        RGBImage2::Pointer colorLabel = convertLabelToRGB(sliceLabel);


        RGBImage2::Pointer sliceOverlay = compositeImages2(colorImage, colorLabel, alpha);

        ImageIO<RGBImage2> io;
        io.WriteImage(output, sliceOverlay);
    }
}
