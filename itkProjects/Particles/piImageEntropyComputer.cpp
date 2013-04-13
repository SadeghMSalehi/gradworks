//
//  piImageEntropy.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/12/13.
//
//

#include "piImageEntropyComputer.h"
#include <itkThreadedImageRegionPartitioner.h>
#include <itkDomainThreader.h>
#include <itkExtractImageFilter.h>
#include "piImageIO.h"

using namespace std;

static pi::ImageIO<pi::RealImage> io;
typedef itk::ExtractImageFilter<pi::RealImage, pi::RealImage> FilterType;


static pi::RealImage::Pointer extractImage(pi::RealImage::Pointer image, const pi::RealImage::RegionType region, pi::RealImage::Pointer output) {
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(image);
    filter->SetExtractionRegion(region);
    try {
        filter->Update();
    } catch (itk::ExceptionObject& ex) {
        ex.Print(cout);
        cout << image->GetBufferedRegion() << endl;
        cout << region << endl;
    }
    output = filter->GetOutput();
    return output;
}

namespace itk {
    typedef std::vector<pi::RealImage::RegionType> RegionStack;

    class ImageNeighborhoodRegionPartitioner: public ThreadedDomainPartitioner<RegionStack> {
    public:
        typedef pi::RealImage::RegionType RegionType;
        typedef pi::RealImage::SizeType SizeType;
        typedef pi::RealImage::IndexType IndexType;

        enum { ImageDimension = pi::RealImage::ImageDimension };

        typedef ImageNeighborhoodRegionPartitioner Self;
        typedef ThreadedDomainPartitioner Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro(Self);
        itkTypeMacro(ImageNeighborhoodRegionPartitioner, ThreadedDomainPartitioner<RegionStack>);

        virtual ThreadIdType PartitionDomain(const ThreadIdType i, const ThreadIdType requestedTotal, const DomainType& completeRegion, DomainType& subRegion) const {
            const uint32_t regionsCount = completeRegion.size();
            const uint32_t regionsPerThread = regionsCount / (requestedTotal);

            if (completeRegion.size() < requestedTotal) {
                subRegion = completeRegion;
                return 1;
            }
            int nBegin = regionsPerThread * i;
            int nEnd = (i + 1 == requestedTotal) ? regionsCount : regionsPerThread * (i + 1);
            for (int j = nBegin; j < nEnd; j++) {
                subRegion.push_back(completeRegion[j]);
            }
            return requestedTotal;
        }

        void CreateCompleteRegion(RegionType maskRegion, RegionType targetRegion, DomainType& completeRegion) {
            SizeType regionsSize = targetRegion.GetSize() - maskRegion.GetSize();
            int nElems = 1;
            for (int j = 0; j < ImageDimension; j++) {
                nElems *= regionsSize[j];
            }
            IndexType loopIndex = targetRegion.GetIndex();
            for (int j = 0; j < nElems; j++) {
                completeRegion.push_back(RegionType(loopIndex, maskRegion.GetSize()));
                for (int k = 0; k < ImageDimension; k++) {
                    loopIndex[k] ++;
                    if (loopIndex[k] + maskRegion.GetSize(k) < targetRegion.GetIndex(k) + targetRegion.GetSize(k)) {
                        break;
                    }
                    loopIndex[k] = targetRegion.GetIndex(k);
                }
            }
        }

    protected:
        ImageNeighborhoodRegionPartitioner() {}
        virtual ~ImageNeighborhoodRegionPartitioner() {}

    private:
        ImageNeighborhoodRegionPartitioner(const ImageNeighborhoodRegionPartitioner&);
        void operator=(const ImageNeighborhoodRegionPartitioner&);
    };

    class ThreadedImageEntropyComputer: public DomainThreader<ImageNeighborhoodRegionPartitioner, pi::ImageHolder> {
    public:
        typedef ThreadedImageEntropyComputer Self;
        typedef SmartPointer<Self> Pointer;
        itkNewMacro(Self);
        itkTypeMacro(ThreadedImageEntropyComputer, DomainThreader);

        void ThreadedExecution(const DomainType& sub, const ThreadIdType threadId) {
            pi::RealImage::Pointer imagePatch;
            pi::ImageEntropyComputer comp;

            for (int i = 0; i < sub.size(); i++) {
                imagePatch = extractImage(m_Associate->sourceImage, sub[i], imagePatch);

                comp.clear();
                comp.setSize(2, imagePatch->GetPixelContainer()->Size());
                comp.addSample(m_Associate->referencePatch->GetBufferPointer());
                comp.addSample(imagePatch->GetBufferPointer());
                comp.computeEntropy();

                pi::DataReal value = comp.entropyValue();
                pi::RealImage::IndexType outputIndex = sub[i].GetIndex();
                m_Associate->output->SetPixel(outputIndex, value);
            }
        }

    protected:
        ThreadedImageEntropyComputer() {}
        virtual ~ThreadedImageEntropyComputer() {}

    private:
        ThreadedImageEntropyComputer(const ThreadedImageEntropyComputer&);
        void operator=(const ThreadedImageEntropyComputer&);
    };
}

namespace pi {
    ImageEntropyComputer::ImageEntropyComputer() {
        _value = 0;
        _m = _n = 0;
    }

    void ImageEntropyComputer::setSize(int mData, int nSamples) {
        _m = mData;
        _n = nSamples;

        _meanStore.set_size(_n);
        _dataStore.set_size(_m, _n);
        _covStore.set_size(_m, _m);
    }

    void ImageEntropyComputer::addSample(DataReal* sampleBuffer) {
        _data.push_back(sampleBuffer);
    }

    void ImageEntropyComputer::clear() {
        _data.clear();
        _value = 0;
    }

    double ImageEntropyComputer::entropyValue() {
        return _value;
    }

    void ImageEntropyComputer::computeEntropy() {
        const int n = _n;
        const int m = _m;

        _meanStore.fill(0);
        for (int j = 0; j < m; j++) {
            DataReal* mean = _meanStore.data_block();
            DataReal* pixels = _data[j];
            for (int i = 0; i < n; i++, mean++) {
                *mean += *pixels;
            }
        }
        _meanStore /= m;


        DataReal* data = _dataStore.data_block();
        for (int j = 0; j < m; j++) {
            DataReal* mean = _meanStore.data_block();
            DataReal* pixels = _data[j];
            for (int i = 0; i < n; i++, data++, mean++) {
                *data = *pixels - *mean;
            }
        }

        DataReal* cov = _covStore.data_block();
        // i-th row
        for (int i = 0; i < _m; i++) {
            // k-th column
            for (int k = 0; k < _m; k++) {
                // j-th element
                for (int j = 0; j < _n; j++) {
                    *cov = _dataStore[i][j] * _dataStore[k][j];
                }
                ++cov;
            }
        }

        vnl_symmetric_eigensystem_compute(_covStore, _V, _D);
        _value = _D.sum();
    }

    RealImage::Pointer ImageEntropyComputer::computeEntropy(RealImage::Pointer source, RealImage::RegionType sourceMask, RealImage::Pointer target, RealImage::RegionType targetRegion) {
        itk::ImageNeighborhoodRegionPartitioner::Pointer part = itk::ImageNeighborhoodRegionPartitioner::New();
        itk::RegionStack regions;

        ImageHolder holder;
        holder.sourceImage = target;
        holder.output = io.CopyImage(target);
        holder.output->FillBuffer(0);
        holder.referencePatch = extractImage(source, sourceMask, holder.referencePatch);
        part->CreateCompleteRegion(sourceMask, targetRegion, regions);

        itk::ThreadedImageEntropyComputer::Pointer comp = itk::ThreadedImageEntropyComputer::New();
//        comp->SetMaximumNumberOfThreads(2);
        comp->Execute(&holder, regions);
        
        return holder.output;
    }

}