//
//  itkMySlice.h
//  itktools
//
//  Created by Joohwi Lee on 9/25/12.
//
//

#ifndef __itktools__itkMySlice__
#define __itktools__itkMySlice__

#include <iostream>

#include "itkObject.h"
#include "itkExtractImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkImageIO.h"


template <class TImage, class TSlice>
class itkMySlice : public itk::Object {
private:
    double m_MinIntensity;
    double m_MaxIntensity;
    double m_MeanIntensity;
    double m_StdIntensity;

public:
    typedef itkMySlice Self;
    typedef itk::Object Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;
    typedef itk::StatisticsImageFilter<TImage> StatisticsImageFilterType;
    typedef typename TSlice::Pointer SlicePointer;

private:
    typename TImage::Pointer m_Volume;
    typename TSlice::Pointer m_Slice[TImage::ImageDimension];

    typename TImage::IndexType _currentSliceIndex;
    typename TImage::SizeType _maxSliceIndex;
    typename StatisticsImageFilterType::Pointer _statisticsImageFilter;

public:
    itkNewMacro(itkMySlice);
    itkTypeMacro(itkMySlice, itk::Object);

    itkGetConstMacro(MinIntensity, double);
    itkGetConstMacro(MaxIntensity, double);
    itkGetConstMacro(MeanIntensity, double);
    itkGetConstMacro(StdIntensity, double);

    typename TSlice::SizeType GetSize(int dir) {
        return m_Slice[dir]->GetBufferedRegion().GetSize();
    }
    
	typename TImage::IndexType GetCurrentSliceIndex() {
		return _currentSliceIndex;
	}

    typename TImage::SizeType GetMaxSliceIndex() {
        return _maxSliceIndex;
    }

    SlicePointer GetSlicePointer(int dir) {
        return m_Slice[dir];
    }

    typename TImage::Pointer GetVolume() {
        return m_Volume;
    }

    void SetVolume(typename TImage::Pointer volume) {
        m_Volume = volume;
        _maxSliceIndex = m_Volume->GetBufferedRegion().GetSize();

        _statisticsImageFilter = StatisticsImageFilterType::New ();
        _statisticsImageFilter->SetInput(volume);
        _statisticsImageFilter->Update();

        m_MeanIntensity = _statisticsImageFilter->GetMean();
        m_StdIntensity = _statisticsImageFilter->GetSigma();
        m_MinIntensity = _statisticsImageFilter->GetMinimum();
        m_MaxIntensity = _statisticsImageFilter->GetMaximum();

        for (int i = 0; i < TImage::ImageDimension; i++) {
            _currentSliceIndex[i] = -1;
            Update(i, _maxSliceIndex[i] / 2);
        }
    }

    /**
     * Extract a slice of direction 'dir' at the sliceIdx.
     * The update will only occurs when the current slice index is not equal to the new slice idx.
     * If the slice is null, then the update must occur.
     *
     */
    void Update(int dir, int sliceIdx) {
        if (_currentSliceIndex[dir] == sliceIdx && m_Slice[dir].IsNotNull()) {
            return;
        }
         _currentSliceIndex[dir] = sliceIdx;

        //cout << "Current Slice Idx[" << dir << "] : " << sliceIdx << endl;
        typedef itk::ExtractImageFilter<TImage, TSlice> SliceExtractor;
        typename SliceExtractor::Pointer slicer = SliceExtractor::New();
        typename TImage::RegionType sliceRegion;
        typename TImage::SizeType sliceSize;
        typename TImage::IndexType sliceIndex;
        
        for (int i = 0; i < TImage::ImageDimension; i++) {
            if (i == dir) {
                sliceSize[i] = 0;
                sliceIndex[i] = sliceIdx;
            } else {
                sliceSize[i] = _maxSliceIndex[i];
                sliceIndex[i] = 0;
            }
        }

        sliceRegion.SetSize(sliceSize);
        sliceRegion.SetIndex(sliceIndex);
        
        slicer->SetInput(m_Volume);
        slicer->SetExtractionRegion(sliceRegion);
        slicer->SetDirectionCollapseToGuess();
        slicer->Update();
        m_Slice[dir] = slicer->GetOutput();

        this->Modified();
    }

    void SaveCurrentSlice(const char* filename, int dir) {
        itkcmds::itkImageIO<TSlice> itkIO;
        itkIO.WriteImageT(filename, m_Slice[dir]);
    }

protected:
    itkMySlice() {
        for (int i = 0; i < TImage::ImageDimension; i++) {
            _currentSliceIndex[i] = -1;
            _maxSliceIndex[i] = -1;
        }
    }

    virtual ~itkMySlice() {

    }

};

#endif /* defined(__itktools__itkMySlice__) */
