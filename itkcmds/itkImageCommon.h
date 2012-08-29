#ifndef __ITK_IMAGE_COMMON_H__
#define __ITK_IMAGE_COMMON_H__

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericTraits.h"
#include "itkNrrdImageIOFactory.h"
#include "itkNiftiImageIOFactory.h"
#include "iostream"


void initReaders() {
	itk::ObjectFactoryBase::RegisterFactory(itk::NrrdImageIOFactory::New());
	itk::ObjectFactoryBase::RegisterFactory(itk::NiftiImageIOFactory::New());
}

template<class T> typename T::Pointer NewImageT(int sx, int sy, int sz, typename T::PixelType fillValue) {
	typename T::Pointer newImage = T::New();
	typename T::SizeType size;
	size[0] = sx;
	size[1] = sy;
	size[2] = sz;

	typename T::RegionType region;
	region.SetSize(size);
	
	typename T::IndexType index;
	index.Fill(0);

	region.SetIndex(index);
	
	newImage->SetLargestPossibleRegion(region);
	newImage->SetBufferedRegion(region);
	newImage->SetRequestedRegion(region);

	
	double spacing[3] = {1, 1, 1};
	double origin[3] = { 0, 0, 0 };
	newImage->SetOrigin(origin);
	newImage->SetSpacing(spacing);

	newImage->Allocate();
	newImage->FillBuffer(fillValue);

	return newImage;
}

template<class T> typename T::Pointer NewImageT(int sx, int sy, int sz) {
  return NewImageT<T>(sx, sy, sz, static_cast<typename T::PixelType>(0));
}

template<class T> void CopyHeaderT(typename T::Pointer src, typename T::Pointer dst) {
	dst->SetSpacing(src->GetSpacing());
	dst->SetOrigin(src->GetOrigin());
	dst->SetDirection(src->GetDirection());
	return;
}

template<class T> typename T::Pointer NewImageT(typename T::Pointer srcImg) {
	typename T::RegionType srcRegion = srcImg->GetRequestedRegion();
	typename T::SizeType srcSize = srcRegion.GetSize();
  typename T::Pointer newImg =  NewImageT<T>(srcSize[0], srcSize[1], srcSize[2], static_cast<typename T::PixelType>(0));
	CopyHeaderT<T>(srcImg, newImg);
	return newImg;
}

template<class T> bool DumpImageT(typename T::Pointer src, typename T::Pointer dst) {
  typename T::RegionType srcRegion = src->GetLargestPossibleRegion();
  typename T::RegionType dstRegion = dst->GetLargestPossibleRegion();

  if (srcRegion != dstRegion) {
    return false;
  }

  itk::ImageRegionConstIterator<T> srcIter(src, srcRegion);
  itk::ImageRegionIterator<T> dstIter(dst, dstRegion);

  for (; !srcIter.IsAtEnd() && !dstIter.IsAtEnd(); ++srcIter, ++dstIter) {
    dstIter.Set(srcIter.Get());
  }
  return true;
}

template<class T> typename T::Pointer ReadImageT(const char* filename, int& ret) {
	initReaders();
	typename itk::ImageFileReader<T>::Pointer reader = itk::ImageFileReader<T>::New();
	reader->SetFileName(filename);

	std::cout << "Reading " << filename;
	std::cout.flush();
	reader->Update();
	std::cout << " done." << std::endl;
	return reader->GetOutput();
}

template<class T> int WriteImageT(const char* filename, typename T::Pointer image, bool compression) {
	typename itk::ImageFileWriter<T>::Pointer writer = itk::ImageFileWriter<T>::New();
	writer->SetFileName(filename);
  if (compression) {
    writer->UseCompressionOn();
  }
	writer->SetInput(image);
  std::cout << "Writing " << filename;
  std::cout.flush();
	writer->Write();
  std::cout << " done." << std::endl;
	return 0;
}

template<class T> int WriteImageT(const char* filename, typename T::Pointer image) {
  WriteImageT<T>(filename, image, true);
	return 0;
}
#endif
