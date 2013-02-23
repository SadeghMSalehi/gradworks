#ifndef __ITK_IMAGE_COMMON_H__
#define __ITK_IMAGE_COMMON_H__

#include "itkImage.h"
#include "itkImageIOBase.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericTraits.h"
#include "itkNrrdImageIOFactory.h"
#include "itkNiftiImageIOFactory.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkExceptionObject.h"
#include "itkImportImageFilter.h"
#include "itkMath.h"
#include "iostream"

using namespace std;

namespace itkcmds {
    static int __noverbose = 0;

	template <typename T>
	class itkImageIO {
	private:
		itk::ImageIOBase::IOPixelType _pixelType;
		itk::ImageIOBase::IOComponentType _componentType;
        
	public:
		typedef typename T::Pointer ImagePointer;
		typedef typename T::PixelType ImagePixel;
		typedef typename T::SizeType ImageSize;
		typedef typename T::RegionType ImageRegion;
		typedef typename T::IndexType ImageIndex;
        typedef itk::TransformFileReader TransformFileReader;
        typedef itk::TransformFileWriter TransformFileWriter;
        typedef TransformFileReader::TransformType TransformType;
        typedef TransformFileReader::TransformListType   TransformListType;
        
        
		itkImageIO() {
			_pixelType = itk::ImageIOBase::SCALAR;
			_componentType = itk::ImageIOBase::DOUBLE;
            
			itk::ObjectFactoryBase::RegisterFactory(itk::NrrdImageIOFactory::New());
			itk::ObjectFactoryBase::RegisterFactory(itk::NiftiImageIOFactory::New());
		}
        
		~itkImageIO() {
		}
        
		typename itk::ImageIOBase::Pointer ReadImageInfo(const char* filename) {
			typename itk::ImageFileReader<T>::Pointer reader = itk::ImageFileReader<T>::New();
			reader->SetFileName(filename);
			return GetImageInfo(reader);
		}
        
		typename itk::ImageIOBase::Pointer GetImageInfo(typename itk::ImageFileReader<T>::Pointer reader) {
			reader->UpdateOutputInformation();
			typename itk::ImageIOBase::Pointer imageIO = reader->GetImageIO();
			_pixelType = imageIO->GetPixelType();
			_componentType = imageIO->GetComponentType();
            return imageIO;
		}
        
		const char* GetPixelTypeString(typename itk::ImageIOBase::IOPixelType px) {
			static const char* pixelTypes[] = { "Unknown Pixel Type", "Scalar", "RGB", "RGBA",
				"Offset", "Vector", "Point", "Covariant Vector", "Symmetric Second Rank Tensor",
				"Diffusion Tensor 3D", "Complex", "Fixed Array", "Matrix" };
			return pixelTypes[px];
		}
        
		const char* GetComponentTypeString(typename itk::ImageIOBase::IOComponentType cx) {
			static const char* componentTypes[] = { "Unknown Component Type", "UCHAR", "CHAR",
				"USHORT", "SHORT", "UINT", "INT", "ULONG", "LONG", "FLOAT", "DOUBLE" };
			return componentTypes[cx];
		}
        
		ImagePointer NewImageT(int sx, int sy, int sz, ImagePixel fillValue) {
			ImagePointer newImage = T::New();
			ImageSize size;
            
			size[0] = sx;
			size[1] = sy;
            if (T::GetImageDimension() == 3) {
                size[2] = sz;
            }
            
			ImageRegion region;
			region.SetSize(size);
            
			ImageIndex index;
			index.Fill(0);
            
			region.SetIndex(index);
            
			newImage->SetLargestPossibleRegion(region);
			newImage->SetBufferedRegion(region);
			newImage->SetRequestedRegion(region);
            
            if (T::GetImageDimension() == 3) {
                double spacing[3] = { 1, 1, 1 };
                double origin[3] = { 0, 0, 0 };
                newImage->SetOrigin(origin);
                newImage->SetSpacing(spacing);
            } else if (T::GetImageDimension() == 2) {
                double spacing[2] = { 1, 1 };
                double origin[2] = { 0, 0 };
                newImage->SetOrigin(origin);
                newImage->SetSpacing(spacing);
            }
            
			newImage->Allocate();
			newImage->FillBuffer(fillValue);
            
			return newImage;
		}
        
		ImagePointer NewImageT(int sx, int sy, int sz) {
			return NewImageT(sx, sy, sz, static_cast<typename T::PixelType>(0));
		}
        
        template <typename S>
		void CopyHeaderT(typename S::Pointer src, ImagePointer dst) {
			dst->SetSpacing(src->GetSpacing());
			dst->SetOrigin(src->GetOrigin());
			dst->SetDirection(src->GetDirection());
			return;
		}
        
		ImagePointer NewImageT(ImagePointer srcImg) {
            return NewImageT<T>(srcImg);
        }
        
        template <class S>
        ImagePointer NewImageT(typename S::Pointer srcImg) {
			typename S::RegionType srcRegion = srcImg->GetRequestedRegion();
			typename S::SizeType srcSize = srcRegion.GetSize();
			typename T::Pointer newImg = NewImageT(srcSize[0], srcSize[1], srcSize[2], static_cast<typename T::PixelType>(0));
			CopyHeaderT<S>(srcImg, newImg);
			return newImg;
		}
        
		bool DumpImageT(ImagePointer src, ImagePointer dst) {
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
        
		ImagePointer ReadImageT(const char* filename) {
            if (!__noverbose) {
                cout << "Reading '" << filename << flush;
            }
            if (CheckExists(filename)) {
                typename itk::ImageFileReader<T>::Pointer reader = itk::ImageFileReader<T>::New();
                reader->SetFileName(filename);
                reader->Update();
                GetImageInfo(reader);
                if (!__noverbose) {
                    std::cout << "' [" << GetComponentTypeString(_componentType) << ", " << GetPixelTypeString(_pixelType) << "]";
                    std::cout << " done." << std::endl;
                }
                return reader->GetOutput();
            } else {
                if (!__noverbose) {
                    cout << "' failed. (file not exist)" << endl;
                }
                return ImagePointer();
            }
		}
        
		int WriteImageT(const char* filename, ImagePointer image, bool compression) {
			typename itk::ImageFileWriter<T>::Pointer writer = itk::ImageFileWriter<T>::New();
			writer->SetFileName(filename);
			if (compression) {
				writer->UseCompressionOn();
			}
			writer->SetInput(image);
			std::cout << "Writing '" << filename;
			std::cout.flush();
            try {
                writer->Write();
            } catch (itk::ExceptionObject& e) {
                e.Print(cout);
            }
			std::cout << "' done." << std::endl;
			return 0;
		}
        
		int WriteImageT(const char* filename, typename T::Pointer image) {
			WriteImageT(filename, image, true);
			return 0;
		}
        
        typename TransformType::Pointer ReadTransform(char* fileName) {
            typename TransformFileReader::Pointer transformFileReader = TransformFileReader::New();
            transformFileReader->SetFileName(fileName);
            // Create the transforms
            transformFileReader->Update();
            TransformListType* transformList = transformFileReader->GetTransformList();
            return transformList->front();
        }
        
        void WriteSingleTransform(char* fileName, typename TransformType::Pointer transform) {
            typename TransformFileWriter::Pointer writer = TransformFileWriter::New();
            writer->SetFileName(fileName);
            writer->AddTransform(transform);
            writer->Update();
        }
        
        void WriteMultipleTransform(char* fileName, typename TransformFileWriter::ConstTransformListType transformList) {
            typename TransformFileWriter::Pointer writer = TransformFileWriter::New();
            writer->SetFileName(fileName);
            typename TransformFileWriter::ConstTransformListType::iterator transformIter = transformList.begin();
            while (transformIter != transformList.end()) {
                writer->AddTransform(*transformIter);
                transformIter++;
            }
            writer->Update();
        }
        
        typename T::Pointer ResampleImageAs(typename T::Pointer image, typename T::Pointer reference) {
            typedef itk::ResampleImageFilter<T,T> FilterType;
            typedef itk::LinearInterpolateImageFunction<T> InterpolatorType;
            typename FilterType::Pointer filter = FilterType::New();
            filter->SetInput(image);
            filter->SetReferenceImage(reference);
            filter->UseReferenceImageOn();
            filter->SetInterpolator(InterpolatorType::New());
            filter->Update();
            return filter->GetOutput();
        }
        
       	bool FileExists(const char* fileName) {
			ifstream ifile(fileName);
			return ifile != NULL;
		}
        
        bool CheckExists(const char* filename) {
            ifstream fin(filename);
            return fin != NULL;
        }
        
        static int GetImageDimension() {
            return T::GetImageDimension();
        }
        
        static typename T::Pointer ImportFromMemory3(typename T::PixelType* src, int num, bool freememory, typename T::SizeType sz) {
            typename T::RegionType region;
            region.SetSize(sz);
            typedef itk::ImportImageFilter<typename T::PixelType, 3> FilterType;
            typename FilterType::Pointer filter = FilterType::New();
            filter->SetRegion(region);
            filter->SetImportPointer(src, num, freememory);
            filter->Update();
            return filter->GetOutput();
        }

        static typename T::Pointer ImportFromMemory2(typename T::PixelType* src, int num, bool freememory, typename T::SizeType sz) {
            typename T::RegionType region;
            region.SetSize(sz);
            typedef itk::ImportImageFilter<typename T::PixelType, 2> FilterType;
            typename FilterType::Pointer filter = FilterType::New();
            filter->SetRegion(region);
            filter->SetImportPointer(src, num, freememory);
            filter->Update();
            return filter->GetOutput();
        }
        
    };
    
}
#endif



