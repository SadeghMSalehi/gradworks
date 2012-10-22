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
#include "itkMath.h"
#include "iostream"

using namespace std;

namespace itkcmds {

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

        bool FileExists(const char* fileName) {
            ifstream ifile(fileName);
            return ifile;
        }

		ImagePointer NewImageT(int sx, int sy, int sz, ImagePixel fillValue) {
			ImagePointer newImage = T::New();
			ImageSize size;
			size[0] = sx;
			size[1] = sy;
			size[2] = sz;

			ImageRegion region;
			region.SetSize(size);

			ImageIndex index;
			index.Fill(0);

			region.SetIndex(index);

			newImage->SetLargestPossibleRegion(region);
			newImage->SetBufferedRegion(region);
			newImage->SetRequestedRegion(region);

			double spacing[3] = { 1, 1, 1 };
			double origin[3] = { 0, 0, 0 };
			newImage->SetOrigin(origin);
			newImage->SetSpacing(spacing);

			newImage->Allocate();
			newImage->FillBuffer(fillValue);

			return newImage;
		}

		ImagePointer NewImageT(int sx, int sy, int sz) {
			return NewImageT(sx, sy, sz, static_cast<typename T::PixelType>(0));
		}

		void CopyHeaderT(ImagePointer src, ImagePointer dst) {
			dst->SetSpacing(src->GetSpacing());
			dst->SetOrigin(src->GetOrigin());
			dst->SetDirection(src->GetDirection());
			return;
		}

		ImagePointer NewImageT(ImagePointer srcImg) {
			typename T::RegionType srcRegion = srcImg->GetRequestedRegion();
			typename T::SizeType srcSize = srcRegion.GetSize();
			typename T::Pointer newImg = NewImageT(srcSize[0], srcSize[1], srcSize[2], static_cast<typename T::PixelType>(0));
			CopyHeaderT(srcImg, newImg);
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
			cout << "Reading " << filename << flush;
			typename itk::ImageFileReader<T>::Pointer reader = itk::ImageFileReader<T>::New();
			reader->SetFileName(filename);
			reader->Update();
			GetImageInfo(reader);
            std::cout << " [" << GetComponentTypeString(_componentType) << ", " << GetPixelTypeString(_pixelType) << "]";
			std::cout << " done." << std::endl;
			return reader->GetOutput();
		}

		int WriteImageT(const char* filename, ImagePointer image, bool compression) {
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
};

template<class T, unsigned int N>
class PointStatistics {
private:
    T _minPoint;
    T _maxPoint;
    T _avgPoint;
    T _sumPoint;
    int n;
public:
    PointStatistics() {
        for (int i = 0; i < N; i++) {
            _sumPoint[i] = 0;
            _minPoint[i] = 1e23;
            _maxPoint[i] = -1e23;
            _avgPoint[i] = 0;
        }
        n = 0;
    }
    ~PointStatistics() {}
    int GetNumberOfPoints() const {
        return n;
    }
    T GetBoundsMinimum() const {
        return _minPoint;
    }

    T GetBoundsMaximum() const {
        return _maxPoint;
    }

    T GetBoundsCenter() const {
        return _avgPoint;
    }

    T GetSumOfPoints() {
        return _sumPoint;
    }

    void AddPoint(T p) {
        ++n;
        for (int i = 0; i < N; i++) {
            _minPoint[i] = _minPoint[i] < p[i] ? _minPoint[i] : p[i];
            _maxPoint[i] = _maxPoint[i] > p[i] ? _maxPoint[i] : p[i];
            _sumPoint[i] += p[i];
            _avgPoint[i] = (_sumPoint[i] / n);
        }
    }
};

template<class T, unsigned int n>
std::ostream& operator<<(std::ostream& out, const PointStatistics<T,n>& s) {
    out << "Center: " << s.GetBoundsCenter() << "; Min: " << s.GetBoundsMinimum() << "; Max: " << s.GetBoundsMaximum() << ";";
}

template<class T, class S>
void ComputeLabelIndexes(typename T::Pointer label, typename S::FixedImageIndexContainer& container) {
    itk::ImageRegionConstIteratorWithIndex<T> labelIter(label, label->GetBufferedRegion());
    for (labelIter.GoToBegin(); !labelIter.IsAtEnd(); ++labelIter) {
        typename T::PixelType label = labelIter.Get();
        if (label > 0) {
            container.push_back(labelIter.GetIndex());
        }
    }
}
}
#endif
