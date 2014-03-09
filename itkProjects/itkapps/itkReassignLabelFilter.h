#include <itkImage.h>
#include <itkUnaryFunctorImageFilter.h>

namespace itk {
    namespace Functor {
        template<class TInput, class TOutput>
        class Reassign {
        private:
          TInput _tLabel;
          TOutput _aLabel;
        public:
            Reassign() : _tLabel(0), _aLabel(0) {}
            ~Reassign() {}

            bool operator!=(const Reassign&) const {
                return false;
            }
            bool operator==(const Reassign& other) const {
                return !(*this != other);
            }
            inline TOutput operator()(const TInput& A) {
              std::cout << int(this) << std::endl;
              std::cout << A << " " << _tLabel << " " << _aLabel << " " << std::endl;
              if (A == _tLabel) {
                  return static_cast<TOutput>(_aLabel);
              } else {
                  return static_cast<TOutput>(A);
              }
            }
            void SetTargetLabel(TInput label) { 
              std::cout << int(this) << std::endl;
              _tLabel = label; }
            void SetAssignLabel(TOutput label) { _aLabel = label; }
        };
    }


template<class TInputImage, class TOutputImage>
    class ITK_EXPORT ReassignLabelFilter
        : public UnaryFunctorImageFilter<TInputImage, TOutputImage,
        Functor::Reassign<typename TInputImage::PixelType, typename TOutputImage::PixelType> >
    {
    public:
        typedef ReassignLabelFilter Self;
        typedef UnaryFunctorImageFilter<TInputImage, TOutputImage, 
                Functor::Reassign<typename TInputImage::PixelType, typename TOutputImage::PixelType> > Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro(Self);

    protected:
        ReassignLabelFilter() {}
        virtual ~ReassignLabelFilter() {}

    private:
        ReassignLabelFilter(const Self&);
        void operator=(const Self&);
    };

};
