#ifndef __itkARGBColormapFunction_h
#define __itkARGBColormapFunction_h

namespace itk
{
    namespace Function
    {

        /**
         * \class ARGBColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT ARGBColormapFunction:public Object
        {
        public:

            typedef ARGBColormapFunction            Self;
            typedef Object                      Superclass;
            typedef SmartPointer< Self >        Pointer;
            typedef SmartPointer< const Self >  ConstPointer;


            /** Run-time type information (and related methods). */
            itkTypeMacro(ColormapFunction, Object);

            typedef TRGBPixel                                      RGBPixelType;
            typedef typename TRGBPixel::ComponentType              RGBComponentType;
            typedef TScalar                                        ScalarType;
            typedef typename NumericTraits< ScalarType >::RealType RealType;

            RGBComponentType m_Alpha;

            itkSetMacro(MinimumRGBComponentValue, RGBComponentType);
            itkGetConstMacro(MinimumRGBComponentValue, RGBComponentType);

            itkSetMacro(MaximumRGBComponentValue, RGBComponentType);
            itkGetConstMacro(MaximumRGBComponentValue, RGBComponentType);

            itkSetMacro(Alpha, RGBComponentType);
            itkGetConstMacro(Alpha, RGBComponentType);
            
            itkSetMacro(MinimumInputValue, ScalarType);
            itkGetConstMacro(MinimumInputValue, ScalarType);

            itkSetMacro(MaximumInputValue, ScalarType);
            itkGetConstMacro(MaximumInputValue, ScalarType);

            virtual bool operator!=(const ARGBColormapFunction &) const
            {
                return false;
            }

            virtual bool operator==(const ARGBColormapFunction & other) const
            {
                return !( *this != other );
            }

            virtual RGBPixelType operator()(const ScalarType &) const = 0;

        protected:
            
            ARGBColormapFunction()
            {
                this->m_MinimumInputValue = NumericTraits< TScalar >::min();
                this->m_MaximumInputValue = NumericTraits< TScalar >::max();
                this->m_MinimumRGBComponentValue = NumericTraits< RGBComponentType >::min();
                this->m_MaximumRGBComponentValue = NumericTraits< RGBComponentType >::max();
                this->m_Alpha = 255;
            }

            ~ARGBColormapFunction() {}

            /**
             * Map [min, max] input values to [0, 1].
             */
            RealType RescaleInputValue(ScalarType v) const
            {
                RealType d = static_cast< RealType >( this->m_MaximumInputValue -
                                                     this->m_MinimumInputValue );
                RealType value = ( static_cast< RealType >( v ) -
                                  static_cast< RealType >( this->m_MinimumInputValue ) ) / d;

                value = vnl_math_max(0.0, value);
                value = vnl_math_min(1.0, value);
                return value;
            }

            /**
             * Map [0, 1] value to [min, max] rgb component values.
             */
            RGBComponentType RescaleRGBComponentValue(RealType v) const
            {
                RealType d = static_cast< RealType >( m_MaximumRGBComponentValue -
                                                     m_MinimumRGBComponentValue );
                const RGBComponentType rescaled = static_cast< RGBComponentType >(
                                                                                  d * v ) + this->m_MinimumRGBComponentValue;

                return rescaled;
            }

            void PrintSelf(std::ostream & os, Indent indent) const
            {
                Superclass::PrintSelf(os, indent);

                os << indent << "Minimum RGB Component Value: "
                << static_cast< typename NumericTraits< RGBComponentType >::PrintType >(
                                                                                        this->GetMinimumRGBComponentValue() ) << std::endl;
                os << indent << "Maximum RGB Component Value: "
                << static_cast< typename NumericTraits< RGBComponentType >::PrintType >(
                                                                                        this->GetMaximumRGBComponentValue() ) << std::endl;
                os << indent << "Minimum Input Value: "
                << static_cast< typename NumericTraits< ScalarType >::PrintType >(
                                                                                  this->GetMinimumInputValue() ) << std::endl;
                os << indent << "Maximum Input Value: "
                << static_cast< typename NumericTraits< ScalarType >::PrintType >(
                                                                                  this->GetMaximumInputValue() ) << std::endl;
            }

        private:
            ARGBColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);  //purposely not implemented

            ScalarType m_MinimumInputValue;
            ScalarType m_MaximumInputValue;

            RGBComponentType m_MinimumRGBComponentValue;
            RGBComponentType m_MaximumRGBComponentValue;

        };

        /**
         * \class AutumnColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * \image html AutumnColormapFunction.png "Autumn colormap."
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT AutumnColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:

            typedef AutumnColormapFunction                 Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;

            /** Method for creation through the object factory. */
            itkNewMacro(Self);

            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;

            virtual RGBPixelType operator()(const TScalar &) const;

        protected:
            AutumnColormapFunction() {}
            ~AutumnColormapFunction() {}
        private:
            AutumnColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);        //purposely not implemented
        };

        template< class TScalar, class TRGBPixel >
        typename AutumnColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        AutumnColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);

            // Apply the color mapping.
            RealType red = 1.0;
            RealType green = value;
            RealType blue = 0.0;

            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);

            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;
            pixel[0] = this->RescaleRGBComponentValue(red);
            pixel[1] = this->RescaleRGBComponentValue(green);
            pixel[2] = this->RescaleRGBComponentValue(blue);

            return pixel;
        }

        /**
         * \class BlueColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * \image html BlueColormapFunction.png "Blue colormap."
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT BlueColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:

            typedef BlueColormapFunction                   Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;

            /** Method for creation through the object factory. */
            itkNewMacro(Self);

            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;

            virtual RGBPixelType operator()(const TScalar &) const;

        protected:
            BlueColormapFunction() {}
            ~BlueColormapFunction() {}
        private:
            BlueColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);      //purposely not implemented
        };

        template< class TScalar, class TRGBPixel >
        typename BlueColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        BlueColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);

            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);

            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;
            pixel[0] = 0;
            pixel[1] = 0;
            pixel[2] = this->RescaleRGBComponentValue(value);

            return pixel;
        }

        /**
         * \class CoolColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * \image html CoolColormapFunction.png "Cool colormap."
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT CoolColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:

            typedef CoolColormapFunction                   Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;

            /** Method for creation through the object factory. */
            itkNewMacro(Self);

            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;

            virtual RGBPixelType operator()(const TScalar &) const;

        protected:
            CoolColormapFunction() {}
            ~CoolColormapFunction() {}
        private:
            CoolColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);      //purposely not implemented
        };

        template< class TScalar, class TRGBPixel >
        typename CoolColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        CoolColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);

            // Apply the color mapping.
            RealType red = value;

            RealType green = 1.0 - value;

            RealType blue = 1.0;

            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);

            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;
            pixel[0] = this->RescaleRGBComponentValue(red);
            pixel[1] = this->RescaleRGBComponentValue(green);
            pixel[2] = this->RescaleRGBComponentValue(blue);

            return pixel;
        }

        /**
         * \class CopperColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \image html CopperColormapFunction.png "Copper colormap."
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT CopperColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:

            typedef CopperColormapFunction                 Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;

            /** Method for creation through the object factory. */
            itkNewMacro(Self);

            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;

            virtual RGBPixelType operator()(const TScalar &) const;

        protected:
            CopperColormapFunction() {}
            ~CopperColormapFunction() {}
        private:
            CopperColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);        //purposely not implemented
        };

        template< class TScalar, class TRGBPixel >
        typename CopperColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        CopperColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);

            // Apply the color map.
            RealType red = 1.2 * value;

            red = vnl_math_min(1.0, red);

            RealType green = 0.8 * value;

            RealType blue = 0.5 * value;

            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);

            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;
            pixel[0] = this->RescaleRGBComponentValue(red);
            pixel[1] = this->RescaleRGBComponentValue(green);
            pixel[2] = this->RescaleRGBComponentValue(blue);

            return pixel;
        }

        /**
         * \class CustomColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT CustomColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:

            typedef CustomColormapFunction                 Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;

            /** Method for creation through the object factory. */
            itkNewMacro(Self);

            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;

            typedef std::vector< RealType > ChannelType;

            virtual RGBPixelType operator()(const TScalar &) const;

            void SetRedChannel(ChannelType red)
            {
                m_RedChannel = red;
            }

            ChannelType GetRedChannel() const
            {
                return m_RedChannel;
            }

            void SetGreenChannel(ChannelType green)
            {
                m_GreenChannel = green;
            }

            ChannelType GetGreenChannel() const
            {
                return m_GreenChannel;
            }

            void SetBlueChannel(ChannelType blue)
            {
                m_BlueChannel = blue;
            }

            ChannelType GetBlueChannel() const
            {
                return m_BlueChannel;
            }

        protected:
            CustomColormapFunction() {}
            ~CustomColormapFunction() {}
        private:
            CustomColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);        //purposely not implemented

            ChannelType m_RedChannel;
            ChannelType m_GreenChannel;
            ChannelType m_BlueChannel;
        };

        template< class TScalar, class TRGBPixel >
        typename CustomColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        CustomColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);

            // Apply the color mapping.
            RealType red = 0.0;

            if ( this->m_RedChannel.size() == 1 || value == 0.0 )
            {
                red = this->m_RedChannel[0];
            }
            else if ( this->m_RedChannel.size() > 1 )
            {
                RealType     size = static_cast< RealType >( this->m_RedChannel.size() );
                unsigned int index = Math::Ceil< unsigned int >( value * ( size - 1.0 ) );
                RealType     p1 = this->m_RedChannel[index];
                RealType     m1 = this->m_RedChannel[index - 1u];
                RealType     d = p1 - m1;
                red = d * ( size - 1.0 ) * ( value - ( index - 1.0 ) / ( size - 1.0 ) )
                + m1;
            }

            RealType green = 0.0;
            if ( this->m_GreenChannel.size() == 1 || value == 0.0 )
            {
                green = this->m_GreenChannel[0];
            }
            else if ( this->m_GreenChannel.size() > 1 )
            {
                RealType     size = static_cast< RealType >( this->m_GreenChannel.size() );
                unsigned int index = Math::Ceil< unsigned int >( value * ( size - 1.0 ) );
                RealType     p1 = this->m_GreenChannel[index];
                RealType     m1 = this->m_GreenChannel[index - 1u];
                RealType     d = p1 - m1;
                green = d * ( size - 1.0 ) * ( value - ( index - 1.0 ) / ( size - 1.0 ) )
                + m1;
            }

            RealType blue = 0.0;
            if ( this->m_BlueChannel.size() == 1 || value == 0.0 )
            {
                blue = this->m_BlueChannel[0];
            }
            else if ( this->m_BlueChannel.size() > 1 )
            {
                RealType     size = static_cast< RealType >( this->m_BlueChannel.size() );
                unsigned int index = Math::Ceil< unsigned int >( value * ( size - 1.0 ) );
                RealType     p1 = this->m_BlueChannel[index];
                RealType     m1 = this->m_BlueChannel[index - 1u];
                RealType     d = p1 - m1;
                blue = d * ( size - 1.0 ) * ( value - ( index - 1.0 ) / ( size - 1.0 ) )
                + m1;
            }

            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);

            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;
            pixel[0] = this->RescaleRGBComponentValue(red);
            pixel[1] = this->RescaleRGBComponentValue(green);
            pixel[2] = this->RescaleRGBComponentValue(blue);

            return pixel;
        }

        /**
         * \class GreenColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \image html GreenColormapFunction.png "Green colormap."
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT GreenColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:

            typedef GreenColormapFunction                  Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;

            /** Method for creation through the object factory. */
            itkNewMacro(Self);

            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;

            virtual RGBPixelType operator()(const TScalar &) const;

        protected:
            GreenColormapFunction() {}
            ~GreenColormapFunction() {}
        private:
            GreenColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);       //purposely not implemented
        };

        template< class TScalar, class TRGBPixel >
        typename GreenColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        GreenColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);

            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);

            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;
            pixel[0] = 0;
            pixel[1] = this->RescaleRGBComponentValue(value);
            pixel[2] = 0;

            return pixel;
        }

        /**
         * \class GreyColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \image html GreyColormapFunction.png "Grey colormap."
         *
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT GreyColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:

            typedef GreyColormapFunction                   Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;

            /** Method for creation through the object factory. */
            itkNewMacro(Self);

            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;

            virtual RGBPixelType operator()(const TScalar &) const;

        protected:
            GreyColormapFunction() {}
            ~GreyColormapFunction() {}
        private:
            GreyColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);      //purposely not implemented
        };

        template< class TScalar, class TRGBPixel >
        typename GreyColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        GreyColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);

            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);

            pixel[0] = this->RescaleRGBComponentValue(value);
            pixel[1] = pixel[0];
            pixel[2] = pixel[0];
            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;
            return pixel;
        }

        /**
         * \class HSVColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \image html HSVColormapFunction.png "HSV colormap."
         *
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT HSVColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:

            typedef HSVColormapFunction                    Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;

            /** Method for creation through the object factory. */
            itkNewMacro(Self);

            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;

            virtual RGBPixelType operator()(const TScalar &) const;

        protected:
            HSVColormapFunction() {}
            ~HSVColormapFunction() {}
        private:
            HSVColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);     //purposely not implemented
        };

        template< class TScalar, class TRGBPixel >
        typename HSVColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        HSVColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);

            // Apply the color mapping.
            // Apply the color mapping.
            RealType red = vnl_math_abs( 5.0 * ( value - 0.5 ) ) - 5.0 / 6.0;

            red = vnl_math_min(red, 1.0);
            red = vnl_math_max(0.0, red);

            RealType green = -vnl_math_abs( 5.0 * ( value - 11.0 / 30.0 ) ) + 11.0 / 6.0;
            green = vnl_math_min(green, 1.0);
            green = vnl_math_max(0.0, green);

            RealType blue = -vnl_math_abs( 5.0 * ( value - 19.0 / 30.0 ) ) + 11.0 / 6.0;
            blue = vnl_math_min(blue, 1.0);
            blue = vnl_math_max(0.0, blue);

            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);

            pixel[0] = this->RescaleRGBComponentValue(red);
            pixel[1] = this->RescaleRGBComponentValue(green);
            pixel[2] = this->RescaleRGBComponentValue(blue);
            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;

            return pixel;
        }

        /**
         * \class HotColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \image html HotColormapFunction.png "Hot colormap."
         *
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT HotColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:

            typedef HotColormapFunction                    Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;

            /** Method for creation through the object factory. */
            itkNewMacro(Self);

            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;

            virtual RGBPixelType operator()(const TScalar &) const;

        protected:
            HotColormapFunction() {}
            ~HotColormapFunction() {}
        private:
            HotColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);     //purposely not implemented
        };

        template< class TScalar, class TRGBPixel >
        typename HotColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        HotColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);

            // Apply the color mapping.
            RealType red   = 63.0 / 26.0 * value - 1.0 / 13.0;

            red = vnl_math_max(0.0, red);
            red = vnl_math_min(1.0, red);

            RealType green = 63.0 / 26.0 * value - 11.0 / 13.0;
            green = vnl_math_max(0.0, green);
            green = vnl_math_min(1.0, green);

            RealType blue  = 4.5 * value - 3.5;
            blue = vnl_math_max(0.0, blue);
            blue = vnl_math_min(1.0, blue);

            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);

            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;
            pixel[0] = this->RescaleRGBComponentValue(red);
            pixel[1] = this->RescaleRGBComponentValue(green);
            pixel[2] = this->RescaleRGBComponentValue(blue);

            return pixel;
        }

        /**
         * \class JetColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \image html JetColormapFunction.png "Jet colormap."
         *
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT JetColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:

            typedef JetColormapFunction                    Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;

            /** Method for creation through the object factory. */
            itkNewMacro(Self);

            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;

            virtual RGBPixelType operator()(const TScalar &) const;

        protected:
            JetColormapFunction() {}
            ~JetColormapFunction() {}
        private:
            JetColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);     //purposely not implemented
        };

        template< class TScalar, class TRGBPixel >
        typename JetColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        JetColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);

            // Apply the color mapping.
            RealType red = -vnl_math_abs( 3.95 * ( value - 0.7460 ) ) + 1.5;

            red = vnl_math_min(red, 1.0);
            red = vnl_math_max(0.0, red);

            RealType green = -vnl_math_abs( 3.95 * ( value - 0.492 ) ) + 1.5;
            green = vnl_math_min(green, 1.0);
            green = vnl_math_max(0.0, green);

            RealType blue = -vnl_math_abs( 3.95 * ( value - 0.2385 ) ) + 1.5;
            blue = vnl_math_min(blue, 1.0);
            blue = vnl_math_max(0.0, blue);

            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);

            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;
            pixel[0] = this->RescaleRGBComponentValue(red);
            pixel[1] = this->RescaleRGBComponentValue(green);
            pixel[2] = this->RescaleRGBComponentValue(blue);

            return pixel;
        }

        /**
         * \class OverUnderColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \image html OverUnderColormapFunction.png "OverUnder colormap."
         *
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT OverUnderColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:

            typedef OverUnderColormapFunction              Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;

            /** Method for creation through the object factory. */
            itkNewMacro(Self);

            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;

            virtual RGBPixelType operator()(const TScalar &) const;

        protected:
            OverUnderColormapFunction() {}
            ~OverUnderColormapFunction() {}
        private:
            OverUnderColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);           //purposely not implemented
        };

        template< class TScalar, class TRGBPixel >
        typename OverUnderColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        OverUnderColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);

            // Apply the color mapping.
            RealType red = value;
            RealType green = value;
            RealType blue = value;

            if ( value == 0.0 )
            {
                // pixel is saturated in the dark
                red = 0.0;
                green = 0.0;
                blue = 1.0;
            }
            else if ( value == 1.0 )
            {
                // pixel is saturated in the white
                red = 1.0;
                green = 0.0;
                blue = 0.0;
            }

            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);
            
            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;
            pixel[0] = this->RescaleRGBComponentValue(red);
            pixel[1] = this->RescaleRGBComponentValue(green);
            pixel[2] = this->RescaleRGBComponentValue(blue);
            
            return pixel;
        }
        
        /**
         * \class RedColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \image html RedColormapFunction.png "Red colormap."
         *
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT RedColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:
            
            typedef RedColormapFunction                    Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;
            
            /** Method for creation through the object factory. */
            itkNewMacro(Self);
            
            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;
            
            virtual RGBPixelType operator()(const TScalar &) const;
            
        protected:
            RedColormapFunction() {}
            ~RedColormapFunction() {}
        private:
            RedColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);     //purposely not implemented
        };
        
        template< class TScalar, class TRGBPixel >
        typename RedColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        RedColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);
            
            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);
            
            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;
            pixel[0] = this->RescaleRGBComponentValue(value);
            pixel[1] = 0;
            pixel[2] = 0;
            
            return pixel;
        }
        
        /**
         * \class SpringColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \image html SpringColormapFunction.png "Spring colormap."
         *
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT SpringColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:
            
            typedef SpringColormapFunction                 Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;
            
            /** Method for creation through the object factory. */
            itkNewMacro(Self);
            
            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;
            
            virtual RGBPixelType operator()(const TScalar &) const;
            
        protected:
            SpringColormapFunction() {}
            ~SpringColormapFunction() {}
        private:
            SpringColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);        //purposely not implemented
        };
        
        template< class TScalar, class TRGBPixel >
        typename SpringColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        SpringColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);
            
            // Apply the color mapping.
            RealType red = 1.0;
            
            RealType green = value;
            
            RealType blue = 1.0 - value;
            
            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);
            
            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;
            pixel[0] = this->RescaleRGBComponentValue(red);
            pixel[1] = this->RescaleRGBComponentValue(green);
            pixel[2] = this->RescaleRGBComponentValue(blue);
            
            return pixel;
        }
        
        /**
         * \class SummerColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \image html SummerColormapFunction.png "Summer colormap."
         *
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT SummerColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:
            
            typedef SummerColormapFunction                 Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;
            
            /** Method for creation through the object factory. */
            itkNewMacro(Self);
            
            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;
            
            virtual RGBPixelType operator()(const TScalar &) const;
            
        protected:
            SummerColormapFunction() {}
            ~SummerColormapFunction() {}
        private:
            SummerColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);        //purposely not implemented
        };
        
        template< class TScalar, class TRGBPixel >
        typename SummerColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        SummerColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);
            
            // Apply the color mapping.
            RealType red = value;
            
            RealType green = 0.5 * value + 0.5;
            
            RealType blue = 0.4;
            
            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);
            
            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;
            pixel[0] = this->RescaleRGBComponentValue(red);
            pixel[1] = this->RescaleRGBComponentValue(green);
            pixel[2] = this->RescaleRGBComponentValue(blue);
            
            return pixel;
        }
        
        /**
         * \class WinterColormapFunction
         * \brief Function object which maps a scalar value into an RGB colormap value.
         *
         * \image html WinterColormapFunction.png "Winter colormap."
         *
         * \author Nicholas Tustison, Hui Zhang, Gaetan Lehmann, Paul Yushkevich
         * and James C. Gee
         *
         * This code was contributed in the Insight Journal paper:
         *
         * "Meeting Andy Warhol Somewhere Over the Rainbow: RGB Colormapping and ITK"
         * http://www.insight-journal.org/browse/publication/285
         * http://hdl.handle.net/1926/1452
         *
         * \ingroup ITKColormap
         */
        template< class TScalar, class TRGBPixel >
        class ITK_EXPORT WinterColormapFunction:
        public ARGBColormapFunction< TScalar, TRGBPixel >
        {
        public:
            
            typedef WinterColormapFunction                 Self;
            typedef ARGBColormapFunction< TScalar, TRGBPixel > Superclass;
            typedef SmartPointer< Self >                   Pointer;
            typedef SmartPointer< const Self >             ConstPointer;
            
            /** Method for creation through the object factory. */
            itkNewMacro(Self);
            
            typedef typename Superclass::RGBPixelType RGBPixelType;
            typedef typename Superclass::ScalarType   ScalarType;
            typedef typename Superclass::RealType     RealType;
            
            virtual RGBPixelType operator()(const TScalar &) const;
            
        protected:
            WinterColormapFunction() {}
            ~WinterColormapFunction() {}
        private:
            WinterColormapFunction(const Self &); //purposely not implemented
            void operator=(const Self &);        //purposely not implemented
        };
        
        template< class TScalar, class TRGBPixel >
        typename WinterColormapFunction< TScalar, TRGBPixel >::RGBPixelType
        WinterColormapFunction< TScalar, TRGBPixel >
        ::operator()(const TScalar & v) const
        {
            // Map the input scalar between [0, 1].
            RealType value = this->RescaleInputValue(v);
            
            // Apply the color map.
            RealType red = 0.0;
            
            RealType green = value;
            
            RealType blue = 1.0 - 0.5 * value;
            
            // Set the rgb components after rescaling the values.
            RGBPixelType pixel;
            NumericTraits<TRGBPixel>::SetLength(pixel, 4);
            
            pixel[3] = ARGBColormapFunction<TScalar, TRGBPixel>::m_Alpha;
            pixel[0] = this->RescaleRGBComponentValue(red);
            pixel[1] = this->RescaleRGBComponentValue(green);
            pixel[2] = this->RescaleRGBComponentValue(blue);
            
            return pixel;
        }
    } // end namespace Function
} // end namespace itk

#endif
