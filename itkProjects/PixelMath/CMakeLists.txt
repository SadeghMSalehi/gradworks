PROJECT(ImageMath)

FIND_PACKAGE(ITK REQUIRED)

IF (ITK_FOUND)
 INCLUDE(${USE_ITK_FILE})
ELSE(ITK_FOUND)
 MESSAGE(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
ENDIF (ITK_FOUND)

set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE( Boost COMPONENTS program_options)

INCLUDE_DIRECTORIES(/opt/local/include)
SET (PixelMath_SRC 
muParserBase.cpp muParserBase.h muParserBytecode.cpp muParserBytecode.h muParserCallback.cpp muParserCallback.h muParser.cpp muParserDef.h muParserDLL.cpp muParserDLL.h muParserError.cpp muParserError.h muParserFixes.h muParser.h muParserInt.cpp muParserInt.h muParserStack.h muParserTemplateMagic.h muParserTest.cpp muParserTest.h muParserToken.h muParserTokenReader.cpp muParserTokenReader.h PixelMathImageFilter.h PixelMathImageFilter.hxx itkImageCommon.h)

ADD_EXECUTABLE( PixelMath PixelMath.cxx ${PixelMath_SRC})
TARGET_LINK_LIBRARIES( PixelMath ${ITK_LIBRARIES} ${Boost_LIBRARIES})
