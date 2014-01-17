PROJECT(PxImageReg)
cmake_minimum_required(VERSION 2.8)

SET(LOCAL_DIR /tools/local CACHE PATH "local library directory")
FIND_PACKAGE(ITK 4 REQUIRED)
FIND_PACKAGE(VTK COMPONENTS REQUIRED)
SET(VTK_REQUIRED_OBJCXX_FLAGS -fobjc-gc)

SET(Boost_NO_BOOST_CMAKE 1)
FIND_PACKAGE(Boost REQUIRED)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ENDIF()

INCLUDE(${USE_ITK_FILE})
INCLUDE(${USE_VTK_FILE})

INCLUDE_DIRECTORIES(${LOCAL_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

LINK_DIRECTORIES(${LOCAL_DIR}/lib)


SET (MY_PARTICLE_WORKS_SOURCES
    itkARGBColorFunction.h
    itkScalarToARGBColormapImageFilter.h
    itkARGBSliceExtractImageFilter.h itkARGBSliceExtractImageFilter.cpp

    piImageDef.h piImageIO.h piMacros.h
    piTools.h piTools.cpp
    piTimer.h piTimer.cpp
    piOptions.cpp piOptions.h
    piConfigFile.h piConfigFile.cpp
    piPx.h piPx.cpp
    piParticleRunner.h piParticleRunner.cpp
    piParticleWarp.h piParticleWarp.cpp
    piImageProc.h piImageProc.cpp
    piEntropyComputer.h piEntropyComputer.cpp

    piImageQA.h piImageQA.cpp
    piImageHistogram.h piImageHistogram.cpp

# Label Fusion
    piLabelFusionRunner.h piLabelFusionRunner.cpp

# Demons Runner
    piDemonsRunner.h piDemonsRunner.cpp

# Metric Evaluation
    piImageRxRunner.h piImageRxRunner.cpp

# ImageToImageMetric
    itk/itkEntropyImageToImageMetricv4.h
    itk/itkEntropyImageToImageMetricv4.hxx
    itk/itkEntropyImageToImageMetricv4GetValueAndDerivativeThreader.h
    itk/itkEntropyImageToImageMetricv4GetValueAndDerivativeThreader.hxx
    itk/itkEntropyImageToImageMetricv4HelperThreader.h
    itk/itkEntropyImageToImageMetricv4HelperThreader.hxx
)

ADD_EXECUTABLE(prun2 particlesRun.cpp ${MY_PARTICLE_WORKS_SOURCES})
SET_TARGET_PROPERTIES(prun2 PROPERTIES DEBUG_POSTFIX "d")
SET_TARGET_PROPERTIES(prun2 PROPERTIES COMPILE_DEFINITIONS "DIMENSIONS=2;ATTR_SIZE=17;BATCH")
TARGET_LINK_LIBRARIES(prun2 ${ITK_LIBRARIES} ${VTK_LIBRARIES} nlopt ${LOCAL_DIR}/lib/libconfig++.a)

ADD_CUSTOM_COMMAND(TARGET prun2 POST_BUILD COMMAND prun2 --help > ${CMAKE_CURRENT_SOURCE_DIR}/../prun.md)

ADD_EXECUTABLE(prun3 particlesRun.cpp ${MY_PARTICLE_WORKS_SOURCES})
SET_TARGET_PROPERTIES(prun3 PROPERTIES DEBUG_POSTFIX "d")
SET_TARGET_PROPERTIES(prun3 PROPERTIES COMPILE_DEFINITIONS "DIMENSIONS=3;ATTR_SIZE=17;BATCH")
TARGET_LINK_LIBRARIES(prun3 ${ITK_LIBRARIES} ${VTK_LIBRARIES} nlopt ${LOCAL_DIR}/lib/libconfig++.a)
