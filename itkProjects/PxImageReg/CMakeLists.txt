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
	piTimer.h piTimer.cpp  
	piOptions.cpp piOptions.h 
	piConfigFile.h piConfigFile.cpp 
	piPx.h piPx.cpp
	piParticleRunner.h piParticleRunner.cpp 
	piParticleWarp.h piParticleWarp.cpp 
	piImageProc.h piImageProc.cpp
<<<<<<< HEAD
	piEntropyComputer.h piEntropyComputer.cpp
    
    piImageQA.h piImageQA.cpp
    piImageHistogram.h piImageHistogram.cpp
)
=======
	piEntropyComputer.h piEntropyComputer.cpp)
>>>>>>> c6853f35b3ff2a6bd8f827d1c47b42a6faf7460c

ADD_EXECUTABLE(prun2 particlesRun.cpp ${MY_PARTICLE_WORKS_SOURCES})
SET_TARGET_PROPERTIES(prun2 PROPERTIES DEBUG_POSTFIX "d")
SET_TARGET_PROPERTIES(prun2 PROPERTIES COMPILE_DEFINITIONS "DIMENSIONS=2;ATTR_SIZE=17;BATCH")
TARGET_LINK_LIBRARIES(prun2 ${ITK_LIBRARIES} ${VTK_LIBRARIES} nlopt ${LOCAL_DIR}/lib/libconfig++.a)
