project(imageParticles)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_XCODE_ATTRIBUTE_GCC_VERSION "com.apple.compilers.llvmgcc42")

find_package(ITK 4 REQUIRED)
find_package(Qt4 REQUIRED)
find_package(Armadillo REQUIRED)

include(${ITK_USE_FILE})
include(${QT_USE_FILE})

include_directories(${ARMADILLO_INCLUDE_DIRS})

set (CMAKE_C_FLAGS -fopenmp)
set (CMAKE_CXX_FLAGS -fopenmp)

#set (IPOPT_LIB_DIR /tools/IpOpt/lib)
#set (IPOPT_LIBRARIES ${IPOPT_LIB_DIR}/libiomp5.dylib ${IPOPT_LIB_DIR}/libipopt.dylib ${IPOPT_LIB_DIR}/libpardiso412-MACOS-X86-64.dylib)
#include_directories(/tools/IpOpt/include)
#include_directories(/usr/llvm-gcc-4.2/lib/gcc/i686-apple-darwin11/4.2.1/include)
set (IPOPT_LIBRARIES )

set(ITKFILES 
	itkScalarToARGBColormapImageFilter.h 
	itkScalarToARGBColormapImageFilter.hxx 
	itkARGBColormapFunction.h
	itkMyFRPROptimizer.h
	itkMyFRPROptimizer.cpp
	itkMyPowellOptimizer.h
	itkMyPowellOptimizer.cpp
	itkMyRegularStepGradientDescentBaseOptimizer.h
	itkMyRegularStepGradientDescentBaseOptimizer.cxx
	itkMyRegularStepGradientDescentOptimizer.h
	itkMyRegularStepGradientDescentOptimizer.cxx
	itkSphereBoundedGradientDescentOptimizer.h
	itkSphereBoundedGradientDescentOptimizer.cxx
)

set(HEADERS mainwindow.h qcustomplot.h)

set(SOURCES main.cpp 
	mainwindow.cpp 
	itkImageIO.h 
	SurfaceEntropyCostFunction.h 
	SurfaceEntropyCostFunction.cpp 
	NBody.h 
	NBody.cpp 
	qcustomplot.cpp
)

set(FORMS particlesUI.ui)

set (ITK_PARTICLES_SOURCES 
	itkParticleAttribute.h
	itkParticleContainer.h
	itkParticleContainerArrayAttribute.h
	itkParticleDomain.h
	itkParticleDomain.hxx
	itkParticleEvents.h
	itkParticleFunctionBasedShapeSpaceData.h
	itkParticleGaussianModeWriter.h
	itkParticleGaussianModeWriter.hxx
	itkParticleEntropyGradientFunction.h
	itkParticleEntropyGradientFunction.hxx
	itkParticleGeneralEntropyGradientFunction.h
	itkParticleGeneralEntropyGradientFunction.hxx
	itkParticleGradientDescentPositionOptimizer.h
	itkParticleGradientDescentPositionOptimizer.hxx
	itkParticleRegionDomain.h
	itkParticleImageDomain.h
	itkParticleImageDomainWithGradients.h
	itkParticleNeighborhood.h
	itkParticleShapeMatrixAttribute.h
	itkParticleSystem.h
	itkParticleSystem.hxx
	itkParticleVectorFunction.h
)

set (IPOPT_FILES SurfaceEntropyNLP.h SurfaceEntropyNLP.cpp)

qt4_wrap_cpp(HEADERS_MOC ${HEADERS})
qt4_wrap_ui(FORMS_HEADERS ${FORMS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_definitions(${QT_DEFINITIONS})

add_executable(imageParticles 
#	${IPOPT_FILES}
	${ITKFILES}
	${HEADERS} 
	${SOURCES} 
	${HEADERS_MOC} 
	${FORMS_HEADERS})

target_link_libraries(imageParticles ${ITK_LIBRARIES} ${QT_LIBRARIES} ${ARMADILLO_LIBRARIES} ${IPOPT_LIBRARIES})
