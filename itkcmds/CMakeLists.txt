project(itkcmds)
cmake_minimum_required(VERSION 2.8)

if(OPENMP_FOUND)
 set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
 set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set (POPT_INCLUDE_DIR "/opt/local/include" CACHE PATH "libpopt include path")
set (POPT_LIB_DIR "/opt/local/lib" CACHE PATH "libpopt library path")

#set(itkModules ITKCommon ITKIONRRD ITKIOPNG ITKIONIFTI ITKImageFilterBase ITKRegistrationCommon ITKIOTransformBase ITKReview)
#find_package(ITK 4.2 COMPONENTS ${itkModules})
#itk_module_config(itkcmds ${itkModules})
#itk_module_use(${itkModules})

find_package(ITK 4 REQUIRED)
include(${ITK_USE_FILE})

set (itkcmds_LIBRARIES ${ITK_LIBRARIES})

set (SRC_METRICS itkMyMetric.h itkMyMetric.hxx itkMyMetricBase.h itkMyMetricBase.hxx)
set (SRC_OPTIMIZERS itkMyPowellOptimizer.h itkMyPowellOptimizer.cxx itkMyFRPROptimizer.h itkMyFRPROptimizer.cpp itkOptimizationReporter.h itkMyRegularStepGradientDescentOptimizer.h itkMyRegularStepGradientDescentOptimizer.cxx itkMyRegularStepGradientDescentBaseOptimizer.h itkMyRegularStepGradientDescentBaseOptimizer.cxx itkMyScaleVersor3DTransformOptimizer.h itkMyScaleVersor3DTransformOptimizer.cxx)
set (SRC_MYLIBS itkImageIO.h itkTimer.h itkTimer.cpp MatrixCode.h MatrixCode.cpp itkMathCode.h)
set (SRC_TRANSFORMS itkMyEulerTransform.h itkMyEulerTransform.hxx)
set (SRC_ALL ${SRC_METRICS} ${SRC_OPTIMIZERS} ${SRC_MYLIBS} ${SRC_TRANSFORMS})

#include_directories(${POPT_INCLUDE_DIR})
#link_directories(${POPT_LIB_DIR})
#add_executable(itkcmds itkcmds.cpp itkinfo.cpp itkImageIO.h ${POPT_INCLUDE_DIR}/popt.h)
#target_link_libraries(itkcmds ${itkcmds_LIBRARIES} popt)

#add_executable(imageRegistration itkImageIO.h imageRegistration.cpp ${SRC_ALL})
#target_link_libraries(imageRegistration ${itkcmds_LIBRARIES})

#add_executable(volumeRegistration itkImageIO.h volumeRegistration.cpp ${SRC_ALL})
#target_link_libraries(volumeRegistration ${itkcmds_LIBRARIES})

#add_executable(multiReg itkImageIO.h itkMultiReg.cpp itkMetaMetrics.h ${SRC_ALL})
#target_link_libraries(multiReg ${itkcmds_LIBRARIES})

add_executable(samplePixel samplePixel.cpp)
target_link_libraries(samplePixel ${ITK_LIBRARIES})
