project(itkcmds)
cmake_minimum_required(VERSION 2.8)

if(OPENMP_FOUND)
 set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
 set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#set(itkModules ITKCommon ITKIONRRD ITKIOPNG ITKIONIFTI ITKImageFilterBase ITKRegistrationCommon ITKIOTransformBase ITKReview)
#find_package(ITK 4.2 COMPONENTS ${itkModules})
#itk_module_config(itkcmds ${itkModules})
#itk_module_use(${itkModules})

find_package(ITK 4 REQUIRED)
include(${ITK_USE_FILE})

set (itkcmds_LIBRARIES ${ITK_LIBRARIES})

add_executable(itkinfo itkinfo.cpp itkImageIO.h)
target_link_libraries(itkinfo ${itkcmds_LIBRARIES})

add_executable(itkImageReg itkImageIO.h itkImageReg.cpp itkTimer.h itkTimer.cpp MatrixCode.cpp itkMyMetric.h itkMyMetric.hxx itkMyFRPROptimizer.h itkMyFRPROptimizer.cpp itkOptimizationReporter.h itkMathCode.h)
target_link_libraries(itkImageReg ${itkcmds_LIBRARIES})

add_executable(multiReg itkImageIO.h itkMultiReg.cpp itkTimer.h itkTimer.cpp itkOptimizationReporter.h itkMetaMetrics.h itkMetaMetrics.cpp itkMyMetric.h itkMyMetric.hxx itkMyFRPROptimizer.h itkMyFRPROptimizer.cpp)
target_link_libraries(multiReg ${itkcmds_LIBRARIES})

add_executable(itkTransform itkTransform.cpp itkMathCode.h MatrixCode.cpp MatrixCode.h itkResampler.h)
target_link_libraries(itkTransform ${itkcmds_LIBRARIES})
