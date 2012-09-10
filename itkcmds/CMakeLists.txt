project(itkcmds)
cmake_minimum_required(VERSION 2.8)

set(itkModules ITKCommon ITKIONRRD ITKIOPNG ITKIONIFTI ITKImageFilterBase ITKRegistrationCommon)
find_package(ITK 4.2 COMPONENTS ${itkModules})

itk_module_config(itkcmds ${itkModules})
itk_module_use(${itkModules})

add_executable(itkinfo itkinfo.cpp)
target_link_libraries(itkinfo ${itkcmds_LIBRARIES})

add_executable(itkImageReg itkImageCommon.h itkImageReg.cpp itkTimer.h itkTimer.cpp MatrixCode.cpp itkMyMetric.h itkMyMetric.hxx itkMyFRPROptimizer.h itkMyFRPROptimizer.cpp)
target_link_libraries(itkImageReg ${itkcmds_LIBRARIES})

add_executable(multiReg itkImageCommon.h itkMultiReg.cpp itkTimer.h itkTimer.cpp itkOptimizationReporter.h itkMetaMetrics.h itkMetaMetrics.cpp itkMyMetric.h itkMyMetric.hxx itkMyFRPROptimizer.h itkMyFRPROptimizer.cpp)
target_link_libraries(multiReg ${itkcmds_LIBRARIES})

add_executable(itkTransform itkTransform.cpp)
target_link_libraries(itkTransform ${itkcmds_LIBRARIES})
