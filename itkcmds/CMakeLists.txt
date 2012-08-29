project(itkcmds)
cmake_minimum_required(VERSION 2.8)

set(itkModules ITKCommon ITKIONRRD ITKIONIFTI)
find_package(ITK 4.2 COMPONENTS ${itkModules})
itk_module_config(itkinfo ${itkModules})
itk_module_use(${itkModules})

add_executable(itkinfo itkinfo.cpp)
target_link_libraries(itkinfo ${itkinfo_LIBRARIES})
