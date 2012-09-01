project(itkcmds)
cmake_minimum_required(VERSION 2.8)

set(itkModules ITKCommon ITKIONRRD ITKIOPNG ITKIONIFTI ITKImageFilterBase ITKRegistrationCommon)
find_package(ITK 4.2 COMPONENTS ${itkModules})

itk_module_config(itkcmds ${itkModules})
itk_module_use(${itkModules})

add_executable(itkinfo itkinfo.cpp)
target_link_libraries(itkinfo ${itkcmds_LIBRARIES})

add_executable(itkRigidReg itkRigidReg.cpp)
target_link_libraries(itkRigidReg ${itkcmds_LIBRARIES})

add_executable(itkTransform itkTransform.cpp)
target_link_libraries(itkTransform ${itkcmds_LIBRARIES})