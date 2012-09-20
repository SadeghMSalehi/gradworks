project(itktools)
cmake_minimum_required(VERSION 2.8)

find_package(ITK 4 REQUIRED)
find_package(Qt4 REQUIRED)

include(${ITK_USE_FILE})

set(itktools_SOURCES main.cpp mainwindow.cpp imageviewer.cpp itkImageIO.h)
set(itktools_QT_HEADERS mainwindow.h imageviewer.h)
set(itktools_FORMS mainwindow.ui imageviewer.ui)

qt4_wrap_cpp(itktools_HEADERS_MOC ${itktools_QT_HEADERS})
qt4_wrap_ui(itktools_FORMS_HEADERS ${itktools_FORMS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})

add_executable(main 
	${itktools_SOURCES} 
	${itktools_HEADERS_MOC} 
	${itktools_FORMS_HEADERS})

target_link_libraries(main ${ITK_LIBRARIES} ${QT_LIBRARIES})
