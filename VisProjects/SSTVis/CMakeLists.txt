PROJECT(Vis)
cmake_minimum_required(VERSION 2.8)

FIND_PACKAGE(Qt4 REQUIRED)
FIND_PACKAGE(Armadillo REQUIRED)

INCLUDE(${QT_USE_FILE})

SET (QT_HEADERS mainwindow.h qcustomplot.h QGraphicsZoomView.h)
SET (QT_SOURCES mainwindow.cpp qcustomplot.cpp QGraphicsZoomView.cpp)
qt4_wrap_cpp(QT_MOC_HEADERS ${QT_HEADERS})
qt4_wrap_ui(QT_UI_FILES mainwindow.ui)

add_definitions(${QT_DEFINITIONS})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

ADD_EXECUTABLE(vis vis.cpp QGraphicsImageItem.h QGraphicsImageItem.cpp qutils.h qutils.cpp ${QT_HEADERS} ${QT_SOURCES} ${QT_MOC_HEADERS} ${QT_UI_FILES})
TARGET_LINK_LIBRARIES(vis ${QT_LIBRARIES} ${Armadillo_LIBRARIES})
