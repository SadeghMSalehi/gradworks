CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(opencvapp)

INCLUDE_DIRECTORIES(/opt/local/include)
LINK_DIRECTORIES(/opt/local/lib)

SET(OPENCV_LIBS cxcore cv highgui cvaux ml)

ADD_EXECUTABLE(run1 run1.cpp)
TARGET_LINK_LIBRARIES(run1 ${OpenCV_LIBS})

