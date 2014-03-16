CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(opencvapp)

#FIND_PACKAGE(OPENCV REQUIRED)

INCLUDE_DIRECTORIES(/opt/local/include)
LINK_DIRECTORIES(/opt/local/lib)
SET(OpenCV_LIBS opencv_highgui opencv_core opencv_ml opencv_video opencv_calib3d)
ADD_EXECUTABLE(run1 run1.cpp)
TARGET_LINK_LIBRARIES(run1 ${OpenCV_LIBS})

