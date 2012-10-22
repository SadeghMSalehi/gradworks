project(pnsc++)
cmake_minimum_required(VERSION 2.8)

find_package(ITK REQUIRED)
find_package(VTK REQUIRED)
find_package(Qt4 REQUIRED)
find_package(Armadillo REQUIRED)

include(${ITK_USE_FILE})
include(${VTK_USE_FILE})
include(${QT_USE_FILE})

set (circleVis_SOURCES main.cpp mainwindow.h mainwindow.cpp PNSBase.h PNSMath.cpp PNSMath.h PNSCostFunction.h PNSCostFunction.cpp itkMyFRPROptimizer.h itkMyFRPROptimizer.cpp itkMyPowellOptimizer.h itkMyPowellOptimizer.cpp ExpLogTransform.h ExpLogTransform.cpp)
set (circleVis_FORMS mainwindow.ui)
set (circleVis_QT_HEADERS mainwindow.h)

qt4_wrap_cpp(circleVis_HEADERS_MOC ${circleVis_QT_HEADERS})
qt4_wrap_ui(circleVis_FORMS_HEADERS ${circleVis_FORMS})

include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories(/tools/Armadillo/armadillo-3.4.3/include)
#set (ARMADILLO_LIBRARIES /tools/Armadillo/Build/libarmadillo.so)

set(MathCode_INCLUDE /tools/gradworks/MathCode/include CACHE TYPE PATH)
set(MathCode_SRC /tools/gradworks/MathCode/src CACHE TYPE PATH)
set(MathCode_FILES ${MathCode_SRC}/MatrixCode.cpp ${MathCode_INCLUDE}/MatrixCode.h)

#set(NEWMAT_SRC newmat/newmat1.cpp newmat/newmat2.cpp newmat/newmat3.cpp newmat/newmat4.cpp newmat/newmat5.cpp newmat/newmat6.cpp newmat/newmat7.cpp newmat/newmat8.cpp newmat/newmat9.cpp newmat/bandmat.cpp newmat/submat.cpp newmat/myexcept.cpp newmat/newmatex.cpp)
#set(NEWRAN_SRC newran/newran.cpp newran/extreal.cpp)
#include_directories(newmat)
#include_directories(newran)


include_directories(${MathCode_INCLUDE})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_definitions(${QT_DEFINITIONS})

add_executable(circleVis
	${MathCode_FILES}
	${circleVis_SOURCES}
	${circleVis_FORMS_HEADERS}
	${circleVis_HEADERS_MOC}
)
target_link_libraries(circleVis ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${QT_LIBRARIES} ${ARMADILLO_LIBRARIES})
