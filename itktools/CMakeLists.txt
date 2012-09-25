project(itktools)
cmake_minimum_required(VERSION 2.8)

find_package(ITK 4 REQUIRED)
find_package(Qt4 REQUIRED)

include(${ITK_USE_FILE})

set(itktools_MYCODE itkImageIO.h itkMyCore.h itkMyCore.cpp itkMyTypes.h
	itkTimer.h itkTimer.cpp itkMySlicer.h 
	MatrixCode.h MatrixCode.cpp itkMathCode.h itkMyRegistration.h
	itkMultiRegMethod.h itkMultiRegMethod.cpp	
)
set(itktools_SOURCES main.cpp mainwindow.cpp)
set(itktools_QT_HEADERS mainwindow.h)
set(itktools_FORMS mainwindow.ui)

set (ants_SRC antsAllocImage.h itkAvantsMutualInformationRegistrationFunction.cxx  
			itkAvantsMutualInformationRegistrationFunction.h    
			itkAvantsPDEDeformableRegistrationFunction.h )

set (itkcmds_SRC_METRICS itkMyMetric.h itkMyMetric.hxx itkMyMetricBase.h itkMyMetricBase.hxx itkMetaMetrics.h itkMetaMetrics.cpp)
set (itkcmds_SRC_OPTIMIZERS itkMyPowellOptimizer.h itkMyPowellOptimizer.cxx 
	itkMyFRPROptimizer.h itkMyFRPROptimizer.cpp itkOptimizationReporter.h 
	itkMyRegularStepGradientDescentOptimizer.h itkMyRegularStepGradientDescentOptimizer.cxx 
	itkMyRegularStepGradientDescentBaseOptimizer.h itkMyRegularStepGradientDescentBaseOptimizer.cxx 
	itkMyScaleVersor3DTransformOptimizer.h itkMyScaleVersor3DTransformOptimizer.cxx)

qt4_wrap_cpp(itktools_HEADERS_MOC ${itktools_QT_HEADERS})
qt4_wrap_ui(itktools_FORMS_HEADERS ${itktools_FORMS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})

add_executable(guireg 
	${ants_SRC}
	${itkcmds_SRC_METRICS}
	${itkcmds_SRC_OPTIMIZERS}
	${itktools_MYCODE}
	${itktools_QT_HEADERS} 
	${itktools_SOURCES} 
	${itktools_HEADERS_MOC} 
	${itktools_FORMS_HEADERS})

target_link_libraries(guireg ${ITK_LIBRARIES} ${QT_LIBRARIES})
