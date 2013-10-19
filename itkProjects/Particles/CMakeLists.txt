PROJECT(ParticleGuidedRegistration)
cmake_minimum_required(VERSION 2.8)

# Tell CMake to run moc when necessary:
set(CMAKE_AUTOMOC ON)

# As moc files are generated in the binary dir, tell CMake
# to always look for includes there:
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Check if environment variable QTDIR is set.
# Extra security for windows environment as well.
if (DEFINED ENV{QTDIR})
    set(CMAKE_PREFIX_PATH $ENV{QTDIR} ${CMAKE_PREFIX_PATH})
else ()
    message (FATAL_ERROR "ERROR: Environment variable QTDIR is not set. Please locate your Qt folder MY_QT5_DIR.")
endif ()


SET (BUILD_2D False CACHE BOOL "Compile 2d version")
SET (BUILD_PEVAL False CACHE BOOL "Compile peval for label evaluation")
SET (BUILD_GUI False CACHE BOOL "Compile GUI tools")
SET (BUILD_MBRAINS False CACHE BOOL "Compile mbrains tools")

#FIND_PACKAGE(PiLibs REQUIRED)

FIND_PACKAGE(ITK 4 REQUIRED)

# QT Library
# Older One
#FIND_PACKAGE(Qt REQUIRED COMPONENTS QtOpenGL)
FIND_PACKAGE(Qt5 REQUIRED COMPONENTS Widgets PrintSupport Concurrent )

FIND_PACKAGE(VTK COMPONENTS REQUIRED)
SET(VTK_REQUIRED_OBJCXX_FLAGS -fobjc-gc)

#FIND_PACKAGE(Armadillo REQUIRED)

SET(Boost_NO_BOOST_CMAKE 1)
FIND_PACKAGE(Boost REQUIRED)

IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ENDIF()

#IF(Aramadillo_FOUND)
#  INCLUDE_DIRECTORIES(${Armadillo_INCLUDE_DIRS})
#ENDIF()

SET(Armadillo_LIBRARIES '')

INCLUDE(${USE_ITK_FILE})
#SET(QT_USE_QTOPENGL TRUE)
#INCLUDE(${QT_USE_FILE})
INCLUDE(${USE_VTK_FILE})
#INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})

#SET (QT_PI_HEADERS particleViewerWindow.h particleAni.h myQGraphicsView.h)
#SET (QT_PI_SOURCES particleViewerWindow.cpp  particleAni.cpp myQGraphicsView.cpp)
#qt4_wrap_cpp(QT_PI_MOC_HEADERS ${QT_PI_HEADERS})
#qt4_wrap_ui(QT_PI_UI_FILES particleViewerWindow.ui particleAni.ui)

#SET (QT_PLUTO_HEADERS piPlutoWindow.h piPlutoCore.h QGraphicsZoomView.h)
#SET (QT_PLUTO_SOURCES piPlutoWindow.cpp piPlutoCore.cpp QGraphicsZoomView.cpp)
#qt4_wrap_cpp(QT_PLUTO_MOC_HEADERS ${QT_PLUTO_HEADERS})
#qt4_wrap_ui(QT_PLUTO_UI_FILES plutowindow.ui)

SET (QT_PSIM2_HEADERS piSimul2.h piSimulCore.h piGroupSimul.h myQGraphicsView.h QGraphicsRectWidget.h QGraphicsZoomView.h qcustomplot.h qfileacceptbutton.h)
SET (QT_PSIM2_SOURCES piSimul2.cpp piSimulCore.cpp piGroupSimul.cpp myQGraphicsView.cpp QGraphicsRectWidget.cpp QGraphicsZoomView.cpp qcustomplot.cpp qfileacceptbutton.cpp)
QT5_WRAP_UI(QT_PSIM2_UIS simul2d.ui)
QT5_ADD_RESOURCES(QT_PSIM2_RCC simul2d.qrc)
#qt4_wrap_cpp(QT_PSIM2_MOC_HEADERS ${QT_PSIM2_HEADERS})
#qt4_wrap_ui(QT_PSIM2_UI_FILES simul2d.ui)
#QT4_ADD_RESOURCES(QT_PSIM2_RCC simul2d.qrc)

#SET (QT_PVIEW_HEADERS pviewAIRWindow.h airSuperZoom.h qcustomplot.h qdualslider.h qgraphicsguideview.h  qgraphicscompositeimageitem.h qgraphicsdrawingimageitem.h qfileacceptbutton.h qgraphicsvolumeview.h airAlgorithmManager.h qdrawingframe.h piqTestModule.h)
#SET (QT_PVIEW_SOURCES pviewAIRWindow.cpp airSuperZoom.cpp qcustomplot.cpp qdualslider.cpp qgraphicsguideview.cpp qgraphicscompositeimageitem.cpp qgraphicsdrawingimageitem.cpp qfileacceptbutton.cpp qgraphicsvolumeview.cpp airAlgorithmManager.cpp qdrawingframe.cpp piqTestModule.cpp)
#QT4_WRAP_CPP(QT_PVIEW_MOC_HEADERS ${QT_PVIEW_HEADERS})
#QT4_WRAP_UI(QT_PVIEW_UI_FILES pviewAIRWindow.ui testDialog.ui)
#QT4_ADD_RESOURCES(QT_PVIEW_RCC pviewAIRWindow.qrc)

#SET (QT_BVIEW_HEADERS bigViewWindow.h myQGraphicsView.h qcustomplot.h qdualslider.h qgraphicsvolumeview.h simpleSliceViewer.h dualImageViewer.h)
#SET (QT_BVIEW_SOURCES bigViewWindow.cpp myQGraphicsView.cpp qcustomplot.cpp qdualslider.cpp qgraphicsvolumeview.cpp simpleSliceViewer.cpp dualImageViewer.cpp)
#QT4_WRAP_CPP(QT_BVIEW_MOC_HEADERS ${QT_BVIEW_HEADERS})
#QT4_WRAP_UI(QT_BVIEW_UI_FILES bigViewWindow.ui simpleSliceViewer.ui dualImageViewer.ui)
#QT4_ADD_RESOURCES(QT_BVIEW_RCC bigViewWindow.qrc simpleSliceViewer.qrc)

add_definitions(${QT_DEFINITIONS})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

SET (MY_Qt_EXTENSIONS QGraphicsGridItem.cpp QGraphicsGridItem.h)

SET (MY_PARTICLE_WORKS_SOURCES piTimer.h piTimer.cpp piParticle.h piParticle.cpp piParticleCore.cpp piParticleCore.h piImageDef.h piImageDef.cpp piParticleBSpline.cpp piParticleBSpline.h piParticleForces.h piParticleForces.cpp piParticleCollision.h piParticleCollision.cpp piOptions.cpp piOptions.h piImageProcessing.cpp piImageProcessing.h piParticleHydroDynamics.cpp piParticleHydroDynamics.h piParticleSystemSolver.h piParticleSystemSolver.cpp piParticleTrainer.h piParticleTrainer.cpp
    piParticleTrace.h piParticleTrace.cpp piParticleTools.h piParticleTools.cpp piEntropyComputer.h piEntropyComputer.cpp piImageHistogram.h piImageHistogram.cpp itkMultiResolutionPyramidImageFilter.h itkMultiResolutionPyramidImageFilter.hxx)

#SET (PI_LIBS_SOURCE itkImageIO.h piTimer.h piTimer.cpp piImageProcessing.cpp piImageProcessing.h piImageDef.h piOptions.h piOptions.cpp piVTK.cpp piVTK.h)

#FILE(GLOB ITK_SOURCES "itk*.h" "itk*.hxx" "itk*.cxx")

#IF (${BUILD_2D})
#ADD_EXECUTABLE(prun2 particlesRun.cpp ${MY_PARTICLE_WORKS_SOURCES})
#SET_TARGET_PROPERTIES(prun2 PROPERTIES DEBUG_POSTFIX "d")
#TARGET_LINK_LIBRARIES(prun2 ${ITK_LIBRARIES} ${Armadillo_LIBRARIES} ${VTK_LIBRARIES})
#ENDIF (${BUILD_2D})

#ADD_EXECUTABLE(prun3 particlesRun.cpp ${MY_PARTICLE_WORKS_SOURCES})
#SET_TARGET_PROPERTIES(prun3 PROPERTIES COMPILE_FLAGS "-DBATCH -DATTR_SIZE=17")
#SET_TARGET_PROPERTIES(prun3 PROPERTIES DEBUG_POSTFIX "d")
#TARGET_LINK_LIBRARIES(prun3 ${ITK_LIBRARIES} ${VTK_LIBRARIES})

#ADD_EXECUTABLE(ptest ptest.cpp ${MY_PARTICLE_WORKS_SOURCES})
#TARGET_LINK_LIBRARIES(ptest ${Armadillo_LIBRARIES})

#ADD_EXECUTABLE(ptool ptool.h ptool.cpp piImageDef.h piImageDef.cpp piImageIO.h piOptions.h piOptions.cpp piMacros.h)
#TARGET_LINK_LIBRARIES(ptool ${ITK_LIBRARIES})

#IF (${BUILD_PEVAL})
#	ADD_EXECUTABLE(peval peval.cpp piImageDef.h piImageProcessing.cpp piImageProcessing.h piOptions.h piOptions.cpp itkImageIO.h)
#	TARGET_LINK_LIBRARIES(peval ${ITK_LIBRARIES})
#ENDIF(${BUILD_PEVAL})

IF (${BUILD_GUI})
	ADD_EXECUTABLE(psim simul2main.cpp ${QT_PSIM2_HEADERS} ${QT_PSIM2_SOURCES} ${QT_PSIM2_RCC} piImageDef.h piImageIO.h piImageSlice.h piImageSlice.cpp piImagePatch.h piImagePatch.cpp piPatchTracking.h piPatchTracking.cpp piImageEntropyComputer.h piImageEntropyComputer.cpp piImageRegistration.h piImageRegistration.cpp piPatchCompare.cpp piPatchCompare.h ${MY_PARTICLE_WORKS_SOURCES} qutils.h qutils.cpp itkSIFTImageFilter.h itkSIFTImageFilter.cpp itkARGBColorFunction.h itkScalarToARGBColormapImageFilter.h QParticlesGraphicsItem.cpp QParticlesGraphicsItem.h QGraphicsDirectionItem.h QGraphicsDirectionItem.cpp QGraphicsImageItem.h QGraphicsImageItem.cpp QGraphicsParticleItems.h QGraphicsParticleItems.cpp QGraphicsPatchItem.h QGraphicsPatchItem.cpp QGraphicsEventItem.h QGraphicsEventItem.cpp itkEntropyImageToImageMetricv4.h itkEntropyImageToImageMetricv4.hxx itkEntropyImageToImageMetricv4GetValueAndDerivativeThreader.h itkEntropyImageToImageMetricv4GetValueAndDerivativeThreader.hxx itkEntropyImageToImageMetricv4HelperThreader.h itkEntropyImageToImageMetricv4HelperThreader.hxx ${QT_PSIM2_UIS})
#${QT_PSIM2_MOC_HEADERS} ${QT_PSIM2_UI_FILES}  ${MY_Qt_EXTENSIONS} ${QT_PSIM2_RCC})
	SET_TARGET_PROPERTIES(psim PROPERTIES COMPILE_FLAGS "-DDIMENSIONS=2 -DATTR_SIZE=17")
	TARGET_LINK_LIBRARIES(psim ${ITK_LIBRARIES} ${QT_LIBRARIES} ${VTK_LIBRARIES})
  QT5_USE_MODULES(psim Widgets PrintSupport Concurrent)
    
#    ADD_EXECUTABLE(pluto plutomain.cpp piPlutoCore.h piPlutoCore.cpp ${QT_PLUTO_HEADERS} ${QT_PLUTO_SOURCES} ${MY_PARTICLE_WORKS_SOURCES} piFitCurve.h piFitCurve.cpp piImagePatch.h piBSplineBasis.h piBSplineBasis.cpp piImageIO.h piOptions.cpp piOptions.h QGraphicsImageItem.h QGraphicsImageItem.cpp QGraphicsPolygonDrawingInteraction.h QGraphicsPolygonDrawingInteraction.cpp qutils.h qutils.cpp ${QT_PLUTO_MOC_HEADERS} ${QT_PLUTO_UI_FILES})
#	SET_TARGET_PROPERTIES(pluto PROPERTIES COMPILE_FLAGS "-DDIMENSIONS=2")
#    TARGET_LINK_LIBRARIES(pluto ${ITK_LIBRARIES} ${QT_LIBRARIES} ${VTK_LIBRARIES})

#	ADD_EXECUTABLE(airview pview.cpp piOptions.cpp piOptions.h piVTK.h piVTK.cpp
#        piImageDef.h piImageIO.h piMacros.h piParticle.h piParticle.cpp
#		piImagePatch.h piParticleTracker.h piParticleTracker.cpp piImageSlice.h piImageSlice.cpp piImageProcessing.cpp piImageProcessing.h
#        piFitCurve.h piFitCurve.cpp piParticleCore.h piParticleCore.cpp piParticleTrace.h piParticleTrace.cpp
#        piParticleForces.h piParticleForces.cpp piParticleSystemSolver.h piParticleSystemSolver.cpp 
#        piEntropyComputer.h piEntropyComputer.cpp piParticleCollision.h piParticleCollision.cpp 
#        piTimer.h piTimer.cpp piParticleBSpline.h piParticleBSpline.cpp piContourSystem.h piContourSystem.cpp
#        piBSplineBasis.h piBSplineBasis.cpp piImageHistogram.cpp piImageHistogram.h qtypedef.h qutils.h qutils.cpp
#        airImageAlgorithm.cpp airImageAlgorithm.h airCLI.h airCLI.cpp
#		${QT_PVIEW_HEADERS} ${QT_PVIEW_SOURCES} ${QT_PVIEW_MOC_HEADERS} ${QT_PVIEW_UI_FILES} ${QT_PVIEW_RCC})
#	TARGET_LINK_LIBRARIES(airview ${ITK_LIBRARIES} ${QT_LIBRARIES} ${VTK_LIBRARIES} QVTK)

#	ADD_EXECUTABLE(Lemon bigview.cpp piImageDef.h 
#		piImageIO.h piOptions.h piOptions.cpp piImageSlice.h piImageSlice.cpp piVolumeDisplay.h 
#		itkARGBSliceExtractImageFilter.h itkARGBSliceExtractImageFilter.cpp
#		qtypedef.h qutils.h qutils.cpp
#		${QT_BVIEW_HEADERS} ${QT_BVIEW_SOURCES} ${QT_BVIEW_MOC_HEADERS} ${QT_BVIEW_UI_FILES} ${QT_BVIEW_RCC})
#	TARGET_LINK_LIBRARIES(Lemon ${ITK_LIBRARIES} ${QT_LIBRARIES} ${VTK_LIBRARIES} QVTK)
#
##	ADD_EXECUTABLE(pani3 piParticleAniMain.cpp piParticleAniMain.h ${QT_PI_HEADERS} ${QT_PI_SOURCES} ${QT_PI_MOC_HEADERS} ${QT_PI_UI_FILES} ${MY_PARTICLE_WORKS_SOURCES} ${MY_Qt_EXTENSIONS} piVTK.cpp piVTK.h)
##	TARGET_LINK_LIBRARIES(pani3 ${ITK_LIBRARIES} ${QT_LIBRARIES} ${VTK_LIBRARIES} QVTK)
ENDIF (${BUILD_GUI})
#
#
#ADD_EXECUTABLE(sgen3 pwShapeGenerator.cpp pwShapeGenerator.h)
#SET_TARGET_PROPERTIES(sgen3 PROPERTIES COMPILE_FLAGS -DDIMENSION3)
#TARGET_LINK_LIBRARIES(sgen3 ${ITK_LIBRARIES})
