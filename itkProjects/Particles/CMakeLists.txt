PROJECT(ParticleGuidedRegistration)
cmake_minimum_required(VERSION 2.8)

SET (BUILD_2D False CACHE BOOL "Compile 2d version")
SET (BUILD_PEVAL False CACHE BOOL "Compile peval for label evaluation")
SET (BUILD_BIG False CACHE BOOL "Compile for distribution")
SET (BUILD_GUI False CACHE BOOL "Compile GUI tools")
SET (BUILD_MBRAINS False CACHE BOOL "Compile mbrains tools")

#FIND_PACKAGE(PiLibs REQUIRED)

FIND_PACKAGE(ITK 4 REQUIRED)
FIND_PACKAGE(Qt REQUIRED COMPONENTS QtOpenGL)

FIND_PACKAGE(VTK COMPONENTS REQUIRED)
SET(VTK_REQUIRED_OBJCXX_FLAGS -fobjc-gc)

FIND_PACKAGE(Armadillo REQUIRED)

SET(Boost_NO_BOOST_CMAKE 1)
FIND_PACKAGE(Boost REQUIRED)

IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ENDIF()

IF(Aramadillo_FOUND)
  INCLUDE_DIRECTORIES(${Armadillo_INCLUDE_DIRS})
ENDIF()

INCLUDE(${USE_ITK_FILE})
SET(QT_USE_QTOPENGL TRUE)
INCLUDE(${QT_USE_FILE})
INCLUDE(${USE_VTK_FILE})
INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})

SET (QT_PI_HEADERS particleViewerWindow.h particleAni.h myQGraphicsView.h)
SET (QT_PI_SOURCES particleViewerWindow.cpp  particleAni.cpp myQGraphicsView.cpp)
qt4_wrap_cpp(QT_PI_MOC_HEADERS ${QT_PI_HEADERS})
qt4_wrap_ui(QT_PI_UI_FILES particleViewerWindow.ui particleAni.ui)

SET (QT_PSIM2_HEADERS piSimul2.h myQGraphicsView.h qcustomplot.h)
SET (QT_PSIM2_SOURCES piSimul2.cpp myQGraphicsView.cpp qcustomplot.cpp)
qt4_wrap_cpp(QT_PSIM2_MOC_HEADERS ${QT_PSIM2_HEADERS})
qt4_wrap_ui(QT_PSIM2_UI_FILES simul2d.ui)

SET (QT_PVIEW_HEADERS pviewAIRWindow.h airSuperZoom.h qcustomplot.h qdualslider.h qgraphicsguideview.h  qgraphicscompositeimageitem.h qgraphicsdrawingimageitem.h qfileacceptbutton.h qgraphicsvolumeview.h airAlgorithmManager.h)
SET (QT_PVIEW_SOURCES pviewAIRWindow.cpp airSuperZoom.cpp qcustomplot.cpp qdualslider.cpp qgraphicsguideview.cpp qgraphicscompositeimageitem.cpp qgraphicsdrawingimageitem.cpp qfileacceptbutton.cpp qgraphicsvolumeview.cpp airAlgorithmManager.cpp)
QT4_WRAP_CPP(QT_PVIEW_MOC_HEADERS ${QT_PVIEW_HEADERS})
QT4_WRAP_UI(QT_PVIEW_UI_FILES pviewAIRWindow.ui testDialog.ui)
QT4_ADD_RESOURCES(QT_PVIEW_RCC pviewAIRWindow.qrc)

add_definitions(${QT_DEFINITIONS})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

SET (MY_Qt_EXTENSIONS QGraphicsGridItem.cpp QGraphicsGridItem.h)

SET (MY_PARTICLE_WORKS_SOURCES itkImageIO.h piTimer.h piTimer.cpp piParticleCore.cpp piParticleCore.h piImageDef.h piParticleBSpline.cpp piParticleBSpline.h piParticleForces.h piParticleForces.cpp piParticleCollision.h piParticleCollision.cpp piOptions.cpp piOptions.h piImageProcessing.cpp piImageProcessing.h piParticleHydroDynamics.cpp piParticleHydroDynamics.h piParticleSystemSolver.h piParticleSystemSolver.cpp piParticleTrace.h piParticleTrace.cpp piParticleTools.h piParticleTools.cpp piVTK.cpp piVTK.h piEntropyComputer.h piEntropyComputer.cpp)

SET (PI_LIBS_SOURCE itkImageIO.h piTimer.h piTimer.cpp piImageProcessing.cpp piImageProcessing.h piImageDef.h piOptions.h piOptions.cpp piVTK.cpp piVTK.h)
SET (PI_VTK_SOURCES piVectorGrid.cpp piVectorGrid.h)

FILE(GLOB ITK_SOURCES "itk*.h" "itk*.hxx" "itk*.cxx")

IF (${BUILD_2D})
ADD_EXECUTABLE(prun2 particlesRun.cpp ${MY_PARTICLE_WORKS_SOURCES})
SET_TARGET_PROPERTIES(prun2 PROPERTIES DEBUG_POSTFIX "d")
TARGET_LINK_LIBRARIES(prun2 ${ITK_LIBRARIES} ${Armadillo_LIBRARIES} ${VTK_LIBRARIES})
ENDIF (${BUILD_2D})

ADD_EXECUTABLE(prun3 particlesRun.cpp ${MY_PARTICLE_WORKS_SOURCES})
SET_TARGET_PROPERTIES(prun3 PROPERTIES COMPILE_FLAGS "-DBATCH")
SET_TARGET_PROPERTIES(prun3 PROPERTIES DEBUG_POSTFIX "d")
TARGET_LINK_LIBRARIES(prun3 ${ITK_LIBRARIES} ${Armadillo_LIBRARIES} ${VTK_LIBRARIES})

ADD_EXECUTABLE(prun3a5 particlesRun.cpp ${MY_PARTICLE_WORKS_SOURCES})
SET_TARGET_PROPERTIES(prun3a5 PROPERTIES COMPILE_FLAGS "-DATTR_SIZE=5")
SET_TARGET_PROPERTIES(prun3a5 PROPERTIES DEBUG_POSTFIX "d")
TARGET_LINK_LIBRARIES(prun3a5 ${ITK_LIBRARIES} ${Armadillo_LIBRARIES} ${VTK_LIBRARIES})

ADD_EXECUTABLE(prun3a5q particlesRun.cpp ${MY_PARTICLE_WORKS_SOURCES})
SET_TARGET_PROPERTIES(prun3a5q PROPERTIES COMPILE_FLAGS "-DATTR_SIZE=5 -DBATCH")
SET_TARGET_PROPERTIES(prun3a5q PROPERTIES DEBUG_POSTFIX "d")
TARGET_LINK_LIBRARIES(prun3a5q ${ITK_LIBRARIES} ${Armadillo_LIBRARIES} ${VTK_LIBRARIES})

IF (${BUILD_BIG})
ADD_EXECUTABLE(prun3a17q particlesRun.cpp ${MY_PARTICLE_WORKS_SOURCES})
SET_TARGET_PROPERTIES(prun3a17q PROPERTIES COMPILE_FLAGS "-DATTR_SIZE=17 -DBATCH")
SET_TARGET_PROPERTIES(prun3a17q PROPERTIES DEBUG_POSTFIX "d")
TARGET_LINK_LIBRARIES(prun3a17q ${ITK_LIBRARIES} ${Armadillo_LIBRARIES} ${VTK_LIBRARIES})
ENDIF(${BUILD_BIG})

ADD_EXECUTABLE(ptest ptest.cpp ${MY_PARTICLE_WORKS_SOURCES})
TARGET_LINK_LIBRARIES(ptest ${ITK_LIBRARIES} ${Armadillo_LIBRARIES} ${VTK_LIBRARIES})

IF (${BUILD_PEVAL})
	ADD_EXECUTABLE(peval peval.cpp piImageDef.h piImageProcessing.cpp piImageProcessing.h piOptions.h piOptions.cpp itkImageIO.h)
	TARGET_LINK_LIBRARIES(peval ${ITK_LIBRARIES})
ENDIF(${BUILD_PEVAL})

IF (${BUILD_GUI})
	ADD_EXECUTABLE(psim simul2main.cpp piSimul2.cpp piSimul2.h piImageSlice.h piImageSlice.cpp itkARGBColorFunction.h itkScalarToARGBColormapImageFilter.h QParticlesGraphicsItem.cpp QParticlesGraphicsItem.h ${QT_PSIM2_HEADERS} ${QT_PSIM2_SOURCES} ${QT_PSIM2_MOC_HEADERS} ${QT_PSIM2_UI_FILES} ${MY_PARTICLE_WORKS_SOURCES} ${MY_Qt_EXTENSIONS})
	SET_TARGET_PROPERTIES(psim PROPERTIES COMPILE_FLAGS "-DDIMENSIONS=2")
	TARGET_LINK_LIBRARIES(psim ${ITK_LIBRARIES} ${QT_LIBRARIES} ${VTK_LIBRARIES} QVTK)

	ADD_EXECUTABLE(airview pview.cpp piOptions.cpp piOptions.h piVTK.h piVTK.cpp piImageDef.h 
		piImageIO.h piImageSlice.h piImageSlice.cpp piImageProcessing.cpp piImageProcessing.h 
		piImageHistogram.cpp piImageHistogram.h qtypedef.h qutils.h qutils.cpp
		${QT_PVIEW_HEADERS} ${QT_PVIEW_SOURCES} ${QT_PVIEW_MOC_HEADERS} ${QT_PVIEW_UI_FILES} ${QT_PVIEW_RCC})
	TARGET_LINK_LIBRARIES(airview ${ITK_LIBRARIES} ${QT_LIBRARIES} ${VTK_LIBRARIES} QVTK)

	ADD_EXECUTABLE(pani3 piParticleAniMain.cpp piParticleAniMain.h ${QT_PI_HEADERS} ${QT_PI_SOURCES} ${QT_PI_MOC_HEADERS} ${QT_PI_UI_FILES} ${MY_PARTICLE_WORKS_SOURCES} ${MY_Qt_EXTENSIONS} piVTK.cpp piVTK.h)
	TARGET_LINK_LIBRARIES(pani3 ${ITK_LIBRARIES} ${QT_LIBRARIES} ${VTK_LIBRARIES} QVTK)
ENDIF (${BUILD_GUI})


# 2D version executables
#ADD_EXECUTABLE(optionTest optionTestMain.cpp piOptions.h piOptions.cpp)

#ADD_EXECUTABLE(particleWin particleViewerWindowMain.cpp ${QT_PI_HEADERS} ${QT_PI_SOURCES} ${QT_PI_MOC_HEADERS} ${QT_PI_UI_FILES} ${MY_PARTICLE_WORKS_SOURCES} ${MY_Qt_EXTENSIONS})
#TARGET_LINK_LIBRARIES(particleWin ${ITK_LIBRARIES} ${QT_LIBRARIES} ${Armadillo_LIBRARIES} QVTK)

#ADD_EXECUTABLE(particleReg particleRegistrationMain.cpp ${MY_PARTICLE_WORKS_SOURCES} ${PI_LIBS_SOURCE})
#TARGET_LINK_LIBRARIES(particleReg ${ITK_LIBRARIES} ${Armadillo_LIBRARIES})

#ADD_EXECUTABLE(particleRun particlesRun.cpp ${MY_PARTICLE_WORKS_SOURCES})
#TARGET_LINK_LIBRARIES(particleRun ${ITK_LIBRARIES} ${Armadillo_LIBRARIES})


# 3D version executables
#ADD_EXECUTABLE(particleWin3 particleViewerWindowMain.cpp ${QT_PI_HEADERS} ${QT_PI_SOURCES} ${QT_PI_MOC_HEADERS} ${QT_PI_UI_FILES} ${MY_PARTICLE_WORKS_SOURCES} ${MY_Qt_EXTENSIONS})
#SET_TARGET_PROPERTIES(particleWin3 PROPERTIES COMPILE_FLAGS -DDIMENSION3)
#TARGET_LINK_LIBRARIES(particleWin3 ${ITK_LIBRARIES} ${QT_LIBRARIES} ${Armadillo_LIBRARIES} QVTK)

#ADD_EXECUTABLE(particleCollision3 particleCollision.cpp ${PI_LIBS_SOURCE} ${MY_PARTICLE_WORKS_SOURCES})
#SET_TARGET_PROPERTIES(particleCollision3 PROPERTIES COMPILE_FLAGS -DDIMENSION3)
#TARGET_LINK_LIBRARIES(particleCollision3 ${ITK_LIBRARIES} ${Armadillo_LIBRARIES})

#ADD_EXECUTABLE(particleConstraint3 particleConstraint.cpp ${PI_LIBS_SOURCE} ${MY_PARTICLE_WORKS_SOURCES} ${PI_VTK_SOURCES})
#SET_TARGET_PROPERTIES(particleConstraint3 PROPERTIES COMPILE_FLAGS -DDIMENSION3)
#TARGET_LINK_LIBRARIES(particleConstraint3 ${ITK_LIBRARIES} ${Armadillo_LIBRARIES} ${VTK_LIBRARIES})



#ADD_EXECUTABLE(preg3 particleRegistrationMain.cpp ${MY_PARTICLE_WORKS_SOURCES} ${PI_LIBS_SOURCE})
#SET_TARGET_PROPERTIES(preg3 PROPERTIES COMPILE_FLAGS -DDIMENSION3)
#TARGET_LINK_LIBRARIES(preg3 ${ITK_LIBRARIES} ${Armadillo_LIBRARIES})

ADD_EXECUTABLE(sgen3 pwShapeGenerator.cpp pwShapeGenerator.h)
SET_TARGET_PROPERTIES(sgen3 PROPERTIES COMPILE_FLAGS -DDIMENSION3)
TARGET_LINK_LIBRARIES(sgen3 ${ITK_LIBRARIES})

#ADD_EXECUTABLE(itkcmds itkcmdsMain.cpp ${MY_PARTICLE_WORKS_SOURCES})
#SET_TARGET_PROPERTIES(itkcmds PROPERTIES COMPILE_FLAGS -DDIMENSION3)
#TARGET_LINK_LIBRARIES(itkcmds ${ITK_LIBRARIES} ${Boost_LIBRARIES})
