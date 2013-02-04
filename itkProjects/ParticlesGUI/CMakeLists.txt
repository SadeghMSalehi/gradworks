PROJECT(ParticlesGUI)
cmake_minimum_required(VERSION 2.8)

FIND_PACKAGE(ITK 4 REQUIRED)
FIND_PACKAGE(Qt4 REQUIRED)
FIND_PACKAGE(VTK REQUIRED)
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
INCLUDE(${QT_USE_FILE})
INCLUDE(${USE_VTK_FILE})
INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})

IF (NOT ITKVtkGlue_LOADED)
#	FIND_PACKAGE(ItkVtkGlue REQUIRED)
#	INCLUDE(${ItkVtkGlue_USE_FILE})
#	SET(Glue ItkVtkGlue)
ENDIF()

SET (QT_HEADERS mainwindow.h dialogPatchCompare.h qcustomplot.h myQGraphicsView.h dialogBSplineVis.h)
SET (QT_SOURCES mainwindow.cpp dialogPatchCompare.cpp qcustomplot.cpp myQGraphicsView.cpp dialogBSplineVis.cpp)
qt4_wrap_cpp(QT_MOC_HEADERS ${QT_HEADERS})
qt4_wrap_ui(QT_UI_FILES mainwindow.ui dialogPatchCompare.ui dialogBSplineVis.ui)

SET (QT_PW_HEADERS particleViewerWindow.h myQGraphicsView.h)
SET (QT_PW_SOURCES particleViewerWindow.cpp  myQGraphicsView.cpp)
qt4_wrap_cpp(QT_PW_MOC_HEADERS ${QT_PW_HEADERS})
qt4_wrap_ui(QT_PW_UI_FILES particleViewerWindow.ui)

add_definitions(${QT_DEFINITIONS})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

SET (MY_SOURCES myImageContainer.h myImageContainer.cpp myParticleAlgorithm.h myParticleAlgorithm.cpp vtkPropScene.h vtkPropScene.cpp PropertyAccess.h PropertyAccess.cpp myImageParticlesAlgorithm.h myImageParticlesAlgorithm.cpp myImplicitSurfaceConstraint.h myImplicitSurfaceConstraint.cpp myEnsembleEntropy.h myEnsembleEntropy.cpp vnlCommon.h vnlCommon.cpp myParticleDynamics.h myParticleDynamics.cpp myParticleSolver.h myParticleSolver.cpp myEventCallback.h myImageTransform.h myImageTransform.cpp myBSplineRegistration.h myBSplineRegistration.cpp QGraphicsGridItem.cpp QGraphicsGridItem.h myParticleKernelFunctions.h myParticleKernelFunctions.cpp myParticlesIO.cpp myParticlesIO.h)

SET (MY_Qt_EXTENSIONS QGraphicsGridItem.cpp QGraphicsGridItem.h)

SET (MY_PARTICLE_WORKS_SOURCES myParticleCore.cpp myParticleCore.h myParticleConstraint.cpp myParticleConstraint.h myImageDef.h myParticleBSpline.cpp myParticleBSpline.h myParticleForces.h myParticleForces.cpp)

FILE(GLOB ITK_SOURCES "itk*.h" "itk*.hxx" "itk*.cxx")

ADD_EXECUTABLE(guiparticles particlesMain.cpp ${MY_SOURCES} ${ITK_SOURCES} ${QT_HEADERS} ${QT_SOURCES} ${QT_MOC_HEADERS} ${QT_UI_FILES} ${MY_PARTICLE_WORKS_SOURCES})
TARGET_LINK_LIBRARIES(guiparticles ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${QT_LIBRARIES} ${Armadillo_LIBRARIES})

ADD_EXECUTABLE(particleWin particleViewerWindowMain.cpp ${QT_PW_HEADERS} ${QT_PW_SOURCES} ${QT_PW_MOC_HEADERS} ${QT_PW_UI_FILES} ${MY_PARTICLE_WORKS_SOURCES} ${MY_Qt_EXTENSIONS})
TARGET_LINK_LIBRARIES(particleWin ${ITK_LIBRARIES} ${QT_LIBRARIES} ${Armadillo_LIBRARIES})

ADD_EXECUTABLE(particlesRun particlesRun.cpp ${MY_PARTICLE_WORKS_SOURCES})
TARGET_LINK_LIBRARIES(particlesRun ${ITK_LIBRARIES} ${Armadillo_LIBRARIES})

ADD_EXECUTABLE(particleRegistration particleRegistrationMain.cpp ${MY_PARTICLE_WORKS_SOURCES})
TARGET_LINK_LIBRARIES(particleRegistration ${ITK_LIBRARIES} ${Armadillo_LIBRARIES})

#ADD_EXECUTABLE(particleViewer particleViewer.cpp ${MY_PARTICLE_WORKS_SOURCES})
#TARGET_LINK_LIBRARIES(particleViewer ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${Armadillo_LIBRARIES})

#ADD_EXECUTABLE(bsplineReg bsplineReg.cxx)
#TARGET_LINK_LIBRARIES(bsplineReg ${ITK_LIBRARIES} ${ItkVtkGlue})
