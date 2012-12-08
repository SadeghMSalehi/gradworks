PROJECT(ParticlesGUI)
cmake_minimum_required(VERSION 2.8)

FIND_PACKAGE(ITK 4 REQUIRED)
FIND_PACKAGE(Qt4 REQUIRED)
FIND_PACKAGE(VTK REQUIRED)

INCLUDE(${USE_ITK_FILE})
INCLUDE(${QT_USE_FILE})
INCLUDE(${USE_VTK_FILE})

SET (QT_HEADERS mainwindow.h dialogPatchCompare.h qcustomplot.h myQGraphicsView.h)
SET (QT_SOURCES mainwindow.cpp dialogPatchCompare.cpp qcustomplot.cpp myQGraphicsView.cpp)
qt4_wrap_cpp(QT_MOC_HEADERS ${QT_HEADERS})
qt4_wrap_ui(QT_UI_FILES mainwindow.ui dialogPatchCompare.ui)
add_definitions(${QT_DEFINITIONS})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

SET (MY_SOURCES myImageContainer.h myImageContainer.cpp myParticleAlgorithm.h myParticleAlgorithm.cpp vtkPropScene.h vtkPropScene.cpp PropertyAccess.h PropertyAccess.cpp myImageParticlesAlgorithm.h myImageParticlesAlgorithm.cpp myImplicitSurfaceConstraint.h myImplicitSurfaceConstraint.cpp myEnsembleEntropy.h myEnsembleEntropy.cpp vnlCommon.h vnlCommon.cpp myParticleDynamics.h myParticleDynamics.cpp myParticleSolver.h myParticleSolver.cpp myEventCallback.h myImageTransform.h myImageTransform.cpp myDeformationField.h myDeformationField.cpp myBSplineRegistration.h myBSplineRegistration.cpp)

FILE(GLOB ITK_SOURCES "itk*.h" "itk*.hxx" "itk*.cxx")

ADD_EXECUTABLE(particlesRun particlesMain.cpp ${MY_SOURCES} ${ITK_SOURCES} ${QT_HEADERS} ${QT_SOURCES} ${QT_MOC_HEADERS} ${QT_UI_FILES})
TARGET_LINK_LIBRARIES(particlesRun ${ITK_LIBRARIES} ${VTK_LIBRARIES} ${QT_LIBRARIES})
