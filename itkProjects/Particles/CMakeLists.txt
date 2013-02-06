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

SET (QT_PW_HEADERS particleViewerWindow.h myQGraphicsView.h)
SET (QT_PW_SOURCES particleViewerWindow.cpp  myQGraphicsView.cpp)
qt4_wrap_cpp(QT_PW_MOC_HEADERS ${QT_PW_HEADERS})
qt4_wrap_ui(QT_PW_UI_FILES particleViewerWindow.ui)

add_definitions(${QT_DEFINITIONS})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

SET (MY_Qt_EXTENSIONS QGraphicsGridItem.cpp QGraphicsGridItem.h)

SET (MY_PARTICLE_WORKS_SOURCES myParticleCore.cpp myParticleCore.h 
	myParticleConstraint.cpp myParticleConstraint.h 
	myImageDef.h 
	myParticleBSpline.cpp myParticleBSpline.h 
	myParticleForces.h myParticleForces.cpp)

FILE(GLOB ITK_SOURCES "itk*.h" "itk*.hxx" "itk*.cxx")

ADD_EXECUTABLE(particleWin particleViewerWindowMain.cpp ${QT_PW_HEADERS} ${QT_PW_SOURCES} ${QT_PW_MOC_HEADERS} ${QT_PW_UI_FILES} ${MY_PARTICLE_WORKS_SOURCES} ${MY_Qt_EXTENSIONS})
TARGET_LINK_LIBRARIES(particleWin ${ITK_LIBRARIES} ${QT_LIBRARIES} ${Armadillo_LIBRARIES})

ADD_EXECUTABLE(particleRun particlesRun.cpp ${MY_PARTICLE_WORKS_SOURCES})
TARGET_LINK_LIBRARIES(particleRun ${ITK_LIBRARIES} ${Armadillo_LIBRARIES})

ADD_EXECUTABLE(particleReg particleRegistrationMain.cpp ${MY_PARTICLE_WORKS_SOURCES})
TARGET_LINK_LIBRARIES(particleReg ${ITK_LIBRARIES} ${Armadillo_LIBRARIES})
