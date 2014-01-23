PROJECT(ktools)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE(ITK 4 REQUIRED)
INCLUDE(${ITK_USE_FILE})

SET(USE_VTK False CACHE BOOL "Use VTK Library and related applications")

SET(BUILD_KCALC False CACHE BOOL "Compile kcalc - a pixel-wise calculation tool")
SET(BUILD_KMBRAIN False CACHE BOOL "Compile kmbrain - a rodent brain localization tool")
SET(BUILD_KTOOLS False CACHE BOOL "Compile ktools - ksegeval, krgb, korient ...")

INCLUDE_DIRECTORIES(muparser)
FILE(GLOB muparser_SRC "muparser/*.cpp")
FILE(GLOB pi_SRC "pi*")

IF(${BUILD_KTOOLS})
ADD_EXECUTABLE(korient korient.cpp piOptions.cpp piOptions.h)
SET_TARGET_PROPERTIES(korient PROPERTIES COMPILE_FLAGS "-DPIXEL_TYPE=float")
TARGET_LINK_LIBRARIES(korient ${ITK_LIBRARIES})

ADD_EXECUTABLE(ksegeval ksegeval.cpp piImageIO.h piOptions.cpp piOptions.h piImageProcessing.cpp piImageProcessing.h)
TARGET_LINK_LIBRARIES(ksegeval ${ITK_LIBRARIES})

ADD_EXECUTABLE(krgb krgb.cpp piImageIO.h piOptions.cpp piOptions.h piImageProcessing.cpp piImageProcessing.h)
TARGET_LINK_LIBRARIES(krgb ${ITK_LIBRARIES})

ADD_EXECUTABLE(kresample kresample.cpp piImageIO.h)
TARGET_LINK_LIBRARIES(kresample ${ITK_LIBRARIES})

ADD_EXECUTABLE(kextract kextract.cpp piImageIO.h piOptions.h piOptions.cpp)
TARGET_LINK_LIBRARIES(kextract ${ITK_LIBRARIES})


ENDIF(${BUILD_KTOOLS})



IF(${BUILD_KCALC})
ADD_EXECUTABLE(kcalc kcalc.cpp piImageIO.h piOptions.cpp piOptions.h ${muparser_SRC} PixelMathImageFilter.h PixelMathImageFilter.hxx)
SET_TARGET_PROPERTIES(kcalc PROPERTIES COMPILE_FLAGS "-DPIXEL_TYPE=float")
TARGET_LINK_LIBRARIES(kcalc ${ITK_LIBRARIES})
ADD_CUSTOM_COMMAND(TARGET kcalc POST_BUILD COMMAND kcalc -h > ${CMAKE_CURRENT_SOURCE_DIR}/../kcalc.md)
ENDIF(${BUILD_KCALC})

IF(${USE_VTK})
	FIND_PACKAGE(VTK REQUIRED)
	INCLUDE(${VTK_USE_FILE})
	SET(BUILD_KVTKTOOLS False CACHE BOOL "Compile ktools with VTK - kgradmap, ksulcus ...")

	IF(${BUILD_KVTKTOOLS})
		ADD_EXECUTABLE(kgradmap kgradmap.cpp piImageIO.h piOptions.cpp piOptions.h piVTK.h piVTK.cpp piImageProcessing.cpp piImageProcessing.h piTimer.h piTimer.cpp)
		TARGET_LINK_LIBRARIES(kgradmap ${ITK_LIBRARIES} ${VTK_LIBRARIES})
	ENDIF(${BUILD_KVTKTOOLS})
    
    
    ADD_EXECUTABLE(ksulcus ksulcus.cpp piImageIO.h piOptions.h piOptions.cpp)
    TARGET_LINK_LIBRARIES(ksulcus ${ITK_LIBRARIES} ${VTK_LIBRARIES})
    
    ADD_EXECUTABLE(kvoxelizer 
        kvoxelizer.h kvoxelizer.cpp
        piOptions.h piOptions.cpp
        vtkio.h vtkio.cpp
    )
    TARGET_LINK_LIBRARIES(kvoxelizer ${ITK_LIBRARIES} ${VTK_LIBRARIES})
    
    ADD_EXECUTABLE(kmesh
        kmesh.h kmesh.cpp
        kimage.h kimage.cpp
        piOptions.h piOptions.cpp
        vtkio.h vtkio.cpp
    )
    TARGET_LINK_LIBRARIES(kmesh ${ITK_LIBRARIES} ${VTK_LIBRARIES})
    ADD_CUSTOM_COMMAND(TARGET kmesh POST_BUILD COMMAND kmesh -h > ${CMAKE_CURRENT_SOURCE_DIR}/kmesh.md)
ENDIF(${USE_VTK})

IF(${BUILD_KMBRAIN})
ADD_EXECUTABLE(kmbrain kmbrain.cpp piImageIO.h piOptions.cpp piOptions.h piImageProcessing.cpp piImageProcessing.h piTimer.h piTimer.cpp)
TARGET_LINK_LIBRARIES(kmbrain ${ITK_LIBRARIES})
ENDIF(${BUILD_KMBRAIN})

ADD_CUSTOM_TARGET(ksync rsync -av --delete --delete-excluded -r /tools/gradworks/itkProjects/ktools /NIRAL/work/joohwi/m2013/ktools --exclude=".*" --no-g --no-o)
