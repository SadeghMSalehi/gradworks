PROJECT(itkapps)
cmake_minimum_required(VERSION 2.4)

if(COMMAND cmake_policy)
cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

FIND_PACKAGE(ITK 4 REQUIRED)
INCLUDE(${USE_ITK_FILE})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${USE_VTK_FILE})

set(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE( Boost COMPONENTS program_options)

FIND_PACKAGE(GenerateCLP REQUIRED)
INCLUDE(${GenerateCLP_USE_FILE})

#SET(MATLAB_ROOT /opt/matlab2010a CACHE PATH "root directory where matlab is installed")
#
#OPTION(BUILD_MATAPPS "Build Matlab related applications")
#IF (BUILD_MATAPPS)
#	INCLUDE_DIRECTORIES(${MATLAB_ROOT}/extern/include)
#	LINK_DIRECTORIES(${MATLAB_ROOT}/bin/glnxa64)
#	SET(MATLAB_LIBS mx mat)
#ENDIF (BUILD_MATAPPS)
#IF (BUILD_MATAPPS)
#	ADD_EXECUTABLE(vtk2mat vtk2mat.cxx)
#	TARGET_LINK_LIBRARIES(vtk2mat ${MATLAB_LIBS} vtkIO)
#ENDIF(BUILD_MATAPPS)
#ADD_LIBRARY(sparse STATIC /tools/sparse/lib)

ADD_EXECUTABLE(itkRandomWalk itkRandomWalker2.cxx)
TARGET_LINK_LIBRARIES(itkRandomWalk ${ITK_LIBRARIES})

ADD_EXECUTABLE(itkImageOrientation itkImageOrientation.cxx)
TARGET_LINK_LIBRARIES(itkImageOrientation ${ITK_LIBRARIES})

ADD_EXECUTABLE(itkRandomWalk2D itkRandomWalk2D.cxx)
TARGET_LINK_LIBRARIES(itkRandomWalk2D ${ITK_LIBRARIES})

ADD_EXECUTABLE(itkNewImage itkNewImage.cxx)
TARGET_LINK_LIBRARIES(itkNewImage ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkImageInverse itkImageInverse.cxx)
TARGET_LINK_LIBRARIES(itkImageInverse ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkEllipseGenerator itkEllipseGenerator.cxx)
TARGET_LINK_LIBRARIES(itkEllipseGenerator ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkGradientFilter itkGradientFilter.cxx)
TARGET_LINK_LIBRARIES(itkGradientFilter ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkHelicityFilter itkHelicityFilter.cxx)
TARGET_LINK_LIBRARIES(itkHelicityFilter ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkVectorSmoothing itkVectorSmoothing.cxx itkPolarCoordinateFilter.h) 
TARGET_LINK_LIBRARIES(itkVectorSmoothing ${ITK_LIBRARIES} vtkIO vtkCommon vtkFiltering) 

ADD_EXECUTABLE(itkSphereMeanSmoothing itkSphereMeanSmoothing.cxx niralSphereMeanComputer.hpp niralSphericalMeanCostFunction.hpp itkMeanOnVectorImageFilter.h)
TARGET_LINK_LIBRARIES(itkSphereMeanSmoothing ${ITK_LIBRARIES} vtkCommon vtkFiltering) 


#INCLUDE_DIRECTORIES(/biomed-resimg/work/joohwi/External/mathpresso/MathPresso)
#LINK_DIRECTORIES(/biomed-resimg/work/joohwi/External/mathpresso/Build)
#SET(MathPresso_LIBRARIES libmathpresso.so)
#ADD_EXECUTABLE(PixelMath PixelMath.cxx)
#TARGET_LINK_LIBRARIES(PixelMath ${ITK_LIBRARIES} ${Boost_LIBRARIES} ${MathPresso_LIBRARIES}) 

INCLUDE_DIRECTORIES(/biomed-resimg/work/joohwi/External/muparser_v2_2_0/include)
LINK_DIRECTORIES(/biomed-resimg/work/joohwi/External/muparser_v2_2_0/lib)
SET(muParser_LIBRARIES libmuparser.a)
ADD_EXECUTABLE(PixelMath PixelMath.cxx PixelMathImageFilter.h)
TARGET_LINK_LIBRARIES(PixelMath ${ITK_LIBRARIES} ${Boost_LIBRARIES} ${muParser_LIBRARIES}) 

ADD_EXECUTABLE(itkImageInfo itkImageInfo.cxx)
TARGET_LINK_LIBRARIES(itkImageInfo ${ITK_LIBRARIES} ${Boost_LIBRARIES}) 

ADD_EXECUTABLE(itkRigidTransform itkRigidTransform.cxx)
TARGET_LINK_LIBRARIES(itkRigidTransform ${ITK_LIBRARIES} ${Boost_LIBRARIES}) 

ADD_EXECUTABLE(itkImageCompress itkImageCompress.cxx)
TARGET_LINK_LIBRARIES(itkImageCompress ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkChangeLabel itkChangeLabel.cxx)
TARGET_LINK_LIBRARIES(itkChangeLabel ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkPCNNPostProcessing itkPCNNPostProcessing.cxx)
TARGET_LINK_LIBRARIES(itkPCNNPostProcessing ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkLabelMapCompress itkLabelMapCompress.cxx)
TARGET_LINK_LIBRARIES(itkLabelMapCompress ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkLabelGeometry itkLabelGeometry.cxx)
TARGET_LINK_LIBRARIES(itkLabelGeometry ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkKmeansFilter itkKmeansFilter.cxx)
TARGET_LINK_LIBRARIES(itkKmeansFilter ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkEdgePreprocessing itkEdgePreprocessing.cxx EdgePreprocessingSettings.cxx)
TARGET_LINK_LIBRARIES(itkEdgePreprocessing ${ITK_LIBRARIES})

SET(dtiCaseMeasure_SRC dtiCaseMeasure.cxx)
GENERATECLP(dtiCaseMeasure_SRC dtiCaseMeasure.xml)
ADD_EXECUTABLE(dtiCaseMeasure ${dtiCaseMeasure_SRC})
TARGET_LINK_LIBRARIES(dtiCaseMeasure ${ITK_LIBRARIES} ${Boost_LIBRARIES}) 

ADD_EXECUTABLE(itkTextToImage itkTextToImage.cxx itkImageCommon.h)
TARGET_LINK_LIBRARIES(itkTextToImage ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkImageToVTI itkImageToVTI.cxx itkImageCommon.h)
TARGET_LINK_LIBRARIES(itkImageToVTI ${ITK_LIBRARIES} vtkCommon vtkFiltering vtkIO ${Boost_LIBRARIES}) 

#ADD_EXECUTABLE(vtkImageToITK vtkImageToITK.cxx itkImageCommon.h)
#TARGET_LINK_LIBRARIES(vtkImageToITK ${ITK_LIBRARIES} vtkCommon vtkFiltering vtkIO ${Boost_LIBRARIES}) 

#ADD_EXECUTABLE(itkImageToVTKVolume itkImageToVTKVolume.cxx itkImageCommon.h)
#TARGET_LINK_LIBRARIES(itkImageToVTKVolume ${ITK_LIBRARIES} vtkCommon vtkFiltering vtkIO) 

ADD_EXECUTABLE(itkSplitCortexByMoments itkSplitCortexByMoments.cxx itkImageCommon.h)
TARGET_LINK_LIBRARIES(itkSplitCortexByMoments ${ITK_LIBRARIES} vtkCommon vtkFiltering vtkIO) 

ADD_EXECUTABLE(itkMomentumComputer itkMomentumComputer.cxx itkImageCommon.h)
TARGET_LINK_LIBRARIES(itkMomentumComputer ${ITK_LIBRARIES} vtkCommon vtkFiltering vtkIO) 

ADD_EXECUTABLE(itkSyncHeader itkSyncHeader.cxx itkImageCommon.h)
TARGET_LINK_LIBRARIES(itkSyncHeader ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkAbsFilter itkAbsFilter.cxx itkImageCommon.h)
TARGET_LINK_LIBRARIES(itkAbsFilter ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkJacobianFilter itkJacobianFilter.cxx itkImageCommon.h)
TARGET_LINK_LIBRARIES(itkJacobianFilter ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkCurlFilter itkCurlFilter.cxx itkImageCommon.h)
TARGET_LINK_LIBRARIES(itkCurlFilter ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkVectorComponent itkVectorComponent.cxx itkImageCommon.h)
TARGET_LINK_LIBRARIES(itkVectorComponent ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkSphericalComponent itkSphericalComponent.cxx itkImageCommon.h)
TARGET_LINK_LIBRARIES(itkSphericalComponent ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkRotationSample itkRotationSample.cxx itkImageCommon.h)
TARGET_LINK_LIBRARIES(itkRotationSample ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(rgseg MainSegmentationByRegionGrowing.cxx itkImageCommon.h itkRegionGrowingAlgorithm.h)
TARGET_LINK_LIBRARIES(rgseg ${ITK_LIBRARIES} ${Boost_LIBRARIES}) 

ADD_EXECUTABLE(inp inp.cxx itkImageCommon.h itkRegionGrowingAlgorithm.h)
TARGET_LINK_LIBRARIES(inp ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(dti_double2float dti_double2float.cxx itkImageCommon.h)
TARGET_LINK_LIBRARIES(dti_double2float ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkGrayDilate itkGrayDilate.cxx)
TARGET_LINK_LIBRARIES(itkGrayDilate ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkGradientMagnitude itkGradientMagnitude.cxx)
TARGET_LINK_LIBRARIES(itkGradientMagnitude ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkIntensityThresholdFilter itkThreshFilter.cxx)
TARGET_LINK_LIBRARIES(itkIntensityThresholdFilter ${ITK_LIBRARIES}) 

ADD_EXECUTABLE(itkRescaleFilter itkRescaleFilter.cxx)
TARGET_LINK_LIBRARIES(itkRescaleFilter ${ITK_LIBRARIES})

ADD_EXECUTABLE(itkAdaptiveHistogramEqualizer itkAdaptiveHistogramEqualizer.cxx)
TARGET_LINK_LIBRARIES(itkAdaptiveHistogramEqualizer ${ITK_LIBRARIES})

ADD_EXECUTABLE(itkVectorProbe itkVectorProbe.cxx)
TARGET_LINK_LIBRARIES(itkVectorProbe ${ITK_LIBRARIES} ${Boost_LIBRARIES} vtkCommon vtkFiltering vtkIO)

ADD_EXECUTABLE(itkVectorSeg itkVectorSeg.cxx)
TARGET_LINK_LIBRARIES(itkVectorSeg ${ITK_LIBRARIES})

ADD_EXECUTABLE(itkReplaceLabel itkReplaceLabel.cxx)
TARGET_LINK_LIBRARIES(itkReplaceLabel ${ITK_LIBRARIES})

#ADD_EXECUTABLE(itkConComp itkConComp.cxx)
#TARGET_LINK_LIBRARIES(itkConComp ${ITK_LIBRARIES})

SET(itkDeformMesh_SRC itkDeformMesh.cxx itkDeformationFieldTransformReader.cxx)
GENERATECLP(itkDeformMesh_SRC itkDeformMesh.xml)
ADD_EXECUTABLE(itkDeformMesh ${itkDeformMesh_SRC})
TARGET_LINK_LIBRARIES(itkDeformMesh ${ITK_LIBRARIES} vtkCommon vtkFiltering vtkIO vtkGraphics)

SET(PCNNSegmentation_SRC PCNNSegmentation.cxx)
GENERATECLP(PCNNSegmentation_SRC PCNNSegmentation.xml)
ADD_EXECUTABLE(PCNNSegmentation ${PCNNSegmentation_SRC})
TARGET_LINK_LIBRARIES(PCNNSegmentation ${ITK_LIBRARIES})

#ADD_EXECUTABLE(itkVectorSegHelper itkVectorSegHelper.cxx)
#TARGET_LINK_LIBRARIES(itkVectorSegHelper ${ITK_LIBRARIES})

ADD_EXECUTABLE(eigen2vti eigen2vti.cxx)
TARGET_LINK_LIBRARIES(eigen2vti ${ITK_LIBRARIES} vtkCommon vtkFiltering vtkIO)

ADD_EXECUTABLE(vtkVectorVis vtkVectorVis.cxx)
TARGET_LINK_LIBRARIES(vtkVectorVis ${ITK_LIBRARIES} vtkCommon vtkRendering vtkGraphics vtkFiltering)

ADD_EXECUTABLE(vtkTensorVis vtkTensorVis.cxx)
TARGET_LINK_LIBRARIES(vtkTensorVis ${ITK_LIBRARIES} vtkCommon vtkRendering vtkGraphics vtkFiltering)

ADD_EXECUTABLE(vtkProjectMesh vtkProjectMesh.cxx)
TARGET_LINK_LIBRARIES(vtkProjectMesh vtkCommon vtkRendering vtkGraphics vtkFiltering)

ADD_EXECUTABLE(vtkSurfaceDistance vtkSurfaceDistance.cxx)
TARGET_LINK_LIBRARIES(vtkSurfaceDistance vtkCommon vtkRendering vtkGraphics vtkFiltering)

SET(itkLabelMapMorphology_SRC itkLabelMapMorphology.cxx)
GENERATECLP(itkLabelMapMorphology_SRC itkLabelMapMorphology.xml)
ADD_EXECUTABLE(itkLabelMapMorphology ${itkLabelMapMorphology_SRC})
TARGET_LINK_LIBRARIES(itkLabelMapMorphology ${ITK_LIBRARIES})


SET(itkSampleMesh_SRC itkSampleMesh.cxx vtkPolyDataToitkMesh.cxx vtkPolyDataConnectivityFilter2.cxx)
GENERATECLP(itkSampleMesh_SRC itkSampleMesh.xml)
ADD_EXECUTABLE(itkSampleMesh ${itkSampleMesh_SRC})
TARGET_LINK_LIBRARIES(itkSampleMesh ${ITK_LIBRARIES} vtkIO vtkCommon vtkFiltering)

ADD_EXECUTABLE(itkcmds itkcmds.cxx itkImageIO.h)
TARGET_LINK_LIBRARIES(itkcmds ${ITK_LIBRARIES})

ADD_EXECUTABLE(itkImageMatch itkImageMatch.cxx itkImageIO.h)
TARGET_LINK_LIBRARIES(itkImageMatch ${ITK_LIBRARIES})
