#include "vtkSphereSource.h"
#include "vtkLineSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkActor.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkUnstructuredGrid.h"
#include "vtkGlyph3D.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkArrowSource.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkCallbackCommand.h"
#include "vtkRendererCollection.h"
#include "vtkLookupTable.h"
#include "vtkColorTransferFunction.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGridWriter.h"

#include "itkVector.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include <vector>
#include <iostream>

#define VTK_CREATE(type, name) \
	type *name = type::New()

//  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

using namespace std;

typedef itk::Vector<double,3> VectorType;
typedef itk::Image<VectorType,3> VectorImageType;
typedef itk::Image<unsigned short,3> UShortImageType;
typedef itk::ImageFileReader<VectorImageType> VectorImageReaderType;
typedef itk::ImageFileReader<UShortImageType> UShortImageReaderType;
typedef itk::ImageRegionConstIteratorWithIndex<UShortImageType> UShortIteratorType;
typedef itk::ImageRegionConstIteratorWithIndex<VectorImageType> VectorIteratorType;
// typedef itk::Vector<UShortImageType::IndexType, 1> IndexVectorType;
typedef std::vector<UShortImageType::IndexType> IndexVectorType;
typedef std::vector<UShortImageType::IndexType>::iterator IndexVectorIteratorType;

int viewSliceNo = 0;
int maxSeed = 0;
IndexVectorType indexVector;
IndexVectorType seedVector;

VectorImageType::Pointer vectImage = NULL;
UShortImageType::Pointer roiImage = NULL;
UShortImageType::Pointer seedImage = NULL;

char roiImageName[256];
char seedImageName[256];

static void constructArrows(vtkRenderer* ren) {
  VTK_CREATE(vtkPoints, points);
  VTK_CREATE(vtkUnstructuredGrid, grid);
  VTK_CREATE(vtkDoubleArray, pointData);
  VTK_CREATE(vtkDoubleArray, colorData);
  VTK_CREATE(vtkLookupTable, lut);
  VTK_CREATE(vtkColorTransferFunction, ctf);


  pointData->SetNumberOfComponents(3);
  colorData->SetNumberOfComponents(1);

  IndexVectorIteratorType iter;
  IndexVectorIteratorType seedIter;
  for (iter = indexVector.begin(); iter != indexVector.end(); iter++) {
    VectorImageType::IndexType roiIdx = *iter;
    roiIdx[1] = viewSliceNo;

    VectorType eigenVector = vectImage->GetPixel(roiIdx);

    if (eigenVector[0] == 0 && eigenVector[1] == 0 && eigenVector[2] == 0) {
      continue;
    }	

    bool foundMatchingSeed = false;
    for (seedIter = seedVector.begin(); seedIter != seedVector.end(); seedIter++) {
      UShortImageType::IndexType seedIdx = *seedIter;
      if (seedIdx[0] == roiIdx[0] && seedIdx[1] == roiIdx[1] && seedIdx[2] == roiIdx[2]) {
        foundMatchingSeed = true;
        int val = seedImage->GetPixel(seedIdx);
        double orderValue = (double(val) / maxSeed * 0.9 + 0.1);
        colorData->InsertNextTuple1(orderValue);
        cout << "found seed " << roiIdx << endl;
        break;
      }
    }

    points->InsertNextPoint(roiIdx[0], roiIdx[1], roiIdx[2]);
    pointData->InsertNextTuple3(eigenVector[0], eigenVector[1], eigenVector[2]);

    if (!foundMatchingSeed) {
      colorData->InsertNextTuple1(0);
    }
  }
  
  grid->SetPoints(points);
  grid->GetPointData()->SetVectors(pointData);
  grid->GetPointData()->SetScalars(colorData);

  VTK_CREATE(vtkUnstructuredGridWriter, gridWriter);
  char filename[256];
  snprintf(filename, 256, "vectorData.%d.vtk", viewSliceNo);

  gridWriter->SetFileName(filename);
  gridWriter->SetInput(grid);
  gridWriter->Write();

  VTK_CREATE(vtkArrowSource, arrowSource);
  VTK_CREATE(vtkGlyph3D, glyph);

  glyph->SetInput(grid);
  glyph->SetSource(arrowSource->GetOutput());
  glyph->SetScaleModeToScaleByVector();
  glyph->SetColorModeToColorByScalar();
  glyph->SetScaleFactor(1);
  glyph->Update();

  VTK_CREATE(vtkPolyDataMapper, glyphMapper);
  glyphMapper->SetInputConnection(glyph->GetOutputPort());

  int lutNum = 256;
  ctf->SetColorSpaceToHSV();
  ctf->HSVWrapOn();
  ctf->AddRGBPoint(0.0, 0.3, 0.3, 0.3);
  ctf->AddRGBPoint(0.00000000001, 0.0, 1.0, 0.0);
  ctf->AddRGBPoint(1.0, 1.0, 0, 0);
  lut->SetNumberOfTableValues(lutNum);

  for (int x = 0; x < 256; x++) {
    double* color = ctf->GetColor(double(x) / double(lutNum));
    lut->SetTableValue(x, color[0], color[1], color[2], 1.0);
  }

  glyphMapper->SetLookupTable(lut);
  glyphMapper->SetScalarRange(0, 1);
  glyphMapper->SetScalarModeToUsePointData();

  VTK_CREATE(vtkActor, glyphActor);
  glyphActor->SetMapper(glyphMapper);

  ren->RemoveAllViewProps();
  ren->AddActor(glyphActor);
  ren->Modified();
}

static void loadROIandSeedImage() {
  cout << "Load Images " << roiImageName << " / " << seedImageName << endl;

	// if roi is provided
  UShortImageReaderType::Pointer roiReader = UShortImageReaderType::New();
  roiReader->SetFileName(roiImageName);
  roiReader->Update();
  roiImage = roiReader->GetOutput();

	// if seed is provided
  UShortImageReaderType::Pointer seedReader = UShortImageReaderType::New();
  seedReader->SetFileName(seedImageName);
  seedReader->Update();
  seedImage = seedReader->GetOutput();

  UShortIteratorType it(roiImage, roiImage->GetRequestedRegion());
  it.GoToBegin();

  indexVector.clear();
  for (it = it.Begin(); !it.IsAtEnd(); ++it) {
    UShortImageType::IndexType idx = it.GetIndex();
    if (it.Get() != 0) {
      // temporary assume axial slice is 2nd element of the index
      if (viewSliceNo == 0) {
        viewSliceNo = idx[1];
      }
      indexVector.push_back(idx);
    }
  }

  UShortIteratorType it2(seedImage, seedImage->GetRequestedRegion());
  it2.GoToBegin();

  seedVector.clear();
  for (it2 = it2.Begin(); !it2.IsAtEnd(); ++it2) {
    UShortImageType::IndexType idx = it2.GetIndex();
    int val = it2.Get();
    if (val > maxSeed) {
      maxSeed = val;
    }
    if (it2.Get() != 0) {
        seedVector.push_back(idx);
    }
  }
}


static void processKeyPressEvent(vtkObject* caller, unsigned long event, void* clientData, void* callData) {
  IndexVectorIteratorType iter;
  vtkRenderer* renderer = NULL;
	switch (event) {
		case vtkCommand::CharEvent: 
			vtkRenderWindowInteractor* interactor = vtkRenderWindowInteractor::SafeDownCast(caller);
			char keyCode = interactor->GetKeyCode();

			switch (keyCode) {
				case '8':
					for (iter = indexVector.begin(); iter != indexVector.end(); iter ++) {
						cout << *iter << endl;
					}
          renderer = interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
          renderer->ResetCamera();
          interactor->Modified();
          interactor->Render();
          break;
				case '9':
          if (viewSliceNo > 0) {
            viewSliceNo --;
          }
          renderer = interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
          constructArrows(renderer);
          interactor->Modified();
          interactor->Render();
					break;
				case '0':
          // need bound check!
          if (viewSliceNo < 300) {
            viewSliceNo ++;
          }
          renderer = interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
          constructArrows(renderer);
          interactor->Modified();
          interactor->Render();
					break;
        case 'l':
          renderer = interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
          constructArrows(renderer);
          loadROIandSeedImage();
          constructArrows(renderer);
          interactor->Modified();
          interactor->Render();
          break;
			}

			cout << "current slice: " << viewSliceNo << endl;
			break;
	}
}

int main (int argc, char* argv[]) {
  if (argc < 4) {
    cout << argv[0] << " in-vector-image [roi-image] [segmentation-image]" << endl;
    return 0;
  } 
  VectorImageReaderType::Pointer vectReader = VectorImageReaderType::New();
  vectReader->SetFileName(argv[1]);
  vectReader->Update();
  vectImage = vectReader->GetOutput();

  // a renderer and render window
  vtkRenderer *ren1 = vtkRenderer::New();
  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren1);
  ren1->SetBackground(1,1,1); // Background color white

  // an interactor
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);

	vtkCallbackCommand* keyPressEvent = vtkCallbackCommand::New();
	keyPressEvent->SetCallback(processKeyPressEvent);
	keyPressEvent->SetClientData(NULL);
	iren->AddObserver(vtkCommand::CharEvent, keyPressEvent, 1.0);

  sprintf(roiImageName, "%s", argv[2]);
  sprintf(seedImageName, "%s", argv[3]);

  loadROIandSeedImage();

  constructArrows(ren1);

  // render an image (lights and cameras are created automatically)
  renWin->Render();

  // begin mouse interaction
  iren->Start();

  // release memory and return
  ren1->Delete();
  renWin->Delete();
  iren->Delete();

  return EXIT_SUCCESS;
} 
