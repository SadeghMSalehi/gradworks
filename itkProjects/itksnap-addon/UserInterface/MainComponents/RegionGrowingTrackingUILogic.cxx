/*=========================================================================

  Program:   ITK-SNAP
  Module:    $RCSfile: RegionGrowingTrackingUILogic.cxx,v $
  Language:  C++
  Date:      $Date: 2009/10/30 16:48:22 $
  Version:   $Revision: 1.0 $
  Copyright (c) 2011 Joohwi Lee
  
  This file is part of ITK-SNAP 

  ITK-SNAP is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  -----

  Copyright (c) 2003 Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information. 

=========================================================================*/

#include "SNAPCommonUI.h"
#include "RegionGrowingTrackingUILogic.h"
#include "UserInterfaceLogic.h"
#include "IRISApplication.h"
#include "GenericImageData.h"
#include "itkOrientedImage.h"

#include <cassert>
#include <string>
#include <iomanip>
#include <sstream>

#include "itkImageCommon.h"
#include "itkRegionGrowingAlgorithm.h"

#include "DCRegionGrowingAlgorithm.h"

#define mag3(x) (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])

typedef itk::Vector<double,3> VectorType;
typedef itk::OrientedImage<VectorType,3> VectorImageType;

std::string m_dtiFileName;
VectorImageType::Pointer m_dtiImage;

std::string m_f1FileName;
niral::FAImageType::Pointer m_f1Image;
bool m_f1Loaded = false;

std::string m_f2FileName;
niral::FAImageType::Pointer m_f2Image;
bool m_f2Loaded = false;
	
bool LoadFeatureImage(const char* inputFileName, std::string &loadedFileName, niral::FAImageType::Pointer &featureImage) {
	if (loadedFileName != std::string(inputFileName)) {
		int ret = 0;
		featureImage = ReadImageT<niral::FAImageType>(inputFileName, ret);
		loadedFileName = std::string(inputFileName);
		return true;
	} else if (strlen(inputFileName) == 0) {
		return false;
	}
	return true;
}

RegionGrowingTrackingUILogic
::RegionGrowingTrackingUILogic()
{

}

void
RegionGrowingTrackingUILogic
::ShowDialog()
{
  assert(m_ParentUI);

  // Get the current image data  
  m_ImageData = m_ParentUI->GetDriver()->GetCurrentImageData();
	m_inputDTI->value(m_dtiFileName.c_str());
	
	Vector3ui cursorPos = m_ParentUI->GetDriver()->GetCursorPosition();
	static char str[255]; 
	sprintf(str, "%d %d %d", cursorPos[0], cursorPos[1], cursorPos[2]);
	
	Fl_Text_Buffer* buf = m_dtiVectorInfo->buffer();
	buf->text(str);
	
  // Show the window
  m_WinRGTrac->show();
  m_ParentUI->CenterChildWindowInMainWindow(m_WinRGTrac);
}

void
RegionGrowingTrackingUILogic
::Register(UserInterfaceLogic *parent_ui)
{ 
  m_ParentUI = parent_ui; 
}


void 
RegionGrowingTrackingUILogic
::OnOkAction()
{
  this->OnApplyAction();
  this->OnCloseAction();
}

void 
RegionGrowingTrackingUILogic
::OnApplyAction()
{
	// region growing parameters
	if (!m_ImageData->IsSegmentationLoaded()) {
		printf("No segmentation is loaded.\n");
		return;
	}

	if (!m_ImageData->IsGreyLoaded()) {
		printf("No grey image is loaded.\n");
		return;
	}
		
	const char* inputDTI = m_inputDTI->value();
	if (m_dtiFileName != std::string(inputDTI)) {
		int ret = 0;
		m_dtiImage = ReadImageT<VectorImageType>(inputDTI, ret);
		if (ret != 0) {
			m_dtiImage = NULL;
			printf("Failed to load '%s'\n", inputDTI);
			return;
		}
		m_dtiFileName = std::string(inputDTI);
	}

	m_f1Loaded = LoadFeatureImage(m_inputF1->value(), m_f1FileName, m_f1Image);
	m_f2Loaded = LoadFeatureImage(m_inputF2->value(), m_f2FileName, m_f2Image);
	
	LabelImageWrapper* labelsWrapper = m_ImageData->GetSegmentation();
	LabelImageWrapper::ImagePointer seedImage = labelsWrapper->GetImage();
	
	VectorImageType::SizeType szDTI = m_dtiImage->GetRequestedRegion().GetSize();
	LabelImageWrapper::ImageType::SizeType szLabel = seedImage->GetLargestPossibleRegion().GetSize();
	
	if (szDTI[0] != szLabel[0] || szDTI[1] != szLabel[1] || szDTI[2] != szLabel[2]) {
		printf("Image sizes do not match\n");
		m_dtiImage = NULL;
		return;
	}

	niral::DCRegionGrowingAlgorithm<LabelImageWrapper::ImageType> dcrga(m_dtiImage, seedImage);
	dcrga.SetInitialStd(m_initialStd->value());
	dcrga.SetCutoffStd(m_cutoffStd->value());
	dcrga.SetFeature1Values(m_f1MinThreshold->value(), m_f1MaxThreshold->value(), m_f1Delta->value());
	dcrga.SetFeature2Values(m_f2MinThreshold->value(), m_f2MaxThreshold->value(), m_f2Delta->value());
	dcrga.SetMaxNumberOfPixels(m_maxNumberOfVoxels->value());	
	if (m_f1Loaded) {
		dcrga.SetFeature1Image(m_f1Image);
	}
	if (m_f2Loaded) {
		dcrga.SetFeature2Image(m_f2Image);
	}
	dcrga.Init();
	
  int outputLabel = 1;
  for (double angleInput = m_angleMinThreshold->value(); angleInput < m_angleMaxThreshold->value(); angleInput += 0.1) {
    double dc = cos(angleInput * M_PI / 180.0);
    cout << "Trying DCt = " << dc << " (angle = " << angleInput << " degree; outputLabel = " << outputLabel << ")" << endl;
    dcrga.SetDC(dc);
    if (m_useIncrementalLabel->value()) {
			dcrga.SetOutputLabel(outputLabel++);
		}
    dcrga.DoRegionGrowing();
		cout << "Number of propagated pixels = " << dcrga.GetVoxelCount() << endl;

    if (dcrga.IsOverPixels() || !dcrga.Reinitialize()) {
      break;
    }   
  }
	m_ImageData->SetSegmentationImage(dcrga.GetOutput());
	
	LabelImageWrapper::ImageType::SizeType szProp = dcrga.GetPropagationImage()->GetBufferedRegion().GetSize();
	cout << "Propagation Image Size: " << szProp[0] << ", " << szProp[1] << ", " << szProp[2] << endl;
	
	WriteImageT<LabelImageWrapper::ImageType>("propagation.gipl.gz", dcrga.GetPropagationImage());
	
	//m_ParentUI->OnImageGeometryUpdate();
	m_ParentUI->RedrawWindows();
}

void 
RegionGrowingTrackingUILogic
::OnCloseAction()
{
  m_WinRGTrac->hide();
}

bool
RegionGrowingTrackingUILogic
::Shown()
{
  return m_WinRGTrac->shown();
}

