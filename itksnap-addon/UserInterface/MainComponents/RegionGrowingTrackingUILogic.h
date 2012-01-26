/*=========================================================================

  Program:   ITK-SNAP
  Module:    $RCSfile: RegionGrowingTrackingUILogic.h,v $
  Language:  C++
  Date:      $Date: 2009/10/30 16:48:24 $
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

#ifndef __RegionGrowingTrackingUILogic_h_
#define __RegionGrowingTrackingUILogic_h_

#include "RegionGrowingTrackingUI.h"
#include <vnl/vnl_matrix_fixed.h>
#include "itkVector.h"
#include "itkOrientedImage.h"
#include "iostream"



class UserInterfaceLogic;
class GenericImageData;


class RegionGrowingTrackingUILogic : public RegionGrowingTrackingUI
{
public:
  // Constructor
  RegionGrowingTrackingUILogic();
  virtual ~RegionGrowingTrackingUILogic() {}
  
  // Callbacks
  void OnOkAction();
  void OnApplyAction();
  void OnCloseAction();

  void ShowDialog();
  void Register(UserInterfaceLogic *parent_ui);

  bool Shown();

private:
		
  UserInterfaceLogic *m_ParentUI;
  GenericImageData *m_ImageData;

};


#endif

