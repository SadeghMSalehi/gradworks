#include "itkDeformationFieldTransformReader.h"
#include "itkHFieldToDeformationFieldImageFilter.h"
#include "itkImageFileReader.h"

namespace itk {

DeformationFieldTransformReader::TransformPointer DeformationFieldTransformReader::GetTransform() {
	return m_txPtr;
}

DeformationImageType::Pointer DeformationFieldTransformReader::ReadDeformationFieldFile() {
  typedef itk::ImageFileReader<DeformationImageType> DeformationImageReader;

  DeformationImageReader::Pointer defreader = DeformationImageReader::New();
  defreader->SetFileName(m_FileName.c_str());

  if(m_HField) {
    typedef itk::HFieldToDeformationFieldImageFilter<DeformationImageType> DeformationConvertType;
    DeformationConvertType::Pointer defconv = DeformationConvertType::New();
    defconv->SetInput(defreader->GetOutput());
//  defconv->SetSpacing(timg->GetSpacing());
    defconv->Update();
    return defconv->GetOutput();

	} else {
    defreader->Update();
    return defreader->GetOutput();
	}
}

void DeformationFieldTransformReader::Update() {
	DeformationImageType::Pointer defImg = ReadDeformationFieldFile();
	DeformationFieldTransformReader::TransformPointer txPtr = DeformationFieldTransformReader::TransformType::New();
	txPtr->SetDeformationField(defImg);
	m_txPtr = txPtr;
	m_defImg = defImg;
}


};
