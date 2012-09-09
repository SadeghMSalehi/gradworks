#ifndef __itkDeformationFieldTransformReader_h__
#define __itkDeformationFieldTransformReader_h__

#include "itkConfigure.h"
#include "itkWarpTransform3D.h"

#include "itkLightProcessObject.h"
#include "metaTransform.h"
#include "itkTransformBase.h"

namespace itk {
class DeformationFieldTransformReader : public LightProcessObject {
public:

  /** SmartPointer typedef support */
  typedef DeformationFieldTransformReader Self;
  typedef SmartPointer<Self> 							Pointer;
  typedef WarpTransform3D<double>					TransformType;
  typedef TransformType::Pointer	        TransformPointer;

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  typedef Object Superclass;
  itkTypeMacro(DeformationFieldTransformReader, LightProcessObject);

  /** Set the filename  */
  itkSetStringMacro(FileName);

  /** Get the filename */
  itkGetStringMacro(FileName);

  /** Write out the transform */
  void Update();

	TransformPointer GetTransform();

	void SetDeformationFieldToHField() {
		m_HField = true;
	}

	void SetDeformationFieldToDisplacement() {
		m_HField = false;
	}

	DeformationImageType::Pointer GetDeformationImage() {
		return m_defImg;
	}

protected:
  DeformationFieldTransformReader(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::string m_FileName;
	bool m_HField;
	TransformPointer m_txPtr;
	DeformationImageType::Pointer m_defImg;

  DeformationFieldTransformReader() {
		m_FileName = "";
		m_HField = false;
	}

	DeformationImageType::Pointer ReadDeformationFieldFile();

  virtual ~DeformationFieldTransformReader() {};

};
};

#endif
