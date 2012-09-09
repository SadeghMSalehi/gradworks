#include <vector>

#include "itkImageCommon.h"
#include "itkDeformMeshCLP.h"
#include "itkWarpTransform3D.h"
#include "itkDeformationFieldTransformReader.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkSmoothPolyDataFilter.h"

using namespace std;

typedef itk::Image<unsigned short, 3> LabelImageType;
typedef itk::DeformationFieldTransformReader TxReaderType;

bool EndsWith(const string& a, const string& b) {
    if (b.size() > a.size()) return false;
    return std::equal(a.begin() + a.size() - b.size(), a.end(), b.begin());
}

vtkPolyData* Read(string &name) {
	if (EndsWith(name, ".vtp")) {
		vtkXMLPolyDataReader* r = vtkXMLPolyDataReader::New();
		r->SetFileName(name.c_str());
		r->Update();
		return r->GetOutput();
	} else if (EndsWith(name, ".vtk")) {
		vtkPolyDataReader* r = vtkPolyDataReader::New();
		r->SetFileName(name.c_str());
		r->Update();
		return r->GetOutput();
	}
	return NULL;
}

void Write(string &name, vtkPolyData* mesh) {
	if (EndsWith(name, ".vtp")) {
		vtkXMLPolyDataWriter* w = vtkXMLPolyDataWriter::New();
		w->SetFileName(name.c_str());
		w->SetInput(mesh);
		w->SetDataModeToAppended();
		w->SetCompressorTypeToZLib();
		w->EncodeAppendedDataOff();
		w->Update();
	} else if (EndsWith(name, ".vtk")) {
		vtkPolyDataWriter* w = vtkPolyDataWriter::New();
		w->SetFileName(name.c_str());
		w->SetInput(mesh);
		w->Write();
	}
}

/**
 * converting image coordinate object to physical coordinate object
 *
 */
void ConvertIndexToPhysicalCoordinate(vtkPolyData* mesh, LabelImageType::Pointer refImg) {
	int nPts = mesh->GetNumberOfPoints();		
	vtkPoints* pts = mesh->GetPoints();	
	LabelImageType::SizeType refSize = refImg->GetBufferedRegion().GetSize();
	LabelImageType::SpacingType refSpacing = refImg->GetSpacing();
	typedef itk::ContinuousIndex<double, 3> ContinuousIndexType;

	double* inPts;
	for (int i = 0; i < nPts; i++) {
		inPts = pts->GetPoint(i);
		ContinuousIndexType inPtsIdx;

		inPtsIdx[0] = (inPts[0] / refSpacing[0]);
		inPtsIdx[1] = (inPts[1] / refSpacing[1]);
		inPtsIdx[2] = (inPts[2] / refSpacing[2]);

		//inPtsIdx[0] = refSize[0] - inPtsIdx[0] - 1;
		//inPtsIdx[1] = refSize[1] - inPtsIdx[1] - 1;

		LabelImageType::PointType physPts;
		refImg->TransformContinuousIndexToPhysicalPoint(inPtsIdx, physPts);

		inPts[0] = physPts[0];
		inPts[1] = physPts[1];
		inPts[2] = physPts[2];

		pts->SetPoint(i, inPts);
	}

	mesh->SetPoints(pts);
}

// As spharm does, convert to index space, and then multiply spacing
void ConvertPhysicalCoordinateToIndex(vtkPolyData* mesh, LabelImageType::Pointer refImg) {
	int nPts = mesh->GetNumberOfPoints();		
	vtkPoints* pts = mesh->GetPoints();	
	LabelImageType::SizeType refSize = refImg->GetBufferedRegion().GetSize();
	LabelImageType::SpacingType refSpacing = refImg->GetSpacing();
	typedef itk::ContinuousIndex<double, 3> ContinuousIndexType;

	double* inPts;
	for (int i = 0; i < nPts; i++) {
		LabelImageType::PointType inPoints;
		inPts = pts->GetPoint(i);
		inPoints[0] = inPts[0];
		inPoints[1] = inPts[1];
		inPoints[2] = inPts[2];

		ContinuousIndexType inPtsIdx;
		refImg->TransformPhysicalPointToContinuousIndex(inPoints, inPtsIdx);

		inPts[0] = inPtsIdx[0] * refSpacing[0];
		inPts[1] = inPtsIdx[1] * refSpacing[1];
		inPts[2] = inPtsIdx[2] * refSpacing[2];

		pts->SetPoint(i, inPts);
	}

	mesh->SetPoints(pts);
}



// main function
int main(int argc, char* argv[]) {
	PARSE_ARGS;

	if ("" == inputDeformation) {
		return 0;
	}

	TxReaderType::Pointer txReader = itk::DeformationFieldTransformReader::New();
	txReader->SetFileName(inputDeformation.c_str());
	if (HField) {
		cout << "convert HField to Displacement Field" << endl;
		txReader->SetDeformationFieldToHField();
	} 

	txReader->Update();
	TxReaderType::TransformPointer txWarp = txReader->GetTransform();

	if (HField) {
		if (saveDisplacementFieldName != "") {
			WriteImageT<DeformationImageType>(saveDisplacementFieldName.c_str(), txReader->GetDeformationImage());
		}
	} else {
		txWarp->SetUseTranslationOn();
	}

	string inMesh = inputMesh;
	string outMesh = outputMesh;

	vtkPolyData* mesh = Read(inMesh);

	// try to convert index-coordinate mesh (SPHARM) to physical coordinate mesh
	LabelImageType::Pointer refLabelImage = 0;
	if (SPHARM) {
		if ("" != referenceLabelImageName) {
			int ret = 0;
			refLabelImage = ReadImageT<LabelImageType>(referenceLabelImageName.c_str(), ret);	
		}

		ConvertIndexToPhysicalCoordinate(mesh, refLabelImage);
		if (processedInputName != "") {
			Write(processedInputName, mesh);
		}
	}

	cout << "reading mesh " << inMesh << endl;
	for (int i = 0; i < mesh->GetNumberOfPoints(); i++) {
		double p[3], pointDisplacement[3];
		TxReaderType::TransformType::InputPointType inPoint;
		TxReaderType::TransformType::OutputPointType outPoint;
		TxReaderType::TransformType::OutputPointType outDisplacement;

		mesh->GetPoint(i, p);
		for (int d = 0; d < 3; d++) {
				inPoint[d] = p[d];
		}

		outPoint = txWarp->TransformPoint(inPoint);
		for (int d = 0; d < 3; d++) {
				p[d] = outPoint[d];
		}

		if (printPoints) {
			cout << inPoint << ":" << outPoint << endl;
		}

		mesh->GetPoints()->SetPoint(i, p);
	} 

	if (smoothingIterations > 0) {
		vtkSmoothPolyDataFilter* smoothFilter = vtkSmoothPolyDataFilter::New();
		smoothFilter->SetInput(mesh);
		smoothFilter->SetNumberOfIterations(smoothingIterations);
		smoothFilter->Update();
		mesh = smoothFilter->GetOutput();
	}

	cout << "writing mesh " << outMesh << endl;
	Write(outMesh, mesh);

	if (outputWithIndex != "") {
		ConvertPhysicalCoordinateToIndex(mesh, refLabelImage);
		Write(outputWithIndex, mesh);
	}

}
