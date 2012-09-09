#ifndef __ntkPointSetViewer_h__
#define __ntkPointSetViewer_h__

#include <vtkGlyph3D.h>
#include <vtkSphereSource.h>
#include <vtkProperty.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>

#include "niralSphereMeanComputer.hpp"

namespace niral {
class SphereMeanComputerViewer {
	private:
		niral::SphereMeanComputer* _smc;
		
		vtkActor* _actor;
		vtkPolyDataMapper* _mapper;
		
		vtkActor* _meanActor;
		vtkPolyDataMapper* _meanMapper;
		
		double _radius;
		
		void deleteObject(vtkObject* obj) {
			if (obj != NULL) {
				obj->Delete();
			}
		}
		
		vtkSphereSource* createSphereSource(double radius) {
			vtkSphereSource* sphereGlyphSource = vtkSphereSource::New();
			sphereGlyphSource->SetPhiResolution(8);
			sphereGlyphSource->SetThetaResolution(8);
			sphereGlyphSource->SetRadius(radius);
			return sphereGlyphSource;			
		}
	
	public:
		SphereMeanComputerViewer() {
			_mapper = vtkPolyDataMapper::New();
			_actor = vtkActor::New();
			_meanMapper = vtkPolyDataMapper::New();
			_meanActor = vtkActor::New();			
			_smc = NULL;
			_radius = 0.01;
		};	
		
		virtual ~SphereMeanComputerViewer() {
			deleteObject(_mapper);
			deleteObject(_actor);
			deleteObject(_meanMapper);
			deleteObject(_meanActor);
		};
		
		vtkActor* GetActor() {
			return _actor;
		}
		
		vtkActor* GetMeanActor() {
			return _meanActor;
		}
		
		void IncreaseGlyphSize() {
			_radius += 0.01;
			MeanModified();
			PointsModified();
		}
		
		void DecreaseGlyphSize() {
			if (_radius > 0.01) {
				_radius -= 0.01;
			}
			MeanModified();
			PointsModified();
		}
		
		void ComputeMean() {
			_smc->ComputeMean();
			MeanModified();
		}

		void ComputeMeanInSphereCoordinate() {
			_smc->ComputeMeanInSphereCoordinate();
			MeanModified();
		}
		
		void SetData(niral::SphereMeanComputer* smc) {
			_smc = smc;
			
			MeanModified();
			PointsModified();
		}
		
		void MeanModified() {
			vtkSphereSource* sphereGlyphSource = createSphereSource(_radius);
			
			double x,y,z;
			_smc->GetMean(x,y,z);
			
			vtkTransform* tx = vtkTransform::New();
			tx->Translate(x,y,z);
			
			vtkTransformPolyDataFilter* txFilter = vtkTransformPolyDataFilter::New();
			txFilter->SetTransform(tx);
			txFilter->SetInput(sphereGlyphSource->GetOutput());
			txFilter->Update();
			
			_meanMapper->SetInput(txFilter->GetOutput());
			_meanActor->SetMapper(_meanMapper);
			
			vtkProperty* prop = _meanActor->GetProperty();
			prop->SetColor(0,1,0);
			
			sphereGlyphSource->Delete();
			tx->Delete();
			txFilter->Delete();
		}
		
		void PointsModified() {
			vtkPolyData* data = vtkPolyData::New();
			data->SetPoints(_smc->GetPoints());
			
			vtkSphereSource* sphereGlyphSource = createSphereSource(_radius);
			
			vtkGlyph3D* glyphFilter = vtkGlyph3D::New();
			glyphFilter->SetSource(sphereGlyphSource->GetOutput());
			glyphFilter->SetInput(data);

			_mapper->SetInput(glyphFilter->GetOutput());
			_actor->SetMapper(_mapper);
			
			vtkProperty* prop = _actor->GetProperty();
			prop->SetColor(1,0,0);
			
			sphereGlyphSource->Delete();
			glyphFilter->Delete();
		}
		
		void LoadData(std::string fileName) {
			vtkPolyDataReader* r = vtkPolyDataReader::New();
			r->SetFileName(fileName.c_str());
			r->Update();
			
			vtkPolyData* p = r->GetOutput();
			_smc->SetPoints(p->GetPoints());
			
			MeanModified();
			PointsModified();
		}
};
};

#endif
