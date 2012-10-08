/*
 * itkMultiRegMethod.cpp
 *
 *  Created on: Sep 24, 2012
 *      Author: joohwi
 */

#include "itkMultiRegMethod.h"
#include "itkContinuousIndex.h"
#include "itkMyTransformFilter.h"
#include "sstream"

void itkMultiRegMethod::Execute(itk::Object* caller,
		const itk::EventObject& event) {
	Execute((const Object*) caller, event);
}

void itkMultiRegMethod::OnIterate(const OptimizerType* optimizer) {
	if (++_numOfIterations % _updateInterval == 0) {
		itk::RealTimeClock::TimeStampType t = _clock->GetTimeInSeconds();
		std::cout << _numOfIterations << "\t" << optimizer->GetValue() << "\t"
				<< optimizer->GetCurrentPosition() << "\t" << (t - _lastTime)
				<< " secs" << std::endl;
		_transformHistory.push_back(optimizer->GetCurrentPosition());

		_lastTime = t;
	}
}

void itkMultiRegMethod::Execute(const itk::Object* object,
		const itk::EventObject& event) {
	if (object == NULL) {
		std::cout << "Null sender is not processed..." << std::endl;
		return;
	}
	if (typeid(event) == typeid(itk::IterationEvent)) {
		const OptimizerType* optimizer =
				dynamic_cast<const OptimizerType*>(object);
		if (optimizer == NULL) {
			cout << "Wrong optimizer type" << endl;
			return;
		}
		OnIterate(optimizer);
	} else if (typeid(event) == typeid(itk::StartEvent)) {
		std::cout << "Optimization has started ..." << std::endl;
		_lastTime = _clock->GetTimeInSeconds();
	} else if (typeid(event) == typeid(itk::EndEvent)) {
		std::cout << "Optimization has finisehd ..." << std::endl;
	}
}

void itkMultiRegMethod::Update() {
	this->Execute((const Object*) NULL, itk::IterationEvent());
}

void itkMultiRegMethod::SetImages(ImageType::Pointer source,
		LabelType::Pointer sourceLabel, ImageType::Pointer target) {
	_movingImage = target;
	_fixedImage = source;
	_fixedImageLabel = sourceLabel;
}


void itkMultiRegMethod::ComputeLabelIndexes() {
    QVector<int> nPixels;
	itk::ImageRegionConstIteratorWithIndex<LabelType> labelIter(
                                                                _fixedImageLabel, _fixedImageLabel->GetBufferedRegion());

	int maxLabel = 0;
	for (labelIter.GoToBegin(); !labelIter.IsAtEnd(); ++labelIter) {
		LabelType::PixelType label = labelIter.Get();
		if (label > 0) {
            if (label > maxLabel) {
                _centerOfRotation.resize(label + 1);
                _labelIndexes.resize(label + 1);
                nPixels.resize(label + 1);
                for (int k = maxLabel; k < label + 1; k++) {
                    _centerOfRotation[k].SetSize(ImageType::ImageDimension);
                    _centerOfRotation[k].Fill(0);
                    nPixels[k] = 0;
                }
                maxLabel = label;
            }
			LabelType::IndexType idx = labelIter.GetIndex();
			_labelIndexes[label].push_back(idx);
            //cout << _centerOfRotation[label] << endl;
			for (int i = 0; i < ImageType::ImageDimension; i++) {
				_centerOfRotation[label][i] = _centerOfRotation[label][i] + idx[i];
			}
			nPixels[label]++;
            // cout << _centerOfRotation[label] << ", " << nPixels[label] << endl;
		}
	}
    cout << _fixedImage->GetBufferedRegion().GetSize() << endl;
	for (int label = 1; label <= maxLabel; label++) {
        itk::ContinuousIndex<double,ImageType::ImageDimension> centerIndex;
        ImageType::PointType centerPoint;
        for (int i = 0; i < ImageType::ImageDimension; i++) {
            centerIndex[i] = _centerOfRotation[label][i] / nPixels[label];
        }

        _fixedImage->TransformContinuousIndexToPhysicalPoint(centerIndex, centerPoint);
        cout << centerIndex << centerPoint << endl;
        for (int i = 0; i < ImageType::ImageDimension; i++) {
            _centerOfRotation[label][i] = centerPoint[i];
        }
	}
}

void itkMultiRegMethod::RunRegistration() {
    ComputeLabelIndexes();

	cout << "Registration starting ..." << endl;

	typedef itk::MyMetric<ImageType,ImageType> SubMetricType;
	typedef MetaMetrics<SubMetricType> MetaMetricType;

	MetaMetricType::Pointer _metaMetrics = MetaMetricType::New();
    OptimizerType::Pointer _optimizer = OptimizerType::New();

    int nSubParams = TransformType::New()->GetNumberOfParameters();

    QVector<ParametersType> params;
    for (unsigned int i = 0; i < _labelIndexes.size(); i++) {
        if (_labelIndexes[i].size() == 0) {
            continue;
        }
        cout << "# of indices: " << _labelIndexes[i].size() << endl;
        SubMetricType::Pointer metric = SubMetricType::New();
        InterpolatorType::Pointer interpolator = InterpolatorType::New();
        TransformType::Pointer transform = TransformType::New();
        metric->SetFixedImage(_fixedImage);
        metric->SetFixedImageIndexes(_labelIndexes[i]);
        transform->SetFixedParameters(_centerOfRotation[i]);
        cout << transform->GetMatrix() << "," << transform->GetOffset() << endl;
        metric->SetMovingImage(_movingImage);
        metric->SetInterpolator(interpolator);
        metric->SetTransform(dynamic_cast<SubMetricType::TransformType*>(transform.GetPointer()));
        metric->UseExactMatchOn();
        metric->Initialize();
        _metaMetrics->AddMetric(metric);
        params.push_back(transform->GetParameters());
    }

    MetaMetricType::ParametersType initialParams;
    initialParams.SetSize(_metaMetrics->GetNumberOfParameters());

    int k = 0;
    for (unsigned int i = 0; i < params.size(); i++) {
        for (int j = 0; j < params[i].GetSize(); j++) {
            initialParams[k++] = params[i][j];
        }
    }

    cout << "Initial Parameters: " << initialParams << endl;


    OptimizerType::ScalesType paramScales;
    paramScales.SetSize(_metaMetrics->GetNumberOfParameters());
    for (int i = 0; i < paramScales.GetSize(); i += nSubParams) {
        paramScales[i] = 30;
        paramScales[i+1] = 30;
        paramScales[i+2] = 30;
        paramScales[i+3] = .5;
        paramScales[i+4] = .5;
        paramScales[i+5] = .5;
        paramScales[i+6] = 160;
        paramScales[i+7] = 160;
        paramScales[i+8] = 160;
    }
    _optimizer->SetScales(paramScales);
    _optimizer->SetCostFunction(_metaMetrics);
    _optimizer->SetInitialPosition(initialParams);
    _optimizer->AddObserver(IterationEvent(), this);
    _optimizer->StartOptimization();

    ParametersType param = _metaMetrics->GetParameters();
    _transformHistory.push_back(param);

    SaveParameterHistory("TransformOutput.txt");
}

LabelType::Pointer itkMultiRegMethod::TransformFixedLabel(int historyId) {
    ParametersType param = _transformHistory[historyId];

    TransformType::Pointer transform = TransformType::New();
    int nParams = transform->GetNumberOfParameters();
    int nTotalParams = param.GetSize();
    int nTransforms = nTotalParams / nParams;

    typedef itkMyTransformFilter<LabelType, LabelType, TransformType::InverseTransformBaseType> MyTransformType;
    //typedef itkMyTransformFilter<LabelType, LabelType, TransformType> MyTransformType;
    MyTransformType::Pointer multiTransform = MyTransformType::New();
    
    for (int i = 0; i < nTransforms; i++) {
        TransformType::Pointer transform = TransformType::New();
        ParametersType p;
        p.SetSize(nParams);
        for (int j = 0; j < nParams; j++) {
            p[j] = param[i*nParams+j];
        }
        transform->SetParameters(p);
        ParametersType c;
        c.SetSize(3);
        for (int j = 0; j < 3; j++) {
            c[j] = _centerOfRotation[i+1][j];
        }
        transform->SetFixedParameters(c);
        multiTransform->AddTransform(transform->GetInverseTransform(), (i+1));
    }


    multiTransform->SetInputMask(_fixedImageLabel);
    multiTransform->SetInput(_fixedImageLabel);
    multiTransform->SetReferenceImage(_fixedImageLabel);
    multiTransform->SetInterpolator(MyTransformType::NNInterpolatorType::New());
    multiTransform->Update();
    LabelType::Pointer transformedOutput = multiTransform->GetOutput();

    return transformedOutput;
}



void itkMultiRegMethod::WriteTransform(const char* filename, int historyId) {

}

void itkMultiRegMethod::LoadParameterHistory(const char *filename) {
    char buff[1024];
    ifstream fin(filename);

    int nCenterOfRotations = 0;
    int nTransforms = 0;
    fin >> nCenterOfRotations;
    fin >> nTransforms;

    // to remove tailing newline
    fin.getline(buff, sizeof(buff));

    cout << "Reading " << nCenterOfRotations << " labels with " << nTransforms << " transforms" << endl;

    _centerOfRotation.clear();
    _centerOfRotation.push_back(ParametersType());
    for (int i = 0; i < nCenterOfRotations; i++) {
        ParametersType p;
        QVector<double> q;
        fin.clear();
        fin.getline(buff, sizeof(buff));
        if (!fin.eof()) {
            double x;
            cout << buff << endl;
            stringstream sin(buff);
            while (sin.good()) {
                sin >> x;
                if (sin.good()) {
                    q.push_back(x);
                }
            }
            p.SetSize(q.size());
            for (int i = 0; i < q.size(); i++) {
                p[i] = q[i];
            }
            _centerOfRotation.push_back(p);
        }
    }

    _transformHistory.clear();
    for (int i = 0; i < nTransforms; i++) {
        ParametersType p;
        QVector<double> q;
        fin.getline(buff, sizeof(buff));
        if (fin.good()) {
            double x;
            stringstream sin(buff);
            while (sin.good()) {
                sin >> x;
                if (sin.good()) {
                    q.push_back(x);
                }
            }
            p.SetSize(q.size());
            for (int i = 0; i < q.size(); i++) {
                p[i] = q[i];
            }
            _transformHistory.push_back(p);
        }
    }

    cout << _centerOfRotation.size() << " center of rotations read" << endl;
    cout << _transformHistory.size() << " transforms read" << endl;
}

std::string itkMultiRegMethod::GetTransformString(int n) {
    stringstream ss;
    ss << "Fixed Parameters: ";
    ss << (_centerOfRotation.size() - 1) << " " << _transformHistory.size() << endl;
    for (int i = 1; i < _centerOfRotation.size(); i++) {
        for (int j = 0; j < _centerOfRotation[i].GetSize(); j++) {
            ss << _centerOfRotation[i][j] << " ";
        }
        ss << endl;
    }
    ss << "Parameters: ";
    for (int j = 0; j < _transformHistory[n].GetSize(); j++) {
        ss << _transformHistory[n][j] << " ";
    }
    ss << endl;
    return ss.str();
}

void itkMultiRegMethod::SaveParameterHistory(const char *filename) {
    ofstream of(filename);
    of << (_centerOfRotation.size() - 1) << " " << _transformHistory.size() << endl;
    for (int i = 1; i < _centerOfRotation.size(); i++) {
        for (int j = 0; j < _centerOfRotation[i].GetSize(); j++) {
            of << _centerOfRotation[i][j] << " ";
        }
        of << endl;
    }

    for (int i = 0; i < _transformHistory.size(); i++) {
        for (int j = 0; j < _transformHistory[i].GetSize(); j++) {
            of << _transformHistory[i][j] << " ";
        }
        of << endl;
    }
    of.close();
}

int itkMultiRegMethod::GetNumberOfTransforms() {
    return _transformHistory.size();
}
