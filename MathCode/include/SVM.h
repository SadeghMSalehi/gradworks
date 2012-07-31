//
//  SVM.h
//  MathCode
//
//  Created by Joohwi Lee on 7/30/12.
//
//

#ifndef MathCode_SVM_h
#define MathCode_SVM_h

#include "MatrixCode.h"

namespace MathCode {

template <int N>
class SVM {
public:
    typedef DVec<float,N> InputType;
    
    InputType _dataSet;
    int _response[N];
    InputType _support;
    InputType _bias;
    
    void setDataSet(InputType* data) {
        _dataSet.copyFrom(data->_V);
    }
    bool validate() {
        int positiveCount = 0;
        for (int i = 0; i < N; i++) {
            if (_response[i]*(_support[i]*_dataSet[i]-_bias[i]) > 0){
                positiveCount ++;
            }
        }
        return positiveCount == N;
    }
};

}
#endif
