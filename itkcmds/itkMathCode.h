//
//  itkMathCode.h
//  itkcmds
//
//  Created by Joohwi Lee on 9/13/12.
//
//

#ifndef itkcmds_itkMathCode_h
#define itkcmds_itkMathCode_h

#include "MatrixCode.h"
#include "itkMatrixOffsetTransformBase.h"

template <typename TImage, typename TTransform>
class itkMathCode {
public:
    typedef MathCode::Mat4 Matrix;
    typedef MathCode::Vec4 Vector;
    typedef typename TTransform::Pointer TransformPointer;


    void transformIndex(Matrix& transform, typename TImage::IndexType& p, typename TImage::IndexType&  q) {
        Vector a(p[0], p[1], p[2], 1), b;
        MathCode::mult(transform, a, b);
        q[0] = b[0];
        q[1] = b[1];
        q[2] = b[2];
    }

    template <typename T>
    void transformContinousIndex(Matrix& transform, typename TImage::IndexType& p, T& q) {
        Vector a(p[0], p[1], p[2], 1), b;
        MathCode::mult(transform, a, b);
        q[0] = b[0];
        q[1] = b[1];
        q[2] = b[2];
    }

    void convertTransformToMatrix(TransformPointer transform, Matrix& matOut) const {
        typedef itk::MatrixOffsetTransformBase<double,3,3> MatrixTransformType;
         MatrixTransformType::Pointer matrixTransform = MatrixTransformType::Pointer(dynamic_cast<MatrixTransformType*>(transform.GetPointer()));

        if (matrixTransform.IsNull()) {
            cout << "can't convert to matrix-based transform: " << endl << transform << endl;
            matOut.identity();
            return;
        }
        
        typename MatrixTransformType::MatrixType itkMat = matrixTransform->GetMatrix();
        typename MatrixTransformType::ParametersType centerOfRotation = transform->GetFixedParameters();


        Matrix fp, fpback, tp, tmp1;
        for (int i = 0; i < 3; i++) {
            fp[i][3] -= centerOfRotation[i];
            fpback[i][3] += centerOfRotation[i];
            for (int j = 0; j < 3; j++) {
                tp[i][j] = itkMat[i][j];
            }
        }

        MathCode::mult(fpback, tp, tmp1);
        MathCode::mult(tmp1, fp, matOut);

        typename MatrixTransformType::OutputVectorType translation = matrixTransform->GetTranslation();
        for (int i = 0; i < 3; i++) {
			matOut[i][3] += translation[i];
        }

        cout << "Center of rotation: " << centerOfRotation << endl;
        cout << "Physical Transform: " << tp << endl;
        cout << "Translation: " << translation << endl;
        cout << "Final Transform: " << matOut << endl;
    }


    void createMat4FromMat(const typename TImage::DirectionType &mat,
                           MathCode::Mat4& matOut) const {
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) {
                matOut[r][c] = mat[r][c];
            }
        }
    }

    void createMat4FromVec(const typename TImage::SpacingType &spacing,
                           MathCode::Mat4& matOut, bool inverse =  false) const {
        for (int i = 0; i < 3; i++) {
            matOut[i][i] = spacing[i];
            if (inverse) {
                matOut[i][i] = 1 / matOut[i][i];
            }
        }
    }

    void createIndex2IndexTransformMatrix(typename TImage::Pointer src, TransformPointer transform,
                                          typename TImage::Pointer ref, MathCode::Mat4& matOut) const {
        MathCode::Mat4 srcDirection, srcSpacing, srcImageToWorld;
        createMat4FromMat(src->GetDirection(), srcDirection);
        createMat4FromVec(src->GetSpacing(), srcSpacing);

        MathCode::mult(srcDirection, srcSpacing, srcImageToWorld);
        for (int i = 0; i < 3; i++) {
            srcImageToWorld[i][3] += src->GetOrigin()[i];
        }

        MathCode::Mat4 srcToRefWorld, srcTransform;
        convertTransformToMatrix(transform, srcTransform);

        MathCode::mult(srcTransform, srcImageToWorld, srcToRefWorld);

        MathCode::Mat4 refDirection, refSpacing, refSpacingDirection, srcImageToRefImage;
        createMat4FromVec(ref->GetSpacing(), refSpacing, true);
        createMat4FromMat(ref->GetInverseDirection(), refDirection);

        MathCode::mult(refSpacing, refDirection, refSpacingDirection);
        MathCode::mult(refSpacingDirection, srcToRefWorld, srcImageToRefImage);
        for (int i = 0; i < 3; i++) {
            srcToRefWorld[i][3] -= ref->GetOrigin()[i];
        }
        matOut.copyFrom(srcImageToRefImage);
    }
};

#endif
