/*
 * HomographyTransform.h
 *
 *  Created on: Jul 28, 2012
 *      Author: joohwile
 */

#ifndef HOMOGRAPHYTRANSFORM_H_
#define HOMOGRAPHYTRANSFORM_H_

#include <iostream>
#include "MathCode.h"
#include "Image.h"

namespace cmath {

union TransformParamType {
    struct {
        float tx;
        float ty;
        float theta;
        float scale;
        float phi;
        float ratio;
        float v1;
        float v2;
    };
    float _x[8];
    TransformParamType() {
        tx = ty = theta = scale = phi = ratio = v1 = v2 = 0;
    }
    TransformParamType(TransformParamType& o) {
        memcpy(_x, o._x, sizeof(_x));
    }
    TransformParamType(float x, float y, float t) {
        tx = x;
        ty = y;
        theta = t;
        scale = phi = ratio = v1 = v2 = 0;
    }
    TransformParamType(float x, float y, float t, float s) {
        tx = x;
        ty = y;
        theta = t;
        scale = s;
        phi = ratio = v1 = v2 = 0;
    }
    TransformParamType(const TransformParamType& p) {
        for (int i = 0; i < 8; i++) {
            _x[i] = p._x[i];
        }
    }
    void coefficient(float a, TransformParamType& c) {
        for (int i = 0; i < 8; i++) {
            c._x[i] = a * _x[i];
        }
    }
    void linearCombination(float a, TransformParamType& b, TransformParamType& c) {
        for (int i = 0; i < 8; i++) {
            c._x[i] = _x[i] + a * b._x[i];
        }
    }
    float& operator[](int n) {
        return _x[n];
    }
    void negate() {
        for (int i = 0; i < 8; i++) {
            _x[i] = -_x[i];
        }
    }
    void multiply(float s) {
        for (int i = 0; i < 8; i++) {
            _x[i] *= s;
        }
    }
    void normal(TransformParamType& out) {
        double m = 0;
        for (int i = 0; i < 8; i++) {
            m += (_x[i] * _x[i]);
        }
        m = sqrt(m);
        for (int i = 0; i < 8; i++) {
            out._x[i] = _x[i] / m;
        }
    }
    void round(TransformParamType& out) {
        out.tx = cmath::round(tx);
        out.ty = cmath::round(ty);
        out.theta = cmath::round(theta);
    }
    float thetaRad() {
        return theta * PI / 180.0f;
    }
    void copyFrom(TransformParamType& in) {
        for (int i = 0; i < 8; i++) {
            _x[i] = in._x[i];
        }
    }
};
std::ostream& operator<<(std::ostream& os, const TransformParamType& param);

class Isometry {
public:
    void createAxisTransform(TransformParamType& param, Mat3& out) {
        float theta = param._x[2] * PI / 180.f;
        float tx = param._x[0];
        float ty = param._x[1];
        out._d[0] = cos(theta);
        out._d[1] = sin(theta);
        out._d[2] = -tx;
        out._d[3] = -sin(theta);
        out._d[4] = cos(theta);
        out._d[5] = -ty;
        out._d[6] = out._d[7] = 0;
        out._d[8] = 1;
    }
    void createPointTransform(TransformParamType& param, Mat3& out) {
        float theta = param._x[2] * PI / 180.f;
        float tx = param._x[0];
        float ty = param._x[1];
        out._d[0] = cos(theta);
        out._d[1] = -sin(theta);
        out._d[2] = tx;
        out._d[3] = sin(theta);
        out._d[4] = cos(theta);
        out._d[5] = ty;
        out._d[6] = out._d[7] = 0;
        out._d[8] = 1;
    }
};

class RotationTransform {
public:
    void createEulerAngleTransform(float x, float y, float z, Mat4& out) {
        float rx = x * PI / 180.0f;
        float ry = y * PI / 180.0f;
        float rz = z * PI / 180.0f;
        float cx = cos(rx), sx = sin(rx);
        float cy = cos(ry), sy = sin(ry);
        float cz = cos(rz), sz = sin(rz);
        out.identity();
        out._d[0] = cy * cz;
        out._d[1] = -cx * sz + sx * sy * cz;
        out._d[2] = sx * sz + cx * sy * cz;
        out._d[4] = cy * sz;
        out._d[5] = cx * cz + sx * sy * sz;
        out._d[6] = -sx * cz + cx * sy * sz;
        out._d[8] = -sy;
        out._d[9] = sx * cy;
        out._d[10] = cx * cy;
    }
};

class SimilarityTransform {
public:
    void createAxisTransform(TransformParamType& param, Mat3& out) {
        float scale = 100.f / param._x[3];
        float theta = param._x[2] * PI / 180.f;
        float tx = param._x[0];
        float ty = param._x[1];
        out._d[0] = scale * cos(theta);
        out._d[1] = sin(theta);
        out._d[2] = -tx;
        out._d[3] = -sin(theta);
        out._d[4] = scale * cos(theta);
        out._d[5] = -ty;
        out._d[6] = out._d[7] = 0;
        out._d[8] = 1;
    }
    void createPointTransform(TransformParamType& param, Mat3& out) {
        float scale = 1; //param._x[3] / 100.f;
        float theta = param._x[2] * PI / 180.f;
        float tx = param._x[0];
        float ty = param._x[1];
        out._d[0] = scale * cos(theta);
        out._d[1] = -sin(theta);
        out._d[2] = tx;
        out._d[3] = sin(theta);
        out._d[4] = scale * cos(theta);
        out._d[5] = ty;
        out._d[6] = out._d[7] = 0;
        out._d[8] = 1;
    }
};

class AffineTransform {
public:
    void createAxisTransform(TransformParamType& param, Mat3& out) {
        float phi = param._x[4];
        float scale = param._x[3] / 100.f;
        float theta = param._x[2];
        float tx = param._x[0];
        float ty = param._x[1];
        out._d[0] = scale * cos(theta);
        out._d[1] = sin(theta);
        out._d[2] = -tx;
        out._d[3] = -sin(theta);
        out._d[4] = scale * cos(theta);
        out._d[5] = -ty;
        out._d[6] = out._d[7] = 0;
        out._d[8] = 1;
    }
    void createPointTransform(TransformParamType& param, Mat3& out) {
        float phi = param._x[4];
        float scale = param._x[3] / 100.f;
        float theta = param._x[2];
        float tx = param._x[0];
        float ty = param._x[1];
        out._d[0] = scale * cos(theta);
        out._d[1] = sin(theta);
        out._d[2] = -tx;
        out._d[3] = -sin(theta);
        out._d[4] = scale * cos(theta);
        out._d[5] = -ty;
        out._d[6] = out._d[7] = 0;
        out._d[8] = 1;
    }
};

class PerspectiveTransform {
    void createPerspectiveTransform(Vec4 cameraPos, Mat4& out) {
        out.identity();
        out._d[3] = -cameraPos._x[0];
        out._d[7] = -cameraPos._x[1];
        out._d[14] = 1 / cameraPos._x[2];
    }
};

class ImageTransform {
public:
    static void createImageToNormalizedCoordinate(int w, int h, Mat4& out) {
        out.identity();
        out._d[0] = 2.0f / w;
        out._d[3] = -1;
        out._d[5] = -2.0f / h;
        out._d[7] = 1;
    }

    static void createNormalizedToImageCoordinate(int w, int h, Mat4& out) {
        out.identity();
        out._d[0] = w / 2.0f;
        out._d[3] = w / 2.0f;
        out._d[5] = -h / 2.0f;
        out._d[7] = h / 2.0f;
    }

    // simple perspecitve projection
    static void createPerspectiveProjection(Mat4& out) {
        out.identity();
        out._d[14] = 1;
        //out._d[15] = 0;
    }

    static void createEulerAngleTransform(float x, float y, float z, Mat4& out) {
        float rx = x * PI / 180.0f;
        float ry = y * PI / 180.0f;
        float rz = z * PI / 180.0f;
        float cx = cos(rx), sx = sin(rx);
        float cy = cos(ry), sy = sin(ry);
        float cz = cos(rz), sz = sin(rz);
        out.identity();
        out._d[0] = cy * cz;
        out._d[1] = -cx * sz + sx * sy * cz;
        out._d[2] = sx * sz + cx * sy * cz;
        out._d[4] = cy * sz;
        out._d[5] = cx * cz + sx * sy * sz;
        out._d[6] = -sx * cz + cx * sy * sz;
        out._d[8] = -sy;
        out._d[9] = sx * cy;
        out._d[10] = cx * cy;
    }

    /*
     static void createSimilarityTransform(TransformParamType& param,
     Mat3& out) {
     float scale = param._x[3] / 100.f;
     float theta = param._x[2] * PI / 180.f;
     float tx = param._x[0];
     float ty = param._x[1];
     out._d[0] = scale * cos(theta);
     out._d[1] = -sin(theta);
     out._d[2] = tx;
     out._d[3] = sin(theta);
     out._d[4] = scale * cos(theta);
     out._d[5] = ty;
     out._d[6] = out._d[7] = 0;
     out._d[8] = 1;
     }
     */

    /**
     * transform ij-coordinate to xy-coordinate
     */
    static void createRigidTransformForImage(TransformParamType& param, int w, int h, Mat3& out) {
        float scale = 1;
        float theta = param._x[2] * PI / 180.f;
        float tx = param._x[0];
        float ty = param._x[1];
        float co = cos(theta);
        float si = sin(theta);
        float hw = w / 2.0f;
        float hh = h / 2.0f;
        out._d[0] = scale * co;
        out._d[1] = scale * -si;
        out._d[2] = scale * (-hw * co + hh * si) + tx;
        out._d[3] = scale * si;
        out._d[4] = scale * co;
        out._d[5] = scale * (-hw * si - hh * co) + ty;
        out._d[6] = out._d[7] = 0;
        out._d[8] = 1;
    }

    /**
     * transform ij-coordinate to xy-coordinate
     */
    static void createSimilarityTransformForImage(TransformParamType& param, int w, int h, Mat3& out) {
        float scale = param._x[3] / 100.f;
        float theta = param._x[2] * PI / 180.f;
        float tx = param._x[0];
        float ty = param._x[1];
        float co = cos(theta);
        float si = sin(theta);
        float hw = w / 2.0f;
        float hh = h / 2.0f;
        out._d[0] = scale * co;
        out._d[1] = scale * -si;
        out._d[2] = scale * (-hw * co + hh * si) + tx;
        out._d[3] = scale * si;
        out._d[4] = scale * co;
        out._d[5] = scale * (-hw * si - hh * co) + ty;
        out._d[6] = out._d[7] = 0;
        out._d[8] = 1;
    }

    static void asImageCenter(Mat3& in, int w, int h) {
        in._d[2] -= (w / 2);
        in._d[5] -= (h / 2);
    }
};

class Homography {
    void createAxisTransform(TransformParamType& param, Mat3& tfm) {

    }
    void createPointTransform(TransformParamType& param, Mat3& tfm) {

    }
};

class CImageResample {
public:
    void rotateFrame(Mat3& rot) {

    }

    bool resampleIJ2XY(const IntImage& src, const IntImage& obj, Mat3 fTx, IntImage& dst) {
        fTx.suppressEpsilion();
        if (!dst.reuseOrCreateEmptyOf(obj.getWidth(), obj.getHeight())) {
            return false;
        }

        Vec2 ij, xy;
        for (ij[1] = 0; ij[1] < obj.getHeight(); ij[1]++) {
            for (ij[0] = 0; ij[0] < obj.getWidth(); ij[0]++) {
                mult(fTx, ij, xy);
                if (src.isInside(xy[0], xy[1])) {
                    dst(ij[0], ij[1]) = src.bilinear(xy[0], xy[1]);
                } else {
                    dst(ij[0], ij[1]) = 0;
                }
            }
        }
        return true;
    }

    template<typename SrcType, typename DstType>
    bool resampleIJ2XYinReal(const CImage<SrcType>& src, const CImage<SrcType>& obj, Mat3 fTx, CImage<DstType>& dst) {
        fTx.suppressEpsilion();
        if (!dst.reuseOrCreateEmptyOf(obj.getWidth(), obj.getHeight())) {
            cout << "Fail to resample" << endl;
            //throw "Resampling Error: destiation not available";
            return false;
        }

        Vec2 ij, xy;
        for (ij[1] = 0; ij[1] < obj.getHeight(); ij[1]++) {
            for (ij[0] = 0; ij[0] < obj.getWidth(); ij[0]++) {
                mult(fTx, ij, xy);
                dst(ij[0], ij[1]) = static_cast<DstType>(src.safe_bilinear(xy[0], xy[1]));
            }
        }
        return true;
    }

    void extract(const IntImage& src, const IntImage& obj, Mat3 fTx, IntImage& dst) {
        fTx.suppressEpsilion();
        dst.createEmptyOf(obj.getWidth(), obj.getHeight());

        Vec2 ij, xy;
        for (xy[1] = 0; xy[1] < obj.getHeight(); xy[1]++) {
            for (xy[0] = 0; xy[0] < obj.getWidth(); xy[0]++) {
                mult(fTx, xy, ij);
                if (src.isInside(ij[0], ij[1])) {
                    dst(xy[0], xy[1]) = src(ij[0], ij[1]);
                } else {
                    dst(xy[0], xy[1]) = 0;
                }
            }
        }
    }

    void resample(IntImage& src, Mat3& fTx, IntImage& dst) {
        Vec2 p[4], q[4], ro, wh;
        p[0].set(0, 0);
        p[1].set(src.getWidth(), 0);
        p[2].set(0, src.getHeight());
        p[3].set(src.getWidth(), src.getHeight());
        mult(fTx, p, q, 4);

        ro = vecmin(q, 4);
        wh = vecmax(q, 4) - ro;

        fTx.translateOrigin(ro);

        cout << fTx << endl;
        Mat3 bTx;
        fTx.inverse(bTx);

        cout << bTx << endl;
        bTx.suppressEpsilion();

        dst.createEmptyOf(wh[0], wh[1]);

        Vec2 ij, xy;
        for (xy[1] = 0; xy[1] < wh[1]; xy[1]++) {
            for (xy[0] = 0; xy[0] < wh[0]; xy[0]++) {
                mult(bTx, xy, ij);
                if (src.isInside(ij[0], ij[1])) {
                    dst(xy[0], xy[1]) = src(ij[0], ij[1]);
                } else {
                    dst(xy[0], xy[1]) = 0;
                }
            }
        }
    }

    void resampleH(IntImage& src, Mat3& fTx, IntImage& dst) {
        Vec3 p[4], q[4];
        p[0].set(0, 0, 1);
        p[1].set(src.getWidth(), 0, 1);
        p[2].set(0, src.getHeight(), 1);
        p[3].set(src.getWidth(), src.getHeight(), 1);
        mult(fTx, p, q, 4);
        scaleZ(q, 4, q);

        for (int i = 0; i < 4; i++) {
            cout << q[i] << endl;
        }

        Vec3 topleft = vecmin(q, 4);
        Vec3 bottomright = vecmax(q, 4);
        Vec3 wh = bottomright - topleft;

        cout << "Bounding box: " << topleft << " ~ " << bottomright << endl;
        cout << "WH: " << wh << endl;
        cout << "Mat: " << fTx << endl;

        fTx.translateOrigin(topleft);

        cout << fTx << endl;
        Mat3 bTx;
        fTx.inverse(bTx);

        cout << bTx << endl;
        bTx.suppressEpsilion();

        dst.createEmptyOf(wh[0], wh[1]);

        Vec2 ij, xy;
        for (xy[1] = 0; xy[1] < wh[1]; xy[1]++) {
            for (xy[0] = 0; xy[0] < wh[0]; xy[0]++) {
                mult(bTx, xy, ij);
                if (src.isInside(ij[0], ij[1])) {
                    dst(xy[0], xy[1]) = src(ij[0], ij[1]);
                } else {
                    dst(xy[0], xy[1]) = 0;
                }
            }
        }
    }

    void resampleH(IntImage& src, Mat4& fTx, IntImage& dst) {
        Vec4 ext(src.getWidth(), src.getHeight(), 0, 1), extX;
        Vec4 ndc(-1, -1, 0, 1);

        Mat4 bTx;
        fTx.inverse(bTx);

        mult(bTx, ndc, extX);
        cout << "Extent: " << extX << endl;

        Mat4 view, txf;
        ImageTransform::createImageToNormalizedCoordinate(dst.getWidth(), dst.getHeight(), view);
        mult(bTx, view, txf);

        for (int j = 0; j < dst.getHeight(); j++) {
            for (int i = 0; i < dst.getWidth(); i++) {
                Vec4 ij(i, j, 0, 1), xy;
                mult(txf, ij, xy);
                if (src.isInside(xy.x(), xy.y())) {
                    dst.set(i, j, src.bilinear(xy.x(), xy.y()));
                }
            }
        }
    }

    void projectH(IntImage& src, Mat4& fTx, IntImage& dst) {
        Mat4 view, txf;
        ImageTransform::createNormalizedToImageCoordinate(dst.getWidth(), dst.getHeight(), view);
        mult(view, fTx, txf);

        for (int j = 0; j < src.getHeight(); j++) {
            for (int i = 0; i < src.getWidth(); i++) {
                Vec4 ij(i, j, 0, 1), xy, xy2;
                mult(txf, ij, xy2);
                xy2.homogenize();
                //mult(view, xy, xy2);
                if (dst.isInside(xy2.x(), xy2.y())) {
                    dst.set(xy2.x(), xy2.y(), src(i, j));
                }
            }
        }
    }

    void resampleTwoPass(IntImage& in, Mat4 M, IntImage& one, IntImage& two) {
        one.createEmptyOf(in.getWidth(), in.getHeight());
        for (int v = 0; v < one.getHeight(); v++) {
            for (int u = 0; u < one.getWidth(); u++) {
                float A = M[0][0];
                float B = M[0][1] * v + M[0][3];
                float C = M[3][0];
                float D = M[3][1] * v + M[3][3];
                float x = (A * u + B) / (C * u + D);
                float y = v;
                if (in.isInside(round(x), round(y))) {
                    one.set(u, v, in.bilinear(x, y));
                }
            }
        }

        two.createEmptyOf(one.getWidth(), one.getHeight());
        for (int v = 0; v < two.getHeight(); v++) {
            for (int x_ = 0; x_ < two.getWidth(); x_++) {
                float x = x_;
                float E = (M[0][1] - M[3][1] * x) / (M[3][0] * x - M[0][0]);
                float F = (M[0][3] - M[3][3] * x) / (M[3][0] * x - M[0][0]);
                float G = M[1][1] + M[1][0] * E;
                float H = M[1][3] + M[1][0] * F;
                float I = M[3][1] + M[3][0] * E;
                float J = M[3][3] + M[3][0] * F;
                float y = (G * v + H) / (I * v + J);
                if (one.isInside(round(x), round(y))) {
                    two.set(x_, v, one.bilinear(x, y));
                }
            }
        }
    }
};
}

#endif /* HOMOGRAPHYTRANSFORM_H_ */
