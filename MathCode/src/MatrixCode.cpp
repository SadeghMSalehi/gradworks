/*
 * CMath.cpp
 *
 *  Created on: Jul 21, 2012
 *      Author: joohwile
 */

#include "MatrixCode.h"
#include "math.h"

namespace MathCode {

float abs(float n) {
    return n > 0 ? n : -n;
}
int round(float n) {
    return (int) (n > 0 ? (n + .5) : (n - .5));
}
float sqrt(float s) {
    return ::sqrt(s);
}
float cos(float x) {
    return ::cos(x);
}
float sin(float x) {
    return ::sin(x);
}
Vec2 max(Vec2* a, int n) {
    Vec2 o = a[0];
    for (int i = 0; i < n; i++) {
        o._V[0] = o._V[0] < a[i]._V[0] ? a[i]._V[0] : o._V[0];
        o._V[1] = o._V[1] < a[i]._V[1] ? a[i]._V[1] : o._V[1];
    }
    return o;
}

Vec3 max(Vec3* a, int n) {
    Vec3 o = a[0];
    for (int i = 0; i < n; i++) {
        o._V[0] = o._V[0] < a[i]._V[0] ? a[i]._V[0] : o._V[0];
        o._V[1] = o._V[1] < a[i]._V[1] ? a[i]._V[1] : o._V[1];
        o._V[2] = o._V[2] < a[i]._V[2] ? a[i]._V[2] : o._V[2];
    }
    return o;
}

Vec4 max(Vec4* a, int n) {
    Vec4 o = a[0];
    for (int i = 0; i < n; i++) {
        o._V[0] = o._V[0] < a[i]._V[0] ? a[i]._V[0] : o._V[0];
        o._V[1] = o._V[1] < a[i]._V[1] ? a[i]._V[1] : o._V[1];
        o._V[2] = o._V[2] < a[i]._V[2] ? a[i]._V[2] : o._V[2];
        o._V[3] = o._V[3] < a[i]._V[3] ? a[i]._V[3] : o._V[3];
    }
    return o;
}

Vec2 min(Vec2* a, int n) {
    Vec2 o = a[0];
    for (int i = 0; i < n; i++) {
        o._V[0] = o._V[0] > a[i]._V[0] ? a[i]._V[0] : o._V[0];
        o._V[1] = o._V[1] > a[i]._V[1] ? a[i]._V[1] : o._V[1];
    }
    return o;
}

Vec3 min(Vec3* a, int n) {
    Vec3 o = a[0];
    for (int i = 0; i < n; i++) {
        o._V[0] = o._V[0] > a[i]._V[0] ? a[i]._V[0] : o._V[0];
        o._V[1] = o._V[1] > a[i]._V[1] ? a[i]._V[1] : o._V[1];
        o._V[2] = o._V[2] > a[i]._V[2] ? a[i]._V[2] : o._V[2];
    }
    return o;
}

Vec4 min(Vec4* a, int n) {
    Vec4 o = a[0];
    for (int i = 0; i < n; i++) {
        o._V[0] = o._V[0] > a[i]._V[0] ? a[i]._V[0] : o._V[0];
        o._V[1] = o._V[1] > a[i]._V[1] ? a[i]._V[1] : o._V[1];
        o._V[2] = o._V[2] > a[i]._V[2] ? a[i]._V[2] : o._V[2];
        o._V[3] = o._V[3] > a[i]._V[3] ? a[i]._V[3] : o._V[3];
    }
    return o;
}

}
