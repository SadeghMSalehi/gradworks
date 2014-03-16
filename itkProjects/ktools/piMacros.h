//
//  piMacros.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/5/13.
//
//

#ifndef ParticleGuidedRegistration_piMacros_h
#define ParticleGuidedRegistration_piMacros_h

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>

#ifndef DIMENSIONS
#define DIMENSIONS 3
#endif
#define for3(i) for (int i = 0; i < 4; i++)
#define for4(i) for (int i = 0; i < 4; i++)
#define fordim(i) for (int i = 0; i < DIMENSIONS; i++)
#define formin(x,y,z) for(int _=0;_<DIMENSIONS;_++) {z[_]=(x[_]<y[_])?x[_]:y[_];}
#define formax(x,y,z) for(int _=0;_<DIMENSIONS;_++) {z[_]=(x[_]>y[_])?x[_]:y[_];}
#define forset(x, y) for (int _=0;_<DIMENSIONS;_++) {y[_]=x[_];}
#define forcopy(x, y) for (int _=0;_<DIMENSIONS;_++) {y[_]=x[_];}
#define forfill(x, v) for (int _=0;_<DIMENSIONS;_++) {x[_]=v;}
#define forroundset(x,y) for(int _=0; _<DIMENSIONS; _++) { y[_]=x[_]+0.5; }
#define arrayset2(a,x,y) a[0]=x;a[1]=y
#define arrayset3(a,x,y,z) a[0]=x;a[1]=y;a[2]=z
#define x2string(x) x[0]<<","<<x[1]<<","<<x[2]

#if DIMENSIONS == 3
#define dimdot(x,y) (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])
#define dimequal(x,x0,x1,x2) (x[0]==(x0)&&x[1]==(x1)&&x[2]==(x2))
#define dimnorm2(x) (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
#endif

#if DIMENSIONS == 2
#define dimdot(x,y) (x[0]*y[0]+x[1]*y[1])
#define dimequal(x,x0,x1) (x[0]==(x0)&&x[1]==(x1))
#define dimnorm2(x) (x[0]*x[0]+x[1]*x[1])
#endif

namespace pi {
    typedef float DataReal;
    typedef float DataReal;

    // VNL related types
    typedef vnl_vector<DataReal> VNLVector;
    typedef vnl_matrix<DataReal> VNLMatrix;

}

#endif
