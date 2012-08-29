/*
 * Marcos.h
 *
 *  Created on: Jul 27, 2012
 *      Author: joohwile
 */

#ifndef MARCOS_H_
#define MARCOS_H_

#define forN(i) for (int i = 0; i < N; i++)
#define forX(i,X) for (int i = 0; i < X; i++)
#define forXY(j,Y,i,X) for (int j = 0; j < Y; j++) for (int i = 0; i < X; i++)
#define __cmath_min(x,y) (x<y?x:y)
#define __cmath_max(x,y) (x>y?x:y)

#ifdef DEBUG
#define _DBG_(stmt) stmt
#else
#define _DBG_(stmt)
#endif
#endif /* MARCOS_H_ */
