/*
 An Algorithm for Automatically Fitting Digitized Curves
 by Philip J. Schneider
 from "Graphics Gems", Academic Press, 1990
 */

/*  fit_cubic.c	*/
/*	Piecewise cubic fitting code	*/


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "GraphicsGems.h"

#define MAXPOINTS	1000		/* Maximum number of points in the input polyline */
namespace gem {
//
//    typedef Point2 *BezierCurve;
//
//    /* Forward declarations */
//    
////    static	void		FitCubic();
////    static	double		*Reparameterize();
////    static	double		NewtonRaphsonRootFind();
////    static	Point2		BezierII();
////    static	double 		B0(), B1(), B2(), B3();
////    static	Vector2		ComputeLeftTangent();
////    static	Vector2		ComputeRightTangent();
////    static	Vector2		ComputeCenterTangent();
////    static	double		ComputeMaxError();
////    static	double		*ChordLengthParameterize();
////    static	BezierCurve	GenerateBezier();
////    static	Vector2		V2AddII();
////    static	Vector2		V2ScaleIII();
////    static	Vector2		V2SubII();
//    static	int		totalOutputPoints = 0;
//    static  Point2* 	bezierBuffer = 0;
//    static  int		pointsAllocated = 0;
//
//    static Vector2 V2SubII(Vector2 a, Vector2 b);
//    static Vector2 V2ScaleIII(Vector2 v, double s);
//    static Vector2 V2AddII(Vector2 a, Vector2 b);
//    static double ComputeMaxError(Point2* d, int first, int last, BezierCurve bezCurve, double* u, int* splitPoint);
//    static double *ChordLengthParameterize(Point2* d, int first, int last);
//    static Vector2 ComputeCenterTangent(Point2* d, int center);
//    static Vector2 ComputeRightTangent(Point2* d, int end);
//    static Vector2 ComputeLeftTangent(Point2 *d, int end);
//    static double B3(double u);
//    static double B2(double u);
//    static double B1(double u);
//    static double B0(double u);
//    static Point2 BezierII(int degree, Point2* V, double t);
//    static double NewtonRaphsonRootFind(BezierCurve Q, Point2 P, double u);
//    static double *Reparameterize(Point2* d, int first, int last, double* u, BezierCurve bezCurve);
//    static BezierCurve  GenerateBezier(Point2* d, int first, int last, double* uPrime, Vector2 tHat1, Vector2 tHat2);
//    static void FitCubic(Point2* d, int first, int last, Vector2 tHat1, Vector2 tHat2, double error);
//    int FitCurve(double* d, int nPts, double error);
//    void CopyResults(double* vector);
//    void OutputBezierCurve(int n, BezierCurve curve);
//    BezierCurve AllocBezier();
//    void DeallocLastBezier();
//
//    void DeallocLastBezier()
//    {
//        printf("DeallocLastBezier\n");
//        pointsAllocated -= 4;
//    }
//
//    BezierCurve AllocBezier()
//    {
//        printf("AllocBezier\n");
//        int offset = pointsAllocated;
//        pointsAllocated += 4;
//        if (!bezierBuffer)
//            bezierBuffer = (Point2 *)malloc(MAXPOINTS * 2 * sizeof(Point2));
//        if (pointsAllocated >= MAXPOINTS * 2)
//            return NULL;
//        return &(bezierBuffer[offset]);
//    }
//
//    void OutputBezierCurve(int n, BezierCurve curve)
//    {
//        printf("outcome Bezier curve: n %d; points: %lf, %lf; %lf, %lf; %lf, %lf; %lf, %lf\n",
//               n, curve[0].x, curve[0].y, curve[1].x, curve[1].y,
//               curve[2].x, curve[2].y, curve[3].x, curve[3].y);
//
//        totalOutputPoints += 8;
//    }
//
//    /**
//     Assume the given pointer has been pre-allocated (which was possible
//     since FitCurve returns the number of coords) and just memcpy it in.
//     Also frees the output bezierBuffer.
//     */
//    void CopyResults(double* vector)
//    {
//        printf("copying from %x to %x: %d bytes\n", bezierBuffer, vector,  totalOutputPoints * sizeof(double));
//        memcpy(vector, bezierBuffer, totalOutputPoints * sizeof(double));
//    }
//
//    /**
//     Freeing the bezierBuffer is optional, because if you are doing
//     multiple calls to FitCurve, it can be re-used.
//     */
//    void FreeResults()
//    {
//        if (bezierBuffer)
//            free(bezierBuffer);
//        bezierBuffer = 0;
//    }
//
//    /*
//     *  FitCurve :
//     *  	Fit a Bezier curve to a set of digitized points
//     */
//    int FitCurve(double* d, int nPts, double error)
//    {
//        Vector2	tHat1, tHat2;	/*  Unit tangent vectors at endpoints */
//        nPts /= 2; 	// we are treating an array of doubles as if it was an array of (x, y) structs
//        // the memory allocation is the same either way, but there are half as many points as coords
//
//        tHat1 = ComputeLeftTangent((Point2*) d, 0);
//        tHat2 = ComputeRightTangent((Point2*) d, nPts - 1);
//        FitCubic((Point2*) d, 0, nPts - 1, tHat1, tHat2, error);
//        printf("returning %d points at address %x\n", totalOutputPoints, (unsigned long)bezierBuffer);
//        return totalOutputPoints;
//    }
//
//    /* for testing the C code by itself
//     main()
//     {
//     static Point2 d[7] = {	//  Digitized points
//     { 0.0, 0.0 },
//     { 0.0, 0.5 },
//     { 1.1, 1.4 },
//     { 2.1, 1.6 },
//     { 3.2, 1.1 },
//     { 4.0, 0.2 },
//     { 4.0, 0.0 },
//     };
//     double	error = 4.0;		//  Squared error
//     int size = 0;
//     FitCurve(d, 7, error, &size);		//  Fit the Bezier curves
//     }
//     */
//
//    /*
//     *  FitCubic :
//     *  	Fit a Bezier curve to a (sub)set of digitized points
//     */
//    static void FitCubic(Point2* d, int first, int last, Vector2 tHat1, Vector2 tHat2, double error)
////    Point2	*d,			/*  Array of digitized points */
////    int		first, int last,	/* Indices of first and last pts in region */
////    Vector2	tHat1, Vector2 tHat2,	/* Unit tangent vectors at endpoints */
////    double	error)		/*  User-defined error squared	   */
//    {
//        BezierCurve	bezCurve; /*Control points of fitted Bezier curve*/
//        double	*u;		/*  Parameter values for point  */
//        double	*uPrime;	/*  Improved parameter values */
//        double	maxError;	/*  Maximum fitting error	 */
//        int		splitPoint;	/*  Point to split point set at	 */
//        int		nPts;		/*  Number of points in subset  */
//        double	iterationError; /*Error below which you try iterating  */
//        int		maxIterations = 4; /*  Max times to try iterating  */
//        Vector2	tHatCenter;   	/* Unit tangent vector at splitPoint */
//        int		i;
//
//        iterationError = error * error;
//        nPts = last - first + 1;
//
//        /*  Use heuristic if region only has two points in it */
//        if (nPts == 2) {
//            double dist = V2DistanceBetween2Points(&d[last], &d[first]) / 3.0;
//
//            bezCurve = AllocBezier();
//            bezCurve[0] = d[first];
//            bezCurve[3] = d[last];
//            V2Add(&bezCurve[0], V2Scale(&tHat1, dist), &bezCurve[1]);
//            V2Add(&bezCurve[3], V2Scale(&tHat2, dist), &bezCurve[2]);
//            OutputBezierCurve(3, bezCurve);
//            return;
//        }
//
//        /*  Parameterize points, and attempt to fit curve */
//        u = ChordLengthParameterize(d, first, last);
//        bezCurve = GenerateBezier(d, first, last, u, tHat1, tHat2);
//
//        /*  Find max deviation of points to fitted curve */
//        maxError = ComputeMaxError(d, first, last, bezCurve, u, &splitPoint);
//        if (maxError < error) {
//            OutputBezierCurve(3, bezCurve);
//            free((void *)u);
//            return;
//        }
//
//
//        /*  If error not too large, try some reparameterization  */
//        /*  and iteration */
//        if (maxError < iterationError) {
//            for (i = 0; i < maxIterations; i++) {
//                uPrime = Reparameterize(d, first, last, u, bezCurve);
//                DeallocLastBezier();
//                bezCurve = GenerateBezier(d, first, last, uPrime, tHat1, tHat2);
//                maxError = ComputeMaxError(d, first, last,
//                                           bezCurve, uPrime, &splitPoint);
//                if (maxError < error) {
//                    OutputBezierCurve(3, bezCurve);
//                    free((void *)u);
//                    free((void *)uPrime);
//                    return;
//                }
//                free((void *)u);
//                u = uPrime;
//            }
//        }
//
//        /* Fitting failed -- split at max error point and fit recursively */
//        free((void *)u);
//        DeallocLastBezier();
//        tHatCenter = ComputeCenterTangent(d, splitPoint);
//        FitCubic(d, first, splitPoint, tHat1, tHatCenter, error);
//        V2Negate(&tHatCenter);
//        FitCubic(d, splitPoint, last, tHatCenter, tHat2, error);
//    }
//
//
//    /*
//     *  GenerateBezier :
//     *  Use least-squares method to find Bezier control points for region.
//     Point2	*d;			/*  Array of digitized points	*/
////    int		first, last;		/*  Indices defining region	*/
////    double	*uPrime;		/*  Parameter values for region */
////    Vector2	tHat1, tHat2;	/*  Unit tangents at endpoints	*/
//    static BezierCurve  GenerateBezier(Point2* d, int first, int last, double* uPrime, Vector2 tHat1, Vector2 tHat2)
//    {
//        int 	i;
//        Vector2 	A[MAXPOINTS][2];	/* Precomputed rhs for eqn	*/
//        int 	nPts;			/* Number of pts in sub-curve */
//        double 	C[2][2];			/* Matrix C		*/
//        double 	X[2];			/* Matrix X			*/
//        double 	det_C0_C1,		/* Determinants of matrices	*/
//        det_C0_X,
//        det_X_C1;
//        double 	alpha_l,		/* Alpha values, left and right	*/
//        alpha_r;
//        Vector2 	tmp;			/* Utility variable		*/
//        BezierCurve	bezCurve;	/* RETURN bezier curve ctl pts	*/
//
//        bezCurve = AllocBezier();
//        nPts = last - first + 1;
//
//
//        /* Compute the A's	*/
//        for (i = 0; i < nPts; i++) {
//            Vector2		v1, v2;
//            v1 = tHat1;
//            v2 = tHat2;
//            V2Scale(&v1, B1(uPrime[i]));
//            V2Scale(&v2, B2(uPrime[i]));
//            A[i][0] = v1;
//            A[i][1] = v2;
//        }
//
//        /* Create the C and X matrices	*/
//        C[0][0] = 0.0;
//        C[0][1] = 0.0;
//        C[1][0] = 0.0;
//        C[1][1] = 0.0;
//        X[0]    = 0.0;
//        X[1]    = 0.0;
//
//        for (i = 0; i < nPts; i++) {
//            C[0][0] += V2Dot(&A[i][0], &A[i][0]);
//            C[0][1] += V2Dot(&A[i][0], &A[i][1]);
//            /*					C[1][0] += V2Dot(&A[i][0], &A[i][1]);*/
//            C[1][0] = C[0][1];
//            C[1][1] += V2Dot(&A[i][1], &A[i][1]);
//
//            tmp = V2SubII(d[first + i],
//                          V2AddII(
//                                  V2ScaleIII(d[first], B0(uPrime[i])),
//                                  V2AddII(
//                                          V2ScaleIII(d[first], B1(uPrime[i])),
//                                          V2AddII(
//                                                  V2ScaleIII(d[last], B2(uPrime[i])),
//                                                  V2ScaleIII(d[last], B3(uPrime[i]))))));
//
//
//            X[0] += V2Dot(&A[i][0], &tmp);
//            X[1] += V2Dot(&A[i][1], &tmp);
//        }
//
//        /* Compute the determinants of C and X	*/
//        det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
//        det_C0_X  = C[0][0] * X[1]    - C[1][0] * X[0];
//        det_X_C1  = X[0]    * C[1][1] - X[1]    * C[0][1];
//
//        /* Finally, derive alpha values	*/
//        alpha_l = (det_C0_C1 == 0) ? 0.0 : det_X_C1 / det_C0_C1;
//        alpha_r = (det_C0_C1 == 0) ? 0.0 : det_C0_X / det_C0_C1;
//
//        /* If alpha negative, use the Wu/Barsky heuristic (see text) */
//        /* (if alpha is 0, you get coincident control points that lead to
//         * divide by zero in any subsequent NewtonRaphsonRootFind() call. */
//        double segLength = V2DistanceBetween2Points(&d[last], &d[first]);
//        double epsilon = 1.0e-6 * segLength;
//        if (alpha_l < epsilon || alpha_r < epsilon)
//        {
//            /* fall back on standard (probably inaccurate) formula, and subdivide further if needed. */
//            double dist = segLength / 3.0;
//            bezCurve[0] = d[first];
//            bezCurve[3] = d[last];
//            V2Add(&bezCurve[0], V2Scale(&tHat1, dist), &bezCurve[1]);
//            V2Add(&bezCurve[3], V2Scale(&tHat2, dist), &bezCurve[2]);
//            return (bezCurve);
//        }
//
//        /*  First and last control points of the Bezier curve are */
//        /*  positioned exactly at the first and last data points */
//        /*  Control points 1 and 2 are positioned an alpha distance out */
//        /*  on the tangent vectors, left and right, respectively */
//        bezCurve[0] = d[first];
//        bezCurve[3] = d[last];
//        V2Add(&bezCurve[0], V2Scale(&tHat1, alpha_l), &bezCurve[1]);
//        V2Add(&bezCurve[3], V2Scale(&tHat2, alpha_r), &bezCurve[2]);
//        return (bezCurve);
//    }
//
//
//    /*
//     *  Reparameterize:
//     *	Given set of points and their parameterization, try to find
//     *   a better parameterization.
//     *
//     */
//    static double *Reparameterize(Point2* d, int first, int last, double* u, BezierCurve bezCurve)
////    Point2	*d;			/*  Array of digitized points	*/
////    int		first, last;		/*  Indices defining region	*/
////    double	*u;			/*  Current parameter values	*/
////    BezierCurve	bezCurve;	/*  Current fitted curve	*/
//    {
//        int 	nPts = last-first+1;
//        int 	i;
//        double	*uPrime;		/*  New parameter values	*/
//
//        uPrime = (double *)malloc(nPts * sizeof(double));
//        for (i = first; i <= last; i++) {
//            uPrime[i-first] = NewtonRaphsonRootFind(bezCurve, d[i], u[i-
//                                                                      first]);
//        }
//        return (uPrime);
//    }
//
//
//
//    /*
//     *  NewtonRaphsonRootFind :
//     *	Use Newton-Raphson iteration to find better root.
//     */
//    static double NewtonRaphsonRootFind(BezierCurve Q, Point2 P, double u)
////    BezierCurve	Q;			/*  Current fitted curve	*/
////    Point2 		P;		/*  Digitized point		*/
////    double 		u;		/*  Parameter value for "P"	*/
//    {
//        double 		numerator, denominator;
//        Point2 		Q1[3], Q2[2];	/*  Q' and Q''			*/
//        Point2		Q_u, Q1_u, Q2_u; /*u evaluated at Q, Q', & Q''	*/
//        double 		uPrime;		/*  Improved u			*/
//        int 		i;
//
//        /* Compute Q(u)	*/
//        Q_u = BezierII(3, Q, u);
//
//        /* Generate control vertices for Q'	*/
//        for (i = 0; i <= 2; i++) {
//            Q1[i].x = (Q[i+1].x - Q[i].x) * 3.0;
//            Q1[i].y = (Q[i+1].y - Q[i].y) * 3.0;
//        }
//
//        /* Generate control vertices for Q'' */
//        for (i = 0; i <= 1; i++) {
//            Q2[i].x = (Q1[i+1].x - Q1[i].x) * 2.0;
//            Q2[i].y = (Q1[i+1].y - Q1[i].y) * 2.0;
//        }
//
//        /* Compute Q'(u) and Q''(u)	*/
//        Q1_u = BezierII(2, Q1, u);
//        Q2_u = BezierII(1, Q2, u);
//
//        /* Compute f(u)/f'(u) */
//        numerator = (Q_u.x - P.x) * (Q1_u.x) + (Q_u.y - P.y) * (Q1_u.y);
//        denominator = (Q1_u.x) * (Q1_u.x) + (Q1_u.y) * (Q1_u.y) +
//        (Q_u.x - P.x) * (Q2_u.x) + (Q_u.y - P.y) * (Q2_u.y);
//        if (denominator == 0.0f) return u;
//
//        /* u = u - f(u)/f'(u) */
//        uPrime = u - (numerator/denominator);
//        return (uPrime);
//    }
//
//
//
//    /*
//     *  Bezier :
//     *  	Evaluate a Bezier curve at a particular parameter value
//     *
//     */
//    static Point2 BezierII(int degree, Point2* V, double t)
////    int		degree;		/* The degree of the bezier curve	*/
////    Point2 	*V;		/* Array of control points		*/
////    double 	t;		/* Parametric value to find point for	*/
//    {
//        int 	i, j;
//        Point2 	Q;	        /* Point on curve at parameter t	*/
//        Point2 	*Vtemp;		/* Local copy of control points		*/
//
//        /* Copy array	*/
//        Vtemp = (Point2 *)malloc((unsigned)((degree+1)
//                                            * sizeof (Point2)));
//        for (i = 0; i <= degree; i++) {
//            Vtemp[i] = V[i];
//        }
//
//        /* Triangle computation	*/
//        for (i = 1; i <= degree; i++) {
//            for (j = 0; j <= degree-i; j++) {
//                Vtemp[j].x = (1.0 - t) * Vtemp[j].x + t * Vtemp[j+1].x;
//                Vtemp[j].y = (1.0 - t) * Vtemp[j].y + t * Vtemp[j+1].y;
//            }
//        }
//
//        Q = Vtemp[0];
//        free((void *)Vtemp);
//        return Q;
//    }
//
//
//    /*
//     *  B0, B1, B2, B3 :
//     *	Bezier multipliers
//     */
//    static double B0(double u)
////    double	u;
//    {
//        double tmp = 1.0 - u;
//        return (tmp * tmp * tmp);
//    }
//
//
//    static double B1(double u)
////    double	u;
//    {
//        double tmp = 1.0 - u;
//        return (3 * u * (tmp * tmp));
//    }
//
//    static double B2(double u)
////    double	u;
//    {
//        double tmp = 1.0 - u;
//        return (3 * u * u * tmp);
//    }
//
//    static double B3(double u)
////    double	u;
//    {
//        return (u * u * u);
//    }
//
//
//
//    /*
//     * ComputeLeftTangent, ComputeRightTangent, ComputeCenterTangent :
//     *Approximate unit tangents at endpoints and "center" of digitized curve
//     */
//    static Vector2 ComputeLeftTangent(Point2 *d, int end)
////    Point2	*d;			/*  Digitized points*/
////    int		end;		/*  Index to "left" end of region */
//    {
//        Vector2	tHat1;
//        tHat1 = V2SubII(d[end+1], d[end]);
//        tHat1 = *V2Normalize(&tHat1);
//        return tHat1;
//    }
//
//    static Vector2 ComputeRightTangent(Point2* d, int end)
////    Point2	*d;			/*  Digitized points		*/
////    int		end;		/*  Index to "right" end of region */
//    {
//        Vector2	tHat2;
//        tHat2 = V2SubII(d[end-1], d[end]);
//        tHat2 = *V2Normalize(&tHat2);
//        return tHat2;
//    }
//    
//    
//    static Vector2 ComputeCenterTangent(Point2* d, int center)
////    Point2	*d;			/*  Digitized points			*/
////    int		center;		/*  Index to point inside region	*/
//    {
//        Vector2	V1, V2, tHatCenter;
//        
//        V1 = V2SubII(d[center-1], d[center]);
//        V2 = V2SubII(d[center], d[center+1]);
//        tHatCenter.x = (V1.x + V2.x)/2.0;
//        tHatCenter.y = (V1.y + V2.y)/2.0;
//        tHatCenter = *V2Normalize(&tHatCenter);
//        return tHatCenter;
//    }
//    
//    
//    /*
//     *  ChordLengthParameterize :
//     *	Assign parameter values to digitized points
//     *	using relative distances between points.
//     */
//
////    Point2	*d;			/* Array of digitized points */
////    int		first, last;		/*  Indices defining region	*/
//    {
//        int		i;
//        double	*u;			/*  Parameterization		*/
//        
//        u = (double *)malloc((unsigned)(last-first+1) * sizeof(double));
//        
//        u[0] = 0.0;
//        for (i = first+1; i <= last; i++) {
//            u[i-first] = u[i-first-1] +
//            V2DistanceBetween2Points(&d[i], &d[i-1]);
//        }
//        
//        for (i = first + 1; i <= last; i++) {
//            u[i-first] = u[i-first] / u[last-first];
//        }
//        
//        return(u);
//    }
//    
//    
//    
//    
//    /*
//     *  ComputeMaxError :
//     *	Find the maximum squared distance of digitized points
//     *	to fitted curve.
//     */
//    static double ComputeMaxError(Point2* d, int first, int last, BezierCure bezCurve, double* u, int* splitPoint)
////    Point2	*d;			/*  Array of digitized points	*/
////    int		first, last;		/*  Indices defining region	*/
////    BezierCurve	bezCurve;		/*  Fitted Bezier curve		*/
////    double	*u;			/*  Parameterization of points	*/
////    int		*splitPoint;		/*  Point of maximum error	*/
//    {
//        int		i;
//        double	maxDist;		/*  Maximum error		*/
//        double	dist;		/*  Current error		*/
//        Point2	P;			/*  Point on curve		*/
//        Vector2	v;			/*  Vector from point to curve	*/
//        
//        *splitPoint = (last - first + 1)/2;
//        maxDist = 0.0;
//        for (i = first + 1; i < last; i++) {
//            P = BezierII(3, bezCurve, u[i-first]);
//            v = V2SubII(P, d[i]);
//            dist = V2SquaredLength(&v);
//            if (dist >= maxDist) {
//                maxDist = dist;
//                *splitPoint = i;
//            }
//        }
//        return (maxDist);
//    }
//
//    static Vector2 V2AddII(Vector2 a, Vector2 b)
////    Vector2 a, b;
//    {
//        Vector2	c;
//        c.x = a.x + b.x;  c.y = a.y + b.y;
//        return (c);
//    }
//    static Vector2 V2ScaleIII(Vector2 v, double s)
////    Vector2	v;
////    double	s;
//    {
//        Vector2 result;
//        result.x = v.x * s;
//        result.y = v.y * s;
//        return (result);
//    }
//    
//    static Vector2 V2SubII(Vector2 a, Vector2 b)
////    Vector2	a, b;
//    {
//        Vector2	c;
//        c.x = a.x - b.x; c.y = a.y - b.y;
//        return (c);
//    }
}