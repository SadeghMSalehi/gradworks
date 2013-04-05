/***************************************************************************
 fitpack.cpp  -  description
 -------------------
 begin                : Wed Feb 27 2002
 copyright            : (C) 2011 by Werner Stille
 email                : stille@uni-freiburg.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <QtGlobal>
#include "piqFitPack.h"

static void fpback(double** a, const double* z, int n, int k, double* c);
static void fpbisp(const double* tx, int nx, const double* ty, int ny,
                   const double* c, int kx, int ky, const double* x, int mx,
                   const double* y, int my, double* z);
static void fpbspl(const double* t, int k, double x, int l, double* h);
static void fpchec(const double* x, int m, const double* t, int n, int k,
                   int* ier);
static void fpcurf(int iopt, const double* x, const double* y,
                   const double* w, int m, double xb, double xe, int k,
                   double s, int nest, double tol, int maxit, int k1, int k2,
                   int* n, double* t, double* c, double* fp, double* fpint,
                   double* z, double** a, double** b, double** g, double** q,
                   int* nrdata, int* ier);
static void fpcuro(double a, double b, double c, double d, double* x, int* n);
static void fpdisc(const double* t, int n, int k2, double** b);
static void fpgivs(double piv, double *ww, double *v_cos, double *v_sin);
static void fpknot(const double* x, double* t, int* n, double* fpint,
                   int* nrdata, int* nrint, int istart);
static void fpintb(const double* t, int n, double* bint, int nk1, double x,
                   double y);
static void fporde(const double* x, const double* y, int m, int kx, int ky,
                   const double* tx, int nx, const double* ty, int ny,
                   int* nummer, int* index, int nreg);
static void fprank(double** a, double *f, int n, int m, double tol,
                   double *c, double* sq, int* rank);
static double fprati(double* p1, double* f1, double p2, double f2,
                     double* p3, double* f3);
static void fprota(double v_cos, double v_sin, double* a, double* b);
static void fpsurf(int iopt, int m, double* x, double* y, const double* z,
                   const double* w, double xb, double xe, double yb,
                   double ye, int kxx, int kyy, double s, int nxest,
                   int nyest, double eta, double tol, int maxit,
                   int nrest, int* nx0, double* tx, int* ny0, double* ty,
                   double* c, double* fp, double* fp0, double* fpint,
                   double* coord, double* f, double* ff, double** a,
                   double** q, double** bx, double** by, double** spx,
                   double** spy, double* h, int* ier);

void fpback(double** a, const double* z, int n, int k, double* c)
{
    /*  subroutine fpback calculates the solution of the system of
     equations a*c = z with a a n x n upper triangular matrix
     of bandwidth k. */
    int i, i1, j, k1, l, m;
    double store;
    k1 = k - 1;
    i = n - 1;
    c[i] = z[i] / a[i][0];
    if (i) {
        for (j = 2; j <= n; j++) {
            store = z[i - 1];
            i1 = (j <= k1) ? (j - 1) : k1;
            m = i;
            for (l = 1; l <= i1; l++)
                store -= c[++m - 1] * a[i - 1][l];
            c[i - 1] = store / a[i - 1][0];
            i--;
        }
    }
}

void fpbisp(const double* tx, int nx, const double* ty, int ny,
            const double* c, int kx, int ky, const double* x, int mx,
            const double* y, int my, double* z)
{
    int i, i1, j, j1, kx1, ky1, l, l1, l2, m, nkx1, nky1;
    double arg, sp, tb, te;
    double **wx, **wy;
    double h[6];
    int* lx = new int[mx];
    int* ly = new int[my];
    wx = new double*[mx];
    kx1 = kx + 1;
    wx[0] = new double[mx * kx1];
    for (i = 1; i < mx; i++)
        wx[i] = &wx[0][i * kx1];
    nkx1 = nx - kx1;
    tb = tx[kx1 - 1];
    te = tx[nkx1];
    l = kx1;
    for (i = 0; i < mx; ++i) {
        arg = x[i];
        if (arg < tb)
            arg = tb;
        if (arg > te)
            arg = te;
        while (!(arg < tx[l] || l == nkx1))
            l++;
        fpbspl(tx, kx, arg, l, h);
        lx[i] = l - kx1;
        for (j = 0; j < kx1; ++j)
            wx[i][j] = h[j];
    }
    ky1 = ky + 1;
    wy = new double*[my];
    wy[0] = new double[my * ky1];
    for (i = 1; i < my; i++)
        wy[i] = &wy[0][i * ky1];
    nky1 = ny - ky1;
    tb = ty[ky1 - 1];
    te = ty[nky1];
    l = ky1;
    for (i = 0; i < my; ++i) {
        arg = y[i];
        if (arg < tb)
            arg = tb;
        if (arg > te)
            arg = te;
        while (!(arg < ty[l] || l == nky1))
            l++;
        fpbspl(ty, ky, arg, l, h);
        ly[i] = l - ky1;
        for (j = 0; j < ky1; ++j)
            wy[i][j] = h[j];
    }
    m = 0;
    for (i = 0; i < mx; i++) {
        l = lx[i] * nky1;
        for (i1 = 0; i1 < kx1; i1++)
            h[i1] = wx[i][i1];
        for (j = 0; j < my; j++) {
            l1 = l + ly[j];
            sp = 0.0;
            for (i1 = 0; i1 < kx1; i1++) {
                l2 = l1;
                for (j1 = 0; j1 < ky1; j1++)
                    sp += c[l2++] * h[i1] * wy[j][j1];
                l1 += nky1;
            }
            z[m++] = sp;
        }
    }
    delete [] wy[0];
    delete [] wy;
    delete [] wx[0];
    delete [] wx;
    delete [] ly;
    delete [] lx;
    return;
}

void fpbspl(const double* t, int k, double x, int l, double* h)
{
    /*  subroutine fpbspl evaluates the (k+1) non-zero b-splines of
     degree k at t(l) <= x < t(l+1) using the stable recurrence
     relation of de boor and cox. */
    int i, j, li, lj;
    double f;
    double hh[5];
    h[0] = 1;
    for (j = 1; j <= k; j++) {
        for (i = 0; i < k; i++)
            hh[i] = h[i];
        h[0] = 0;
        for (i = 0; i < j; i++) {
            li = l + i;
            lj = li - j;
            f = hh[i] / (t[li] - t[lj]);
            h[i] += f * (t[li] - x);
            h[i + 1] = f * (x - t[lj]);
        }
    }
}

void fpchec(const double* x, int m, const double* t, int n, int k, int* ier)
{
    /*  subroutine fpchec verifies the number and the position of the knots
     t(j),j=1,2,...,n of a spline of degree k, in relation to the number
     and the position of the data points x(i),i=1,2,...,m. if all of the
     following conditions are fulfilled, the error parameter ier is set
     to zero. if one of the conditions is violated ier is set to ten.
     1) k+1 <= n-k-1 <= m
     2) t(1) <= t(2) <= ... <= t(k+1)
     t(n-k) <= t(n-k+1) <= ... <= t(n)
     3) t(k+1) < t(k+2) < ... < t(n-k)
     4) t(k+1) <= x(i) <= t(n-k)
     5) the conditions specified by schoenberg and whitney must hold
     for at least one subset of data points, i.e. there must be a
     subset of data points y(j) such that
     t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1 */
    int i, j, k1, k2, l, nk1, nk2, nk3;
    k1 = k + 1;
    k2 = k1 + 1;
    nk1 = n - k1;
    nk2 = nk1 + 1;
    *ier = 10;
    /*  check condition no 1 */
    if (nk1 < k1 || nk1 > m)
        return;
    /*  check condition no 2 */
    j = n - 1;
    for (i = 0; i < k; i++) {
        if (t[i] > t[i + 1])
            return;
        if (t[j] < t[j - 1])
            return;
        j--;
    }
    /*  check condition no 3 */
    for (i = k1; i < nk2; ++i)
        if (t[i] <= t[i - 1])
            return;
    /*  check condition no 4 */
    if (x[0] < t[k] || x[m - 1] > t[nk1])
        return;
    /*  check condition no 5 */
    if (x[0] >= t[k1] || x[m - 1] <= t[nk1 - 1])
        return;
    i = 1;
    l = k2;
    nk3 = nk1 - 1;
    if (nk3 >= 2) {
        for (j = 1; j < nk3; j++) {
            double tj = t[j];
            double tl = t[l++];
            do {
                i++;
                if (i >= m)
                    return;
            } while (x[i - 1] <= tj);
            if (x[i - 1] >= tl)
                return;
        }
    }
    *ier = 0;
    return;
}

void fpcurf(int iopt, const double* x, const double* y, const double* w, int m,
            double xb, double xe, int k, double s, int nest, double tol,
            int maxit, int k1, int k2, int* n, double* t, double* c,
            double* fp, double* fpint, double* z, double** a, double** b,
            double** g, double** q, int* nrdata, int* ier)
{
    int i, ich1, ich3, it, iter, i1, i2, j, k3, l, l0, mk1,v_new,
    nmin, npl1, nrint, n8;
    int nk1 = 0;
    int nmax = 0;
    int nplus = 0;
    double v_cos, d1, fpart, f1, f2, f3, p, pinv, piv, p1, p2, p3, rn,
    v_sin, store, term, wi, xi, yi;
    double acc = 0;
    double fpms = 0;
    double fpold = 0;
    double fp0 = 0;
    double h[7];
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     part 1: determination of the number of knots and their position     c
     **************************************************************      c
     given a set of knots we compute the least-squares spline sinf(x),   c
     and the corresponding sum of squared residuals fp=f(p=inf).         c
     if iopt=-1 sinf(x) is the requested approximation.                  c
     if iopt=0 or iopt=1 we check whether we can accept the knots:       c
     if fp <=s we will continue with the current set of knots.         c
     if fp > s we will increase the number of knots and compute the    c
     corresponding least-squares spline until finally fp<=s.        c
     the initial choice of knots depends on the value of s and iopt.   c
     if s=0 we have spline interpolation; in that case the number of   c
     knots equals nmax = m+k+1.                                        c
     if s > 0 and                                                      c
     iopt=0 we first compute the least-squares polynomial of         c
     degree k; n = nmin = 2*k+2                                      c
     iopt=1 we start with the set of knots found at the last         c
     call of the routine, except for the case that s > fp0; then     c
     we compute directly the least-squares polynomial of degree k.   c
     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /*  determine nmin, the number of knots for polynomial approximation. */
    nmin = k1 << 1;
    if (iopt < 0)
        goto L60;
    /*  calculation of acc, the absolute tolerance for the root of f(p)=s. */
    acc = tol * s;
    /*  determine nmax, the number of knots for spline interpolation. */
    nmax = m + k1;
    if (s <= 0) {
        /*  if s=0, s(x) is an interpolating spline.
         test whether the required storage space exceeds the available one. */
        *n = nmax;
        if (nmax > nest)
            goto L420;
        /*  find the position of the interior knots in case of interpolation. */
    L10:
        mk1 = m - k1;
        if (mk1 == 0)
            goto L60;
        k3 = k / 2;
        i = k2;
        j = k3 + 2;
        if ((k3 << 1) != k) {
            for (l = 0; l < mk1; ++l) {
                t[i - 1] = x[j - 1];
                ++i;
                ++j;
            }
            goto L60;
        }
        for (l = 0; l < mk1; ++l) {
            t[i - 1] = (x[j - 1] + x[j - 2]) * 0.5;
            ++i;
            ++j;
        }
        goto L60;
    }
    /*  if s>0 our initial choice of knots depends on the value of iopt.
     if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
     polynomial of degree k which is a spline without interior knots.
     if iopt=1 and fp0>s we start computing the least squares spline
     according to the set of knots found at the last call of the routine. */
    if (iopt && (*n != nmin)) {
        fp0 = fpint[*n - 1];
        fpold = fpint[*n - 2];
        nplus = nrdata[*n - 1];
        if (fp0 > s)
            goto L60;
    }
    *n = nmin;
    fpold = 0;
    nplus = 0;
    nrdata[0] = m - 2;
    /*  main loop for the different sets of knots. m is a save upper bound
     for the number of trials. */
L60:
    for (iter = 0; iter < m; ++iter) {
        if (*n == nmin)
            *ier = -2;
        /*  find nrint, tne number of knot intervals. */
        nrint = *n - nmin + 1;
        /*  find the position of the additional knots which are needed for
         the b-spline representation of s(x). */
        nk1 = *n - k1;
        i = *n;
        for (j = 0; j < k1; ++j) {
            t[j] = xb;
            t[--i] = xe;
        }
        /*  compute the b-spline coefficients of the least-squares spline
         sinf(x). the observation matrix a is built up row by row and
         reduced to upper triangular form by givens transformations.
         at the same time fp=f(p=inf) is computed. */
        *fp = 0;
        /*  initialize the observation matrix a. */
        for (i = 0; i < nk1; ++i) {
            z[i] = 0;
            for (j = 0; j < k1; ++j)
                a[i][j] = 0;
        }
        l = k1;
        for (it = 0; it < m; ++it) {
            /*  fetch the current data point x(it),y(it). */
            xi = x[it];
            wi = w[it];
            yi = y[it] * wi;
            /*  search for knot interval t(l) <= xi < t(l+1). */
            while (!(xi < t[l] || l == nk1))
                ++l;
            /*  evaluate the (k+1) non-zero b-splines at xi and store them in q. */
            fpbspl(t, k, xi, l, h);
            for (i = 0; i < k1; ++i) {
                q[it][i] = h[i];
                h[i] *= wi;
            }
            /*  rotate the new row of the observation matrix into triangle. */
            j = l - k1;
            for (i = 1; i <= k1; ++i) {
                ++j;
                piv = h[i - 1];
                if (piv) {
                    /*  calculate the parameters of the givens transformation. */
                    fpgivs(piv, a[j - 1], &v_cos, &v_sin);
                    /*  transformations to right hand side. */
                    fprota(v_cos, v_sin, &yi, &z[j - 1]);
                    if (i == k1)
                        break;
                    i2 = 1;
                    for (i1 = i; i1 < k1; ++i1)
                    /*  transformations to left hand side. */
                        fprota(v_cos, v_sin, &h[i1], &a[j - 1][i2++]);
                }
            }
            /*  add contribution of this row to the sum of squares of residual
             right hand sides. */
            *fp += yi * yi;
        }
        if (*ier == -2)
            fp0 = *fp;
        fpint[*n - 1] = fp0;
        fpint[*n - 2] = fpold;
        nrdata[*n - 1] = nplus;
        /*  backward substitution to obtain the b-spline coefficients. */
        fpback(a, z, nk1, k1, c);
        /*  test whether the approximation sinf(x) is an acceptable solution. */
        if (iopt < 0)
            return;
        fpms = *fp - s;
        if (fabs(fpms) < acc)
            return;
        /*  if f(p=inf) < s accept the choice of knots. */
        if (fpms < 0)
            goto L250;
        /*  if n = nmax, sinf(x) is an interpolating spline. */
        if (*n == nmax)
            goto L430;
        /*  increase the number of knots.
         if n=nest we cannot increase the number of knots because of
         the storage capacity limitation. */
        if (*n == nest)
            goto L420;
        /*  determine the number of knots nplus we are going to add. */
        if (*ier == 0) {
            npl1 = nplus << 1;
            rn = (double) nplus;
            if ((fpold - *fp) > acc)
                npl1 = (int) (rn * fpms / (fpold - *fp));
            nplus = qMin(nplus * 2, qMax(qMax(npl1, nplus / 2), 1));
        } else {
            nplus = 1;
            *ier = 0;
        }
        fpold = *fp;
        /*  compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
         t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint. */
        fpart = 0;
        i = 0;
        l = k2;
        v_new = 0;
        for (it = 1; it <= m; ++it) {
            if (!(x[it - 1] < t[l - 1] || l > nk1)) {
                v_new = 1;
                ++l;
            }
            term = 0;
            l0 = l - k2;
            for (j = 0; j < k1; ++j)
                term += c[l0++] * q[it - 1][j];
            d1 = w[it - 1] * (term - y[it - 1]);
            term = d1 * d1;
            fpart += term;
            if (v_new == 0)
                continue;
            store = term * 0.5;
            fpint[i++] = fpart - store;
            fpart = store;
            v_new = 0;
        }
        fpint[nrint - 1] = fpart;
        for (l = 0; l < nplus; ++l) {
            /*  add a new knot. */
            fpknot(x, t, n, fpint, nrdata, &nrint, 1);
            /*  if n=nmax we locate the knots as for interpolation. */
            if (*n == nmax)
                goto L10;
            /*  test whether we cannot further increase the number of knots. */
            if (*n == nest)
                break;
        }
        /*  restart the computations with the new set of knots. */
    }
    /*  test whether the least-squares kth degree polynomial is a solution
     of our approximation problem. */
L250:
    if (*ier == -2)
        return;
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     part 2: determination of the smoothing spline sp(x).                c
     ***************************************************                 c
     we have determined the number of knots and their position.          c
     we now compute the b-spline coefficients of the smoothing spline    c
     sp(x). the observation matrix a is extended by the rows of matrix   c
     b expressing that the kth derivative discontinuities of sp(x) at    c
     the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     c
     ponding weights of these additional rows are set to 1/p.            c
     iteratively we then have to determine the value of p such that      c
     f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that    c
     the least-squares kth degree polynomial corresponds to p=0, and     c
     that the least-squares spline corresponds to p=infinity. the        c
     iteration process which is proposed here, makes use of rational     c
     interpolation. since f(p) is a convex and strictly decreasing       c
     function of p, it can be approximated by a rational function        c
     r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  c
     ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      c
     to calculate the new value of p such that r(p)=s. convergence is    c
     guaranteed by taking f1>0 and f3<0.                                 c
     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /*  evaluate the discontinuity jump of the kth derivative of the
     b-splines at the knots t(l),l=k+2,...n-k-1 and store in b. */
    fpdisc(t, *n, k2, b);
    /*  initial value for p. */
    p1 = 0;
    f1 = fp0 - s;
    p3 = -1;
    f3 = fpms;
    p = 0;
    for (i = 0; i < nk1; ++i)
        p += a[i][0];
    rn = (double) nk1;
    p = rn / p;
    ich1 = 0;
    ich3 = 0;
    n8 = *n - nmin;
    /*  iteration process to find the root of f(p) = s. */
    for (iter = 1; iter <= maxit; ++iter) {
        /*  the rows of matrix b with weight 1/p are rotated into the
         triangularised observation matrix a which is stored in g. */
        pinv = 1 / p;
        for (i = 0; i < nk1; ++i) {
            c[i] = z[i];
            g[i][k2 - 1] = 0;
            for (j = 0; j < k1; ++j)
                g[i][j] = a[i][j];
        }
        for (it = 1; it <= n8; ++it) {
            /*  the row of matrix b is rotated into triangle by givens transformation */
            for (i = 0; i < k2; ++i)
                h[i] = b[it - 1][i] * pinv;
            yi = 0;
            for (j = it; j <= nk1; ++j) {
                piv = h[0];
                /*  calculate the parameters of the givens transformation. */
                fpgivs(piv, g[j - 1], &v_cos, &v_sin);
                /*  transformations to right hand side. */
                fprota(v_cos, v_sin, &yi, &c[j - 1]);
                if (j == nk1)
                    goto L300;
                i2 = k1;
                if (j > n8)
                    i2 = nk1 - j;
                for (i = 1; i <= i2; ++i) {
                    /*  transformations to left hand side. */
                    fprota(v_cos, v_sin, &h[i], &g[j - 1][i]);
                    h[i - 1] = h[i];
                }
                h[i2] = 0;
            }
        L300:
            ;
        }
        /*  backward substitution to obtain the b-spline coefficients. */
        fpback(g, c, nk1, k2, c);
        /*  computation of f(p). */
        *fp = 0;
        l = k2;
        for (it = 0; it < m; ++it) {
            while (!(x[it] < t[l - 1] || l > nk1))
                ++l;
            l0 = l - k2;
            term = 0;
            for (j = 0; j < k1; ++j)
                term += c[l0++] * q[it][j];
            d1 = w[it] * (term - y[it]);
            *fp += d1 * d1;
        }
        /*  test whether the approximation sp(x) is an acceptable solution. */
        fpms = *fp - s;
        if (fabs(fpms) < acc)
            return;
        /*  test whether the maximal number of iterations is reached. */
        if (iter == maxit)
            goto L400;
        /*  carry out one more step of the iteration process. */
        p2 = p;
        f2 = fpms;
        if (ich3 == 0) {
            if ((f2 - f3) <= acc) {
                /*  our initial choice of p is too large. */
                p3 = p2;
                f3 = f2;
                p *= 0.04;
                if (p <= p1)
                    p = p1 * 0.9 + p2 * 0.1;
                goto L360;
            }
            if (f2 < 0)
                ich3 = 1;
        }
        if (!ich1) {
            if ((f1 - f2) <= acc) {
                /*  our initial choice of p is too small */
                p1 = p2;
                f1 = f2;
                p /= 0.04;
                if (p3 < 0)
                    goto L360;
                if (p >= p3)
                    p = p2 * 0.1 + p3 * 0.9;
                goto L360;
            }
            if (f2 > 0)
                ich1 = 1;
        }
        /*  test whether the iteration process proceeds as theoretically
         expected. */
        if ((f2 >= f1) || (f2 <= f3))
            goto L410;
        /*  find the new value for p. */
        p = fprati(&p1, &f1, p2, f2, &p3, &f3);
    L360:
        ;
    }
    /*  error codes and messages. */
L400:
    *ier = 3;
    return;
L410:
    *ier = 2;
    return;
L420:
    *ier = 1;
    return;
L430:
    *ier = -1;
    return;
}

void fpcuro(double a, double b, double c, double d, double* x, int* n)
{
    /*  subroutine fpcuro finds the real zeros of a cubic polynomial
     p(x) = a*x**3+b*x**2+c*x+d.

     calling sequence:
     call fpcuro(a,b,c,d,x,n)

     input parameters:
     a,b,c,d: real values, containing the coefficients of p(x).

     output parameters:
     x      : real array,length 3, which contains the real zeros of p(x)
     n      : integer, giving the number of real zeros of p(x). */
    int i;
    double disc, step, f, q, r, u, y, a1, b1, c1, d1, p3, u1, u2, df, r1, r2;
    /*  set constants */
    static const double ovfl = 1e4;
    static const double e3 = 1.0 / 3;
    static const double pi3 = M_PI / 3;
    a1 = fabs(a);
    b1 = fabs(b);
    c1 = fabs(c);
    d1 = fabs(d);
    /*  test whether p(x) is a third degree polynomial. */
    r1 = qMax(b1, c1);
    if (qMax(r1, d1) < a1 * ovfl) {
        /*  p(x) is a third degree polynomial. */
        b1 = b / a * e3;
        c1 = c / a;
        d1 = d / a;
        q = c1 * e3 - b1 * b1;
        r = b1 * b1 * b1 + (d1 - b1 * c1) * 0.5;
        disc = q * q * q + r * r;
        if (disc > 0) {
            u = sqrt(disc);
            u1 = -r + u;
            u2 = -r - u;
            *n = 1;
            r1 = pow(fabs(u1), e3);
            r2 = pow(fabs(u2), e3);
            x[0] = copysign(r1, u1) + copysign(r2, u2) - b1;
        } else {
            u = sqrt((fabs(q)));
            if (r < 0)
                u = -u;
            p3 = atan2(sqrt(-disc), (fabs(r))) * e3;
            u2 = u + u;
            *n = 3;
            x[0] = -u2 * cos(p3) - b1;
            x[1] = u2 * cos(pi3 - p3) - b1;
            x[2] = u2 * cos(pi3 + p3) - b1;
        }
    } else
    /*  test whether p(x) is a second degree polynomial. */
        if (qMax(c1, d1) < b1 * ovfl) {
            /*  p(x) is a second degree polynomial. */
            disc = c * c - 4 * b * d;
            *n = 0;
            if (disc < 0)
                return;
            *n = 2;
            u = sqrt(disc);
            b1 = b + b;
            x[0] = (-c + u) / b1;
            x[1] = (-c - u) / b1;
        } else
        /*  test whether p(x) is a first degree polynomial. */
            if (d1 >= c1 * ovfl) {
                /*  p(x) is a constant function. */
                *n = 0;
                return;
            } else {
                /*  p(x) is a first degree polynomial. */
                *n = 1;
                x[0] = -d / c;
            }
    /*  apply a newton iteration to improve the accuracy of the roots. */
    for (i = 0; i < *n; ++i) {
        y = x[i];
        f = ((a * y + b) * y + c) * y + d;
        df = (3 * a * y + 2 * b) * y + c;
        step = 0;
        if (fabs(f) < fabs(df) * 0.1)
            step = f / df;
        x[i] = y - step;
    }
    return;
}

void fpdisc(const double* t, int n, int k2, double** b)
{
    /*  subroutine fpdisc calculates the discontinuity jumps of the kth
     derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1) */
    int i, ik, j, jk, k, k1, l, lj, lk, lmk, lp, nk1, nrint;
    double an, fac, prod;
    double h[12];
    k1 = k2 - 1;
    k = k1 - 1;
    nk1 = n - k1;
    nrint = nk1 - k;
    an = nrint;
    fac = an / (t[nk1] - t[k1 - 1]);
    for (l = k2; l <= nk1; ++l) {
        lmk = l - k1;
        for (j = 0; j < k1; ++j) {
            ik = j + k1;
            lj = l + j;
            lk = lj - k2;
            h[j] = t[l - 1] - t[lk];
            h[ik] = t[l - 1] - t[lj];
        }
        lp = lmk;
        for (j = 0; j < k2; ++j) {
            jk = j;
            prod = h[j];
            for (i = 1; i <= k; ++i)
                prod *= h[++jk] * fac;
            lk = lp + k1;
            b[lmk - 1][j] = (t[lk - 1] - t[lp++ - 1]) / prod;
        }
    }
}

void fpgivs(double piv, double *ww, double *v_cos, double *v_sin)
{
    /*  subroutine fpgivs calculates the parameters of a givens
     transformation. */
    double dd, store, d1;
    store = fabs(piv);
    d1 = *ww / piv;
    d1 = d1 * d1;
    if (store >= *ww)
        dd = store * sqrt(1 + d1);
    else
        dd = *ww * sqrt(1 + 1 / d1);
    *v_cos = *ww / dd;
    *v_sin = piv / dd;
    *ww = dd;
}

void fpintb(const double* t, int n, double* bint, int nk1, double x, double y)
{
    /*  subroutine fpintb calculates integrals of the normalized b-splines
     nj,k+1(x) of degree k, defined on the set of knots t(j),j=1,2,...n.
     it makes use of the formulae of gaffney for the calculation of
     indefinite integrals of b-splines.

     calling sequence:
     call fpintb(t,n,bint,nk1,x,y)

     input parameters:
     t    : real array,length n, containing the position of the knots.
     n    : integer value, giving the number of knots.
     nk1  : integer value, giving the number of b-splines of degree k,
     defined on the set of knots ,i.e. nk1 = n-k-1.
     x,y  : real values, containing the end points of the integration
     interval.
     output parameter:
     bint : array,length nk1, containing the integrals of the b-splines. */
    int i, j, k, l, j1, k1, ib, li, lj, lk, it, min_v;
    int ia = 0;
    double a, b, d1, f, ak, arg;
    double aint[6], h[6], h1[6];
    k1 = n - nk1;
    ak = (double) k1;
    k = k1 - 1;
    for (i = 0; i < nk1; ++i)
        bint[i] = 0.;
    /*  the integration limits are arranged in increasing order. */
    d1 = x - y;
    if (d1 == 0)
        return;
    if (d1 < 0) {
        a = x;
        b = y;
        min_v = 0;
    } else {
        a = y;
        b = x;
        min_v = 1;
    }
    if (a < t[k])
        a = t[k];
    if (b > t[nk1])
        b = t[nk1];
    /*  using the expression of gaffney for the indefinite integral of a
     b-spline we find that
     bint(j) = (t(j+k+1)-t(j))*(res(j,b)-res(j,a))/(k+1)
     where for t(l) <= x < t(l+1)
     res(j,x) = 0, j=1,2,...,l-k-1
     = 1, j=l+1,l+2,...,nk1
     = aint(j+k-l+1), j=l-k,l-k+1,...,l
     = sumi((x-t(j+i))*nj+i,k+1-i(x)/(t(j+k+1)-t(j+i)))
     i=0,1,...,k */
    l = k1;
    /*  set arg = a. */
    arg = a;
    for (it = 0; it < 2; ++it) {
        /*  search for the knot interval t(l) <= arg < t(l+1). */
        while (!(arg < t[l] || l == nk1))
            l++;
        /*  calculation of aint(j), j=1,2,...,k+1.
         initialization. */
        for (j = 0; j < k1; ++j)
            aint[j] = 0.;
        aint[0] = (arg - t[l - 1]) / (t[l] - t[l - 1]);
        h1[0] = 1.0;
        for (j = 1; j <= k; ++j) {
            /*  evaluation of the non-zero b-splines of degree j at arg,i.e.
             h(i+1) = nl-j+i,j(arg), i=0,1,...,j. */
            h[0] = 0.0;
            for (i = 0; i < j; ++i) {
                li = l + i;
                lj = li - j;
                f = h1[i] / (t[li] - t[lj]);
                h[i] += f * (t[li] - arg);
                h[i + 1] = f * (arg - t[lj]);
            }
            /*  updating of the integrals aint. */
            j1 = j + 1;
            for (i = 0; i < j1; ++i) {
                li = l + i;
                lj = li - j1;
                aint[i] += h[i] * (arg - t[lj]) / (t[li] - t[lj]);
                h1[i] = h[i];
            }
        }
        if (it == 1)
            break;
        /*  updating of the integrals bint */
        lk = l - k;
        ia = lk;
        --lk;
        for (i = 0; i < k1; ++i)
            bint[lk++] = -aint[i];
        /*  set arg = b. */
        arg = b;
    }
    /*  updating of the integrals bint. */
    lk = l - k;
    ib = lk - 1;
    --lk;
    for (i = 0; i < k1; ++i)
        bint[lk++] += aint[i];
    if (ib >= ia)
        for (i = ia - 1; i < ib; ++i)
            bint[i] += 1;
    /*  the scaling factors are taken into account. */
    f = 1.0 / ak;
    for (i = 0; i < nk1; ++i)
        bint[i] *= (t[i + k1] - t[i]) * f;
    /*  the order of the integration limits is taken into account. */
    if (min_v == 0)
        return;
    for (i = 0; i < nk1; ++i)
        bint[i] = -bint[i];
    return;
}

void fpknot(const double* x, double* t, int* n, double* fpint, int* nrdata,
            int* nrint, int istart)
{
    /*  subroutine fpknot locates an additional knot for a spline of degree
     k and adjusts the corresponding parameters,i.e.
     t     : the position of the knots.
     n     : the number of knots.
     nrint : the number of knotintervals.
     fpint : the sum of squares of residual right hand sides
     for each knot interval.
     nrdata: the number of data points inside each knot interval.
     istart indicates that the smallest data point at which the new knot
     may be added is x(istart+1) */
    int ihalf, j, jbegin, jj, jk, jpoint, k, next, nrx;
    int maxbeg = 0;
    int maxpt = 0;
    int number = 0;
    double an, am, fpmax;
    k = (*n - *nrint - 1) / 2;
    /*  search for knot interval t(number+k) <= x <= t(number+k+1) where
     fpint(number) is maximal on the condition that nrdata(number)
     not equals zero. */
    fpmax = 0;
    jbegin = istart;
    for (j = 0; j < *nrint; ++j) {
        jpoint = nrdata[j];
        if ((fpmax < fpint[j]) && jpoint) {
            fpmax = fpint[j];
            number = j + 1;
            maxpt = jpoint;
            maxbeg = jbegin;
        }
        jbegin = jbegin + jpoint + 1;
    }
    /*  let coincide the new knot t(number+k+1) with a data point x(nrx)
     inside the old knot interval t(number+k) <= x <= t(number+k+1). */
    ihalf = maxpt / 2 + 1;
    nrx = maxbeg + ihalf;
    next = number + 1;
    if (next <= *nrint)
    /*  adjust the different parameters. */
        for (j = next; j <= *nrint; ++j) {
            jj = next + *nrint - j;
            fpint[jj] = fpint[jj - 1];
            nrdata[jj] = nrdata[jj - 1];
            jk = jj + k;
            t[jk] = t[jk - 1];
        }
    nrdata[number - 1] = ihalf - 1;
    nrdata[next - 1] = maxpt - ihalf;
    am = (double) maxpt;
    an = (double) nrdata[number - 1];
    fpint[number - 1] = fpmax * an / am;
    an = (double) nrdata[next - 1];
    fpint[next - 1] = fpmax * an / am;
    jk = next + k;
    t[jk - 1] = x[nrx - 1];
    ++(*n);
    ++(*nrint);
}

void fporde(const double* x, const double* y, int m, int kx, int ky,
            const double* tx, int nx, const double* ty, int ny, int* nummer,
            int* index, int nreg)
{
    /*  subroutine fporde sorts the data points (x(i),y(i)),i=1,2,...,m
     according to the panel tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong
     to. for each panel a stack is constructed  containing the numbers
     of data points lying inside; index(j),j=1,2,...,nreg points to the
     first data point in the jth panel while nummer(i),i=1,2,...,m gives
     the number of the next data point in the panel. */
    int i, im, k, kx1, ky1, k1, l, l1, nk1x, nk1y, num, nyy;
    double xi, yi;
    kx1 = kx + 1;
    ky1 = ky + 1;
    nk1x = nx - kx1;
    nk1y = ny - ky1;
    nyy = nk1y - ky;
    for (i = 0; i < nreg; i++)
        index[i] = 0;
    for (im = 0; im < m; im++) {
        xi = x[im];
        yi = y[im];
        l = kx1;
        l1 = l + 1;
        while (!(xi < tx[l] || l == nk1x)) {
            l = l1;
            l1 = l + 1;
        }
        k = ky1;
        k1 = k + 1;
        while (!(yi < ty[k] || k == nk1y)) {
            k = k1;
            k1 = k + 1;
        }
        num = (l - kx1) * nyy + k - ky - 1;
        nummer[im] = index[num];
        index[num] = im + 1;
    }
}

void fprank(double** a, double *f, int n, int m, double tol, double *c,
            double* sq, int* rank)
{
    /*  subroutine fprank finds the minimum norm solution of a least-
     squares problem in case of rank deficiency.

     input parameters:
     a : array, which contains the non-zero elements of the observation
     matrix after triangularization by givens transformations.
     f : array, which contains the transformed right hand side.
     n : integer,wich contains the dimension of a.
     m : integer, which denotes the bandwidth of a.
     tol : real value, giving a threshold to determine the rank of a.

     output parameters:
     c : array, which contains the minimum norm solution.
     sq : real value, giving the contribution of reducing the rank
     to the sum of squared residuals.
     rank : integer, which contains the rank of matrix a. */

    int i, ii, ij, i1, i2, j, jj, j1, j2, j3 = 0, k, kk, m1, nl;
    double v_cos, fac, piv, v_sin, yi, store, stor1, stor2, stor3;
    double** aa = new double*[n];
    double* ff = new double[n];
    double* h = new double[m];
    aa[0] = new double[n * m];
    for (i = 1; i < n; i++)
        aa[i] = aa[i - 1] + m;
    m1 = m - 1;
    /*  the rank deficiency nl is considered to be the number of sufficient
     small diagonal elements of a. */
    nl = 0;
    *sq = 0;
    for (i = 1; i <= n; i++) {

        if (a[i - 1][0] > tol)
            goto L90;
        /*  if a sufficient small diagonal element is found, we put it to
         zero. the remainder of the row corresponding to that zero diagonal
         element is then rotated into triangle by givens rotations.
         the rank deficiency is increased by one. */
        nl++;
        if (i == n)
            goto L90;
        yi = f[i - 1];
        for (j = 1; j <= m1; j++)
            h[j - 1] = a[i - 1][j];
        h[m - 1] = 0;
        i1 = i + 1;
        for (ii = i1; ii <= n; ii++) {
            i2 = qMin(n - ii, m1);
            piv = h[0];
            if (piv == 0)
                goto L30;
            fpgivs(piv, &a[ii - 1][0], &v_cos, &v_sin);
            fprota(v_cos, v_sin, &yi, &f[ii - 1]);
            if (i2 == 0)
                goto L70;
            for (j = 1; j <= i2; j++) {
                j1 = j + 1;
                fprota(v_cos, v_sin, &h[j1 - 1], &a[ii - 1][j1 - 1]);
                h[j - 1] = h[j1 - 1];
            }
            goto L50;
        L30:
            if (i2 == 0)
                goto L70;
            for (j = 0; j < i2; j++)
                h[j] = h[j + 1];
        L50:
            h[i2] = 0;
        }
        /*  add to the sum of squared residuals the contribution of deleting
         the row with small diagonal element. */
    L70:
        *sq += yi * yi;
    L90:
        ;
    }
    /*  rank denotes the rank of a. */
    *rank = n - nl;
    /*  let b denote the (rank*n) upper trapezoidal matrix which can be
     obtained from the (n*n) upper triangular matrix a by deleting
     the rows and interchanging the columns corresponding to a zero
     diagonal element. if this matrix is factorized using givens
     transformations as  b = (r) (u)  where
     r is a (rank*rank) upper triangular matrix,
     u is a (rank*n) orthonormal matrix
     then the minimal least-squares solution c is given by c = b' v,
     where v is the solution of the system  (r) (r)' v = g  and
     g denotes the vector obtained from the old right hand side f, by
     removing the elements corresponding to a zero diagonal element of a.
     initialization. */
    for (i = 1; i <= *rank; i++)
        for (j = 1; j <= m; j++)
            aa[i - 1][j - 1] = 0;
    /*  form in aa the upper triangular matrix obtained from a by
     removing rows and columns with zero diagonal elements. form in ff
     the new right hand side by removing the elements of the old right
     hand side corresponding to a deleted row. */
    ii = 0;
    for (i = 1; i <= n; i++) {
        if (a[i - 1][0] <= tol)
            goto L120;
        ii++;
        ff[ii - 1] = f[i - 1];
        aa[ii - 1][0] = a[i - 1][0];
        jj = ii;
        kk = 1;
        j = i;
        j1 = qMin(j - 1, m1);
        if (j1 == 0)
            goto L120;
        for (k = 1; k <= j1; k++) {
            j--;
            if (a[j - 1][0] <= tol)
                goto L110;
            kk++;
            jj--;
            aa[jj - 1][kk - 1] = a[j - 1][k];
        L110:
            ;
        }
    L120:
        ;
    }
    /*  form successively in h the columns of a with a zero diagonal element. */
    ii = 0;
    for (i = 1; i <= n; i++) {
        ii++;
        if (a[i - 1][0] > tol)
            goto L200;
        ii--;
        if (ii == 0)
            goto L200;
        jj = 1;
        j = i;
        j1 = qMin(j - 1, m1);
        for (k = 1; k <= j1; k++) {
            j--;
            if (a[j - 1][0] <= tol)
                goto L130;
            h[jj - 1] = a[j - 1][k];
            jj++;
        L130:
            ;
        }
        for (kk = jj; kk <= m; kk++)
            h[kk - 1] = 0;
        /*  rotate this column into aa by givens transformations. */
        jj = ii;
        for (i1 = 1; i1 <= ii; i1++) {
            j1 = qMin(jj - 1, m1);
            piv = h[0];
            if (piv != 0)
                goto L160;
            if (j1 == 0)
                goto L200;
            for (j2 = 1; j2 <= j1; j2++) {
                j3 = j2 + 1;
                h[j2 - 1] = h[j3 - 1];
            }
            goto L180;
        L160:
            fpgivs(piv, &aa[jj - 1][0], &v_cos, &v_sin);
            if (j1 == 0)
                goto L200;
            kk = jj;
            for (j2 = 1; j2 <= j1; j2++) {
                j3 = j2 + 1;
                kk--;
                fprota(v_cos, v_sin, &h[j3 - 1], &aa[kk - 1][j3 - 1]);
                h[j2 - 1] = h[j3 - 1];
            }
        L180:
            jj--;
            h[j3 - 1] = 0;
        }
    L200:
        ;
    }
    /*  solve the system (aa) (f1) = ff */
    ff[*rank - 1] /= aa[*rank - 1][0];
    i = *rank - 1;
    if (i == 0)
        goto L230;
    for (j = 2; j <= *rank; j++) {
        store = ff[i - 1];
        i1 = qMin(j - 1, m1);
        k = i;
        for (ii = 1; ii <= i1; ii++) {
            k++;
            stor1 = ff[k - 1];
            stor2 = aa[i - 1][ii];
            store -= stor1 * stor2;
        }
        stor1 = aa[i - 1][0];
        ff[i - 1] = store / stor1;
        i--;
    }
    /*  solve the system  (aa)' (f2) = f1 */
L230:
    ff[0] /= a[0][0];
    if (*rank == 1)
        goto L260;
    for (j = 2; j <= *rank; j++) {
        store = ff[j - 1];
        i1 = qMin(j - 1, m1);
        k = j;
        for (ii = 1; ii <= i1; ii++) {
            k--;
            stor1 = ff[k - 1];
            stor2 = aa[k - 1][ii];
            store -= stor1 * stor2;
        }
        stor1 = aa[j - 1][0];
        ff[j - 1] = store / stor1;
    }
    /*  premultiply f2 by the transpoze of a. */
L260:
    k = 0;
    for (i = 1; i <= n; i++) {
        store = 0;
        if (a[i - 1][0] > tol)
            k++;
        j1 = qMin(i, m);
        kk = k;
        ij = i + 1;
        for (j = 1; j <= j1; j++) {
            ij--;
            if (a[ij - 1][0] <= tol)
                goto L270;
            stor1 = a[ij - 1][j - 1];
            stor2 = ff[kk - 1];
            store += stor1 * stor2;
            kk--;
        L270:
            ;
        }
        c[i - 1] = store;
    }
    /*  add to the sum of squared residuals the contribution of putting
     to zero the small diagonal elements of matrix (a). */
    stor3 = 0;
    for (i = 1; i <= n; i++) {
        if (a[i - 1][0] > tol)
            goto L310;
        store = f[i - 1];
        i1 = qMin(n - i, m1);
        if (i1 == 0)
            goto L300;
        for (j = 1; j <= i1; j++) {
            ij = i + j;
            stor1 = c[ij - 1];
            stor2 = a[i - 1][j];
            store -= stor1 * stor2;
        }
    L300:
        /*    fac = a[i - 1][0] * c[i - 1]; */
        stor1 = a[i - 1][0];
        stor2 = c[i - 1];
        stor1 *= stor2;
        stor3 += stor1 * (stor1 - store - store);
    L310:
        ;
    }
    fac = stor3;
    *sq += fac;
    delete [] aa[0];
    delete [] h;
    delete [] ff;
    delete [] aa;
}

double fprati(double* p1, double* f1, double p2, double f2, double* p3,
              double* f3)
{
    /*  given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
     gives the value of p such that the rational interpolating function
     of the form r(p) = (u*p+v)/(p+w) equals zero at p. */
    double p, h1, h2, h3;
    if (*p3 > 0) {
        /*  value of p in case p3 ^= infinity. */
        h1 = *f1 * (f2 - *f3);
        h2 = f2 * (*f3 - *f1);
        h3 = *f3 * (*f1 - f2);
        p = -(*p1 * p2 * h3 + p2 * *p3 * h1 + *p3 * *p1 * h2) /
        (*p1 * h1 + p2 * h2 + *p3 * h3);
    } else
    /*  value of p in case p3 = infinity. */
        p = (*p1 * (*f1 - *f3) * f2 - p2 * (f2 - *f3) * *f1) / ((*f1 - f2) * *f3);
    /*  adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0. */
    if (f2 < 0.0) {
        *p3 = p2;
        *f3 = f2;
    } else {
        *p1 = p2;
        *f1 = f2;
    }
    return p;
}

void fprota(double v_cos, double v_sin, double* a, double* b)
{
    /*  subroutine fprota applies a givens rotation to a and b. */
    double stor1, stor2;
    stor1 = *a;
    stor2 = *b;
    *b = v_cos * stor2 + v_sin * stor1;
    *a = v_cos * stor1 - v_sin * stor2;
}

void fpsurf(int iopt, int m, double* x, double* y, const double* z,
            const double* w, double xb, double xe, double yb, double ye,
            int kxx, int kyy, double s, int nxest, int nyest, double eta,
            double tol, int maxit, int nrest, int* nx0, double* tx, int* ny0,
            double* ty, double* c, double* fp, double* fp0, double* fpint,
            double* coord, double* f, double* ff, double** a, double** q,
            double** bx, double** by, double** spx, double** spy, double* h,
            int* ier)
{
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /* part 1: determination of the number of knots and their position.     c */
    /* ****************************************************************     c */
    /* given a set of knots we compute the least-squares spline sinf(x,y),  c */
    /* and the corresponding weighted sum of squared residuals fp=f(p=inf). c */
    /* if iopt=-1  sinf(x,y) is the requested approximation.                c */
    /* if iopt=0 or iopt=1 we check whether we can accept the knots:        c */
    /*   if fp <=s we will continue with the current set of knots.          c */
    /*   if fp > s we will increase the number of knots and compute the     c */
    /*      corresponding least-squares spline until finally  fp<=s.        c */
    /* the initial choice of knots depends on the value of s and iopt.      c */
    /*   if iopt=0 we first compute the least-squares polynomial of degree  c */
    /*     kx in x and ky in y; nx=nminx=2*kx+2 and ny=nminy=2*ky+2.        c */
    /*     fp0=f(0) denotes the corresponding weighted sum of squared       c */
    /*     residuals                                                        c */
    /*   if iopt=1 we start with the knots found at the last call of the    c */
    /*     routine, except for the case that s>=fp0; then we can compute    c */
    /*     the least-squares polynomial directly.                           c */
    /* eventually the independent variables x and y (and the corresponding  c */
    /* parameters) will be switched if this can reduce the bandwidth of the c */
    /* system to be solved.                                                 c */
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    /*  ichang denotes whether(1) or not(-1) the directions have been inter-
     changed. */
    int i, iband = 0, iband1 = 0, iband3, iband4, ichang, ich1, ich3, ii, in,
    irot, iter, i1, i2, j, jrot, jxy, j1, kx, kx1, kx2, ky, ky1, ky2, l, la,
    lf, lh, lx, ly, l1, l2, n, ncof = 0, nk1x = 0, nk1y = 0, nminx, nminy,
    nreg = 0, nrint, num, num1, nx, nxe, nxx, ny, nye, nyy = 0, n1, rank;
    double acc = 0.0, arg, v_cos, dmax, fac1, fac2, fpmax, fpms = 0.0, f1, f2,
    f3, hxi, p, pinv, piv, p1, p2, p3, sigma, v_sin, sq, store, wi, x0,
    x1, y0, y1, zi, eps, d1;
    int* index = new int[nrest];
    int* nummer = new int[m];
    double hx[6], hy[6];
    ichang = -1;
    x0 = xb;
    x1 = xe;
    y0 = yb;
    y1 = ye;
    kx = kxx;
    ky = kyy;
    kx1 = kx + 1;
    ky1 = ky + 1;
    nxe = nxest;
    nye = nyest;
    eps = sqrt(eta);
    if (iopt < 0)
        goto L20;
    /*  calculation of acc, the absolute tolerance for the root of f(p)=s. */
    acc = tol * s;
    if (iopt && (*fp0 > s))
        goto L20;
    /*  initialization for the least-squares polynomial. */
    nminx = kx1 << 1;
    nminy = ky1 << 1;
    nx = nminx;
    ny = nminy;
    *ier = -2;
    goto L30;
L20:
    nx = *nx0;
    ny = *ny0;
    /*  main loop for the different sets of knots. m is a save upper bound
     for the number of trials. */
L30:
    for (iter = 1; iter <= m; iter++) {
        /*  find the position of the additional knots which are needed for the
         b-spline representation of s(x,y). */
        l = nx - 1;
        for (i = 0; i < kx1; i++) {
            tx[i] = x0;
            tx[l--] = x1;
        }
        l = ny - 1;
        for (i = 0; i < ky1; i++) {
            ty[i] = y0;
            ty[l--] = y1;
        }
        /*  find nrint, the total number of knot intervals and nreg, the number
         of panels in which the approximation domain is subdivided by the
         intersection of knots. */
        nxx = nx - (kx1 << 1) + 1;
        nyy = ny - (ky1 << 1) + 1;
        nrint = nxx + nyy;
        nreg = nxx * nyy;
        /*  find the bandwidth of the observation matrix a.
         if necessary, interchange the variables x and y, in order to obtain
         a minimal bandwidth. */
        iband1 = kx * (ny - ky1) + ky;
        l = ky * (nx - kx1) + kx;
        if (iband1 > l) {
            iband1 = l;
            ichang = -ichang;
            for (i = 0; i < m; i++) {
                store = x[i];
                x[i] = y[i];
                y[i] = store;
            }
            store = x0;
            x0 = y0;
            y0 = store;
            store = x1;
            x1 = y1;
            y1 = store;
            n = qMin(nx, ny);
            for (i = 0; i < n; i++) {
                store = tx[i];
                tx[i] = ty[i];
                ty[i] = store;
            }
            n1 = n + 1;
            i = nx - ny;
            if (i == 0)
                goto L120;
            if (i > 0)
                goto L100;
            for (i = n1 - 1; i < ny; i++)
                tx[i] = ty[i];
            goto L120;
        L100:
            for (i = n1 - 1; i < nx; i++)
                ty[i] = tx[i];
        L120:
            l = nx;
            nx = ny;
            ny = l;
            l = nxe;
            nxe = nye;
            nye = l;
            l = nxx;
            nxx = nyy;
            nyy = l;
            l = kx;
            kx = ky;
            ky = l;
            kx1 = kx + 1;
            ky1 = ky + 1;
        }
        iband = iband1 + 1;
        /*  arrange the data points according to the panel they belong to. */
        fporde(x, y, m, kx, ky, tx, nx, ty, ny, nummer, index, nreg);
        /*  find ncof, the number of b-spline coefficients. */
        nk1x = nx - kx1;
        nk1y = ny - ky1;
        ncof = nk1x * nk1y;
        /*  initialize the observation matrix a. */
        for (i = 0; i < ncof; i++) {
            f[i] = 0;
            for (j = 0; j < iband; j++)
                a[i][j] = 0;
        }
        /*  initialize the sum of squared residuals. */
        *fp = 0;
        /*  fetch the data points in the new order. main loop for the
         different panels. */
        for (num = 1; num <= nreg; num++) {
            /*  fix certain constants for the current panel; jrot records the column
             number of the first non-zero element in a row of the observation
             matrix according to a data point of the panel. */
            num1 = num - 1;
            lx = num1 / nyy;
            l1 = lx + kx1;
            ly = num1 - lx * nyy;
            l2 = ly + ky1;
            jrot = lx * nk1y + ly;
            /*  test whether there are still data points in the panel. */
            in = index[num1];
        L150:
            if (in) {
                /*  fetch a new data point. */
                wi = w[in - 1];
                zi = z[in - 1] * wi;
                /*  evaluate for the x-direction, the (kx+1) non-zero b-splines at x(in). */
                fpbspl(tx, kx, x[in - 1], l1, hx);
                /*  evaluate for the y-direction, the (ky+1) non-zero b-splines at y(in). */
                fpbspl(ty, ky, y[in - 1], l2, hy);
                /*  store the value of these b-splines in spx and spy respectively. */
                for (i = 0; i < kx1; i++)
                    spx[in - 1][i] = hx[i];
                for (i = 0; i < ky1; i++)
                    spy[in - 1][i] = hy[i];
                /*  initialize the new row of observation matrix. */
                for (i = 0; i < iband; i++)
                    h[i] = 0;
                /*  calculate the non-zero elements of the new row by making the cross
                 products of the non-zero b-splines in x- and y-direction. */
                i1 = 0;
                for (i = 0; i < kx1; i++) {
                    hxi = hx[i];
                    j1 = i1;
                    for (j = 0; j < ky1; j++)
                        h[j1++] = hxi * hy[j] * wi;
                    i1 += nk1y;
                }
                /*  rotate the row into triangle by givens transformations . */
                irot = jrot - 1;
                for (i = 1; i <= iband; i++) {
                    irot++;
                    piv = h[i - 1];
                    if (piv) {
                        /*  calculate the parameters of the givens transformation. */
                        fpgivs(piv, &a[irot][0], &v_cos, &v_sin);
                        /*  apply that transformation to the right hand side. */
                        fprota(v_cos, v_sin, &zi, &f[irot]);
                        if (i == iband)
                            break;
                        /*  apply that transformation to the left hand side. */
                        i2 = 1;
                        for (j = i; j < iband; j++)
                            fprota(v_cos, v_sin, &h[j], &a[irot][i2++]);
                    }
                }
                /*  add the contribution of the row to the sum of squares of residual
                 right hand sides. */
                *fp += zi * zi;
                /*  find the number of the next data point in the panel. */
                in = nummer[in - 1];
                goto L150;
            }
        }
        /*  find dmax, the maximum value for the diagonal elements in the reduced
         triangle. */
        dmax = 0;
        for (i = 0; i < ncof; i++)
            dmax = qMax(a[i][0], dmax);
        /*  check whether the observation matrix is rank deficient. */
        sigma = eps * dmax;
        for (i = 0; i < ncof; i++)
            if (a[i][0] <= sigma)
                goto L280;
        /*  backward substitution in case of full rank. */
        fpback(a, f, ncof, iband, c);
        rank = ncof;
        for (i = 0; i < ncof; i++)
            q[i][0] = a[i][0] / dmax;
        goto L300;
        /*  in case of rank deficiency, find the minimum norm solution.
         check whether there is sufficient working space */
    L280:
        for (i = 0; i < ncof; i++) {
            ff[i] = f[i];
            for (j = 0; j < iband; j++)
                q[i][j] = a[i][j];
        }
        fprank(q, ff, ncof, iband, sigma, c, &sq, &rank);
        for (i = 0; i <= ncof; i++)
            q[i][0] = q[i][0] / dmax;
        /*  add to the sum of squared residuals, the contribution of reducing
         the rank. */
        *fp += sq;
    L300:
        if (*ier == -2)
            *fp0 = *fp;
        /*  test whether the least-squares spline is an acceptable solution. */
        if (iopt < 0)
            goto L820;
        fpms = *fp - s;
        if (fabs(fpms) <= acc) {
            if (*fp <= 0)
                goto L815;
            else
                goto L820;
        }
        /*  test whether we can accept the choice of knots. */
        if (fpms < 0)
            goto L430;
        /*  test whether we cannot further increase the number of knots. */
        if (ncof > m)
            goto L790;
        *ier = 0;
        /*  search where to add a new knot.
         find for each interval the sum of squared residuals fpint for the
         data points having the coordinate belonging to that knot interval.
         calculate also coord which is the same sum, weighted by the position
         of the data points considered. */
        for (i = 0; i < nrint; i++) {
            fpint[i] = 0;
            coord[i] = 0;
        }
        for (num = 1; num <= nreg; num++) {
            num1 = num - 1;
            lx = num1 / nyy;
            ly = num1 - lx * nyy;
            l2 = ly + nxx;
            jrot = lx * nk1y + ly;
            in = index[num1];
            while (in) {
                store = 0;
                i1 = jrot;
                for (i = 0; i < kx1; i++) {
                    hxi = spx[in - 1][i];
                    j1 = i1;
                    for (j = 0; j < ky1; j++)
                        store += hxi * spy[in - 1][j] * c[j1++];
                    i1 += nk1y;
                }
                d1 = w[in - 1] * (z[in - 1] - store);
                store = d1 * d1;
                fpint[lx] += store;
                coord[lx] += store * x[in - 1];
                fpint[l2] += store;
                coord[l2] += store * y[in - 1];
                in = nummer[in - 1];
            }
        }
        /*  find the interval for which fpint is maximal on the condition that
         there still can be added a knot. */
    L370:
        l = 0;
        fpmax = 0;
        l1 = 1;
        l2 = nrint;
        if (nx == nxe)
            l1 = nxx + 1;
        if (ny == nye)
            l2 = nxx;
        if (l1 > l2)
            goto L810;
        for (i = l1; i <= l2; i++)
            if (fpmax < fpint[i - 1]) {
                l = i;
                fpmax = fpint[i - 1];
            }
        /*  test whether we cannot further increase the number of knots. */
        if (l == 0)
            goto L785;
        /*  calculate the position of the new knot. */
        arg = coord[l - 1] / fpint[l - 1];
        /*  test in what direction the new knot is going to be added. */
        if (l <= nxx) {
            /*  addition in the x-direction. */
            jxy = l + kx1;
            fpint[l - 1] = 0;
            fac1 = tx[jxy - 1] - arg;
            fac2 = arg - tx[jxy - 2];
            if (fac1 > 10 * fac2 || fac2 > 10 * fac1)
                goto L370;
            j = nx;
            for (i = jxy; i <= nx; i++) {
                tx[j] = tx[j - 1];
                --j;
            }
            tx[jxy - 1] = arg;
            nx++;
        } else {
            /*  addition in the y-direction. */
            jxy = l + ky1 - nxx;
            fpint[l - 1] = 0;
            fac1 = ty[jxy - 1] - arg;
            fac2 = arg - ty[jxy - 2];
            if (fac1 > 10 * fac2 || fac2 > 10 * fac1)
                goto L370;
            j = ny;
            for (i = jxy; i <= ny; i++) {
                ty[j] = ty[j - 1];
                --j;
            }
            ty[jxy-1] = arg;
            ny++;
        }
        /*  restart the computations with the new set of knots. */
    }
    /*  test whether the least-squares polynomial is a solution of our
     approximation problem. */
L430:
    if (*ier == -2)
        goto L830;
    /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     part 2: determination of the smoothing spline sp(x,y)                c
     *****************************************************                c
     we have determined the number of knots and their position. we now    c
     compute the b-spline coefficients of the smoothing spline sp(x,y).   c
     the observation matrix a is extended by the rows of a matrix,        c
     expressing that sp(x,y) must be a polynomial of degree kx in x and   c
     ky in y. the corresponding weights of these additional rows are set  c
     to 1./p.  iteratively we than have to determine the value of p       c
     such that f(p)=sum((w(i)*(z(i)-sp(x(i),y(i))))**2) be = s.           c
     we already know that the least-squares polynomial corresponds to     c
     p=0  and that the least-squares spline corresponds to p=infinity.    c
     the iteration process which is proposed here makes use of rational   c
     interpolation. since f(p) is a convex and strictly decreasing        c
     function of p, it can be approximated by a rational function r(p)=   c
     (u*p+v)/(p+w). three values of p(p1,p2,p3) with corresponding values c
     of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the c
     new value of p such that r(p)=s. convergence is guaranteed by taking c
     f1 > 0 and f3 < 0.                                                   c
     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    kx2 = kx1 + 1;
    /*  test whether there are interior knots in the x-direction. */
    if (nk1x == kx1)
        goto L440;
    /*  evaluate the discotinuity jumps of the kx-th order derivative of
     the b-splines at the knots tx(l),l=kx+2,...,nx-kx-1. */
    fpdisc(tx, nx, kx2, bx);
L440:
    ky2 = ky1 + 1;
    /*  test whether there are interior knots in the y-direction. */
    if (nk1y == ky1)
        goto L450;
    /*  evaluate the discontinuity jumps of the ky-th order derivative of
     the b-splines at the knots ty(l),l=ky+2,...,ny-ky-1. */
    fpdisc(ty, ny, ky2, by);
    /*  initial value for p. */
L450:
    p1 = 0;
    f1 = *fp0 - s;
    p3 = -1;
    f3 = fpms;
    p = 0;
    for (i = 0; i < ncof; i++)
        p += a[i][0];
    p = ncof / p;
    /*  find the bandwidth of the extended observation matrix. */
    iband3 = kx1 * nk1y;
    iband4 = iband3 + 1;
    ich1 = 0;
    ich3 = 0;
    /*  iteration process to find the root of f(p)=s. */
    for (iter = 1; iter <= maxit; iter++) {
        pinv = 1 / p;
        /*  store the triangularized observation matrix into q. */
        for (i = 0; i < ncof; i++) {
            ff[i] = f[i];
            for (j = 0; j < iband; j++)
                q[i][j] = a[i][j];
            for (j = iband; j < iband4; j++)
                q[i][j] = 0;
        }
        if (nk1y == ky1)
            goto L560;
        /*  extend the observation matrix with the rows of a matrix, expressing
         that for x=cst. sp(x,y) must be a polynomial in y of degree ky. */
        for (i = ky2; i <= nk1y; i++) {
            ii = i - ky1;
            for (j = 1; j <= nk1x; j++) {
                /*  initialize the new row. */
                for (l = 0; l < iband; l++)
                    h[l] = 0;
                /*  fill in the non-zero elements of the row. jrot records the column
                 number of the first non-zero element in the row. */
                for (l = 0; l < ky2; l++)
                    h[l] = by[ii - 1][l] * pinv;
                zi = 0;
                jrot = (j - 1) * nk1y + ii;
                /*  rotate the new row into triangle by givens transformations without
                 square roots. */
                for (irot = jrot; irot <= ncof; irot++) {
                    piv = h[0];
                    i2 = qMin(iband1, ncof - irot);
                    if (piv) {
                        /*  calculate the parameters of the givens transformation. */
                        fpgivs(piv, q[irot - 1], &v_cos, &v_sin);
                        /*  apply that givens transformation to the right hand side. */
                        fprota(v_cos, v_sin, &zi, &ff[irot - 1]);
                        if (i2)
                        /*  apply that givens transformation to the left hand side. */
                            for (l = 1; l <= i2; l++)
                                fprota(v_cos, v_sin, &h[l], &q[irot - 1][l]);
                    }
                    if (i2 > 0) {
                        for (l = 0; l < i2; l++)
                            h[l] = h[l + 1];
                        h[i2] = 0;
                    }
                }
            }
        }
    L560:
        if (nk1x != kx1) {
            /*  extend the observation matrix with the rows of a matrix expressing
             that for y=cst. sp(x,y) must be a polynomial in x of degree kx. */
            for (i = kx2; i <= nk1x; i++) {
                ii = i - kx1;
                for (j = 1; j <= nk1y; j++) {
                    /*  initialize the new row */
                    for (l = 0; l < iband4; l++)
                        h[l] = 0;
                    /*  fill in the non-zero elements of the row. jrot records the column
                     number of the first non-zero element in the row. */
                    j1 = 0;
                    for (l = 0; l < kx2; l++) {
                        h[j1] = bx[ii - 1][l] * pinv;
                        j1 += nk1y;
                    }
                    zi = 0;
                    jrot = (i - kx2) * nk1y + j;
                    /*  rotate the new row into triangle by givens transformations . */
                    for (irot = jrot; irot <= ncof; irot++) {
                        piv = h[0];
                        i2 = qMin(iband3, ncof - irot);
                        if (piv) {
                            /*  calculate the parameters of the givens transformation. */
                            fpgivs(piv, q[irot - 1], &v_cos, &v_sin);
                            /*  apply that givens transformation to the right hand side. */
                            fprota(v_cos, v_sin, &zi, &ff[irot - 1]);
                            if (i2)
                            /*  apply that givens transformation to the left hand side. */
                                for (l = 1; l <= i2; l++)
                                    fprota(v_cos, v_sin, &h[l], &q[irot - 1][l]);
                        }
                        if (i2 > 0) {
                            for (l = 0; l < i2; l++)
                                h[l] = h[l + 1];
                            h[i2] = 0;
                        }
                    }
                }
            }
        }
        /*  find dmax, the maximum value for the diagonal elements in the
         reduced triangle. */
        dmax = 0;
        for (i = 0; i < ncof; i++)
            dmax = qMax(q[i][0], dmax);
        /*  check whether the matrix is rank deficient. */
        sigma = eps * dmax;
        for (i = 0; i < ncof; i++)
            if (q[i][0] <= sigma)
                goto L670;
        /*  backward substitution in case of full rank. */
        fpback(q, ff, ncof, iband4, c);
        rank = ncof;
        goto L675;
        /*  in case of rank deficiency, find the minimum norm solution. */
    L670:
        lf = 1;
        lh = lf + ncof;
        la = lh + iband4;
        fprank(q, ff, ncof, iband4, sigma, c, &sq, &rank);
    L675:
        for (i = 0; i < ncof; i++)
            q[i][0] = q[i][0] / dmax;
        /*  compute f(p). */
        *fp = 0;
        for (num = 1; num <= nreg; num++) {
            num1 = num - 1;
            lx = num1 / nyy;
            ly = num1 - lx * nyy;
            jrot = lx * nk1y + ly;
            in = index[num - 1];
            while (in) {
                store = 0;
                i1 = jrot;
                for (i = 1; i <= kx1; i++) {
                    hxi = spx[in - 1][i - 1];
                    j1 = i1;
                    for (j = 0; j < ky1; j++)
                        store += hxi * spy[in - 1][j] * c[j1++];
                    i1 += nk1y;
                }
                d1 = w[in - 1] * (z[in - 1] - store);
                *fp += d1 * d1;
                in = nummer[in - 1];
            }
        }
        /*  test whether the approximation sp(x,y) is an acceptable solution. */
        fpms = *fp - s;
        if (fabs(fpms) <= acc)
            goto L820;
        /*  test whether the maximum allowable number of iterations has been
         reached. */
        if (iter == maxit)
            goto L795;
        /*  carry out one more step of the iteration process. */
        p2 = p;
        f2 = fpms;
        if (ich3 == 0) {
            if ((f2 - f3) <= acc) {
                /*  our initial choice of p is too large. */
                p3 = p2;
                f3 = f2;
                p *= 0.04;
                if (p <= p1)
                    p = p1 * 0.9 + p2 * 0.1;
                continue;
            }
            if (f2 < 0)
                ich3 = 1;
        }
        if (ich1 != 0)
            goto L760;
        if (f1 - f2 > acc)
            goto L750;
        /*  our initial choice of p is too small */
        p1 = p2;
        f1 = f2;
        p /= 0.04;
        if (p3 < 0)
            continue;
        if (p >= p3)
            p = p2 * 0.1 + p3 * 0.9;
        continue;
    L750:
        if (f2 > 0)
            ich1 = 1;
        /*  test whether the iteration process proceeds as theoretically */
        /*  expected. */
    L760:
        if (f2 >= f1 || f2 <= f3)
            goto L800;
        /*  find the new value of p. */
        p = fprati(&p1, &f1, p2, f2, &p3, &f3);
    }
    /*  error codes and messages. */
L785:
    *ier = 5;
    goto L830;
L790:
    *ier = 4;
    goto L830;
L795:
    *ier = 3;
    goto L830;
L800:
    *ier = 2;
    goto L830;
L810:
    *ier = 1;
    goto L830;
L815:
    *ier = -1;
    *fp = 0;
L820:
    if (ncof != rank)
        *ier = -rank;
    /*  test whether x and y are in the original order. */
L830:
    if (ichang < 0)
        goto L930;
    /*  if not, interchange x and y once more. */
    l1 = 0;
    for (i = 0; i < nk1x; i++) {
        l2 = i;
        for (j = 0; j < nk1y; j++) {
            f[l2] = c[l1++];
            l2 += nk1x;
        }
    }
    for (i = 0; i < ncof; i++)
        c[i] = f[i];
    for (i = 0; i < m; i++) {
        store = x[i];
        x[i] = y[i];
        y[i] = store;
    }
    n = qMin(nx, ny);
    for (i = 0; i < n; i++) {
        store = tx[i];
        tx[i] = ty[i];
        ty[i] = store;
    }
    n1 = n + 1;
    i = nx - ny;
    if (i < 0)
        goto L880;
    else
        if (i == 0)
            goto L920;
        else
            goto L900;
L880:
    for (i = n1 - 1; i < ny; i++)
        tx[i] = ty[i];
    goto L920;
L900:
    for (i = n1 - 1; i < nx; i++)
        ty[i] = tx[i];
L920:
    l = nx;
    nx = ny;
    ny = l;
L930:
    if (iopt >= 0) {
        *nx0 = nx;
        *ny0 = ny;
    }
    delete [] nummer;
    delete [] index;
    return;
}

FitPack::FitPack()
{
}

FitPack::~FitPack()
{
}

void FitPack::bispev(const double* tx, int nx, const double* ty, int ny,
                     const double* c, int kx, int ky, const double* x, int mx,
                     const double* y, int my, double *z, int* ier)
{
    /*  subroutine bispev evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,...
     ,my a bivariate spline s(x,y) of degrees kx and ky, given in the
     b-spline representation.

     calling sequence:
     call bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,
     * iwrk,kwrk,ier)

     input parameters:
     tx    : real array, length nx, which contains the position of the
     knots in the x-direction.
     nx    : integer, giving the total number of knots in the x-direction
     ty    : real array, length ny, which contains the position of the
     knots in the y-direction.
     ny    : integer, giving the total number of knots in the y-direction
     c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the
     b-spline coefficients.
     kx,ky : integer values, giving the degrees of the spline.
     x     : real array of dimension (mx).
     before entry x(i) must be set to the x co-ordinate of the
     i-th grid point along the x-axis.
     tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx.
     mx    : on entry mx must specify the number of grid points along
     the x-axis. mx >=1.
     y     : real array of dimension (my).
     before entry y(j) must be set to the y co-ordinate of the
     j-th grid point along the y-axis.
     ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my.
     my    : on entry my must specify the number of grid points along
     the y-axis. my >=1.
     wrk   : real array of dimension lwrk. used as workspace.
     lwrk  : integer, specifying the dimension of wrk.
     lwrk >= mx*(kx+1)+my*(ky+1)
     iwrk  : integer array of dimension kwrk. used as workspace.
     kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my.

     output parameters:
     z     : real array of dimension (mx*my).
     on succesful exit z(my*(i-1)+j) contains the value of s(x,y)
     at the point (x(i),y(j)),i=1,...,mx;j=1,...,my.
     ier   : integer error flag
     ier=0 : normal return
     ier=10: invalid input data (see restrictions)

     restrictions:
     mx >=1, my >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my
     tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx
     ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my

     other subroutines required:
     fpbisp,fpbspl

     references :
     de boor c : on calculating with b-splines, j. approximation theory
     6 (1972) 50-62.
     cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
     applics 10 (1972) 134-149.
     dierckx p. : curve and surface fitting with splines, monographs on
     numerical analysis, oxford university press, 1993.

     author :
     p.dierckx
     dept. computer science, k.u.leuven
     celestijnenlaan 200a, b-3001 heverlee, belgium.
     e-mail : Paul.Dierckx@cs.kuleuven.ac.be

     latest update : march 1987 */
    /*  C++ translation by Werner Stille. */

    /*  before starting computations a data check is made. if the input data
     are invalid control is immediately repassed to the calling program. */
    int i, i1;
    *ier = 10;
    i1 = mx - 1;
    if (i1 < 0)
        return;
    if (i1)
        for (i = 1; i < mx; ++i)
            if (x[i] < x[i - 1])
                return;
    i1 = my - 1;
    if (i1 < 0)
        return;
    if (i1)
        for (i = 1; i < my; ++i)
            if (y[i] < y[i - 1])
                return;
    *ier = 0;
    fpbisp(tx, nx, ty, ny, c, kx, ky, x, mx, y, my, z);
    return;
}

void FitPack::curfit(int iopt, int m, const double* x, const double* y,
                     const double* w, double xb, double xe, int k, double s,
                     int nest, int* n, double* t, double* c, double* fp,
                     double* wrk, int lwrk, int* iwrk, int* ier)
{
    /*  given the set of data points (x(i),y(i)) and the set of positive
     numbers w(i),i=1,2,...,m,subroutine curfit determines a smooth spline
     approximation of degree k on the interval xb <= x <= xe.
     if iopt=-1 curfit calculates the weighted least-squares spline
     according to a given set of knots.
     if iopt>=0 the number of knots of the spline s(x) and the position
     t(j),j=1,2,...,n is chosen automatically by the routine. the smooth-
     ness of s(x) is then achieved by minimalizing the discontinuity
     jumps of the k-th derivative of s(x) at the knots t(j),j=k+2,k+3,...,
     n-k-1. the amount of smoothness is determined by the condition that
     f(p)=sum((w(i)*(y(i)-s(x(i))))**2) be <= s, with s a given non-
     negative constant, called the smoothing factor.
     the fit s(x) is given in the b-spline representation (b-spline coef-
     ficients c(j),j=1,2,...,n-k-1) and can be evaluated by means of
     subroutine splev.

     calling sequence:
     call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,
     * lwrk,iwrk,ier)

     parameters:
     iopt  : integer flag. on entry iopt must specify whether a weighted
     least-squares spline (iopt=-1) or a smoothing spline (iopt=
     0 or 1) must be determined. if iopt=0 the routine will start
     with an initial set of knots t(i)=xb, t(i+k+1)=xe, i=1,2,...
     k+1. if iopt=1 the routine will continue with the knots
     found at the last call of the routine.
     attention: a call with iopt=1 must always be immediately
     preceded by another call with iopt=1 or iopt=0.
     unchanged on exit.
     m     : integer. on entry m must specify the number of data points.
     m > k. unchanged on exit.
     x     : real array of dimension at least (m). before entry, x(i)
     must be set to the i-th value of the independent variable x,
     for i=1,2,...,m. these values must be supplied in strictly
     ascending order. unchanged on exit.
     y     : real array of dimension at least (m). before entry, y(i)
     must be set to the i-th value of the dependent variable y,
     for i=1,2,...,m. unchanged on exit.
     w     : real array of dimension at least (m). before entry, w(i)
     must be set to the i-th value in the set of weights. the
     w(i) must be strictly positive. unchanged on exit.
     see also further comments.
     xb,xe : real values. on entry xb and xe must specify the boundaries
     of the approximation interval. xb<=x(1), xe>=x(m).
     unchanged on exit.
     k     : integer. on entry k must specify the degree of the spline.
     1<=k<=5. it is recommended to use cubic splines (k=3).
     the user is strongly dissuaded from choosing k even,together
     with a small s-value. unchanged on exit.
     s     : real.on entry (in case iopt>=0) s must specify the smoothing
     factor. s >=0. unchanged on exit.
     for advice on the choice of s see further comments.
     nest  : integer. on entry nest must contain an over-estimate of the
     total number of knots of the spline returned, to indicate
     the storage space available to the routine. nest >=2*k+2.
     in most practical situation nest=m/2 will be sufficient.
     always large enough is  nest=m+k+1, the number of knots
     needed for interpolation (s=0). unchanged on exit.
     n     : integer.
     unless ier =10 (in case iopt >=0), n will contain the
     total number of knots of the spline approximation returned.
     if the computation mode iopt=1 is used this value of n
     should be left unchanged between subsequent calls.
     in case iopt=-1, the value of n must be specified on entry.
     t     : real array of dimension at least (nest).
     on succesful exit, this array will contain the knots of the
     spline,i.e. the position of the interior knots t(k+2),t(k+3)
     ...,t(n-k-1) as well as the position of the additional knots
     t(1)=t(2)=...=t(k+1)=xb and t(n-k)=...=t(n)=xe needed for
     the b-spline representation.
     if the computation mode iopt=1 is used, the values of t(1),
     t(2),...,t(n) should be left unchanged between subsequent
     calls. if the computation mode iopt=-1 is used, the values
     t(k+2),...,t(n-k-1) must be supplied by the user, before
     entry. see also the restrictions (ier=10).
     c     : real array of dimension at least (nest).
     on succesful exit, this array will contain the coefficients
     c(1),c(2),..,c(n-k-1) in the b-spline representation of s(x)
     fp    : real. unless ier=10, fp contains the weighted sum of
     squared residuals of the spline approximation returned.
     wrk   : real array of dimension at least (m*(k+1)+nest*(7+3*k)).
     used as working space. if the computation mode iopt=1 is
     used, the values wrk(1),...,wrk(n) should be left unchanged
     between subsequent calls.
     lwrk  : integer. on entry,lwrk must specify the actual dimension of
     the array wrk as declared in the calling (sub)program.lwrk
     must not be too small (see wrk). unchanged on exit.
     iwrk  : integer array of dimension at least (nest).
     used as working space. if the computation mode iopt=1 is
     used,the values iwrk(1),...,iwrk(n) should be left unchanged
     between subsequent calls.
     ier   : integer. unless the routine detects an error, ier contains a
     non-positive value on exit, i.e.
     ier=0  : normal return. the spline returned has a residual sum of
     squares fp such that abs(fp-s)/s <= tol with tol a relat-
     ive tolerance set to 0.001 by the program.
     ier=-1 : normal return. the spline returned is an interpolating
     spline (fp=0).
     ier=-2 : normal return. the spline returned is the weighted least-
     squares polynomial of degree k. in this extreme case fp
     gives the upper bound fp0 for the smoothing factor s.
     ier=1  : error. the required storage space exceeds the available
     storage space, as specified by the parameter nest.
     probably causes : nest too small. if nest is already
     large (say nest > m/2), it may also indicate that s is
     too small
     the approximation returned is the weighted least-squares
     spline according to the knots t(1),t(2),...,t(n). (n=nest)
     the parameter fp gives the corresponding weighted sum of
     squared residuals (fp>s).
     ier=2  : error. a theoretically impossible result was found during
     the iteration proces for finding a smoothing spline with
     fp = s. probably causes : s too small.
     there is an approximation returned but the corresponding
     weighted sum of squared residuals does not satisfy the
     condition abs(fp-s)/s < tol.
     ier=3  : error. the maximal number of iterations maxit (set to 20
     by the program) allowed for finding a smoothing spline
     with fp=s has been reached. probably causes : s too small
     there is an approximation returned but the corresponding
     weighted sum of squared residuals does not satisfy the
     condition abs(fp-s)/s < tol.
     ier=10 : error. on entry, the input data are controlled on validity
     the following restrictions must be satisfied.
     -1<=iopt<=1, 1<=k<=5, m>k, nest>2*k+2, w(i)>0,i=1,2,...,m
     xb<=x(1)<x(2)<...<x(m)<=xe, lwrk>=(k+1)*m+nest*(7+3*k)
     if iopt=-1: 2*k+2<=n<=min(nest,m+k+1)
     xb<t(k+2)<t(k+3)<...<t(n-k-1)<xe
     the schoenberg-whitney conditions, i.e. there
     must be a subset of data points xx(j) such that
     t(j) < xx(j) < t(j+k+1), j=1,2,...,n-k-1
     if iopt>=0: s>=0
     if s=0 : nest >= m+k+1
     if one of these conditions is found to be violated,control
     is immediately repassed to the calling program. in that
     case there is no approximation returned.

     further comments:
     by means of the parameter s, the user can control the tradeoff
     between closeness of fit and smoothness of fit of the approximation.
     if s is too large, the spline will be too smooth and signal will be
     lost ; if s is too small the spline will pick up too much noise. in
     the extreme cases the program will return an interpolating spline if
     s=0 and the weighted least-squares polynomial of degree k if s is
     very large. between these extremes, a properly chosen s will result
     in a good compromise between closeness of fit and smoothness of fit.
     to decide whether an approximation, corresponding to a certain s is
     satisfactory the user is highly recommended to inspect the fits
     graphically.
     recommended values for s depend on the weights w(i). if these are
     taken as 1/d(i) with d(i) an estimate of the standard deviation of
     y(i), a good s-value should be found in the range (m-sqrt(2*m),m+
     sqrt(2*m)). if nothing is known about the statistical error in y(i)
     each w(i) can be set equal to one and s determined by trial and
     error, taking account of the comments above. the best is then to
     start with a very large value of s ( to determine the least-squares
     polynomial and the corresponding upper bound fp0 for s) and then to
     progressively decrease the value of s ( say by a factor 10 in the
     beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
     approximation shows more detail) to obtain closer fits.
     to economize the search for a good s-value the program provides with
     different modes of computation. at the first call of the routine, or
     whenever he wants to restart with the initial set of knots the user
     must set iopt=0.
     if iopt=1 the program will continue with the set of knots found at
     the last call of the routine. this will save a lot of computation
     time if curfit is called repeatedly for different values of s.
     the number of knots of the spline returned and their location will
     depend on the value of s and on the complexity of the shape of the
     function underlying the data. but, if the computation mode iopt=1
     is used, the knots returned may also depend on the s-values at
     previous calls (if these were smaller). therefore, if after a number
     of trials with different s-values and iopt=1, the user can finally
     accept a fit as satisfactory, it may be worthwhile for him to call
     curfit once more with the selected value for s but now with iopt=0.
     indeed, curfit may then return an approximation of the same quality
     of fit but with fewer knots and therefore better if data reduction
     is also an important objective for the user.

     other subroutines required:
     fpback,fpbspl,fpchec,fpcurf,fpdisc,fpgivs,fpknot,fprati,fprota

     references:
     dierckx p. : an algorithm for smoothing, differentiation and integ-
     ration of experimental data using spline functions,
     j.comp.appl.maths 1 (1975) 165-184.
     dierckx p. : a fast algorithm for smoothing data on a rectangular
     grid while using spline functions, siam j.numer.anal.
     19 (1982) 1286-1304.
     dierckx p. : an improved algorithm for curve fitting with spline
     functions, report tw54, dept. computer science,k.u.
     leuven, 1981.
     dierckx p. : curve and surface fitting with splines, monographs on
     numerical analysis, oxford university press, 1993.

     author:
     p.dierckx
     dept. computer science, k.u. leuven
     celestijnenlaan 200a, b-3001 heverlee, belgium.
     e-mail : Paul.Dierckx@cs.kuleuven.ac.be

     creation date : may 1979
     latest update : march 1987 */
    /*  C++ translation by Werner Stille. */
    int i, ia, ib, ifp, ig, iq, iz, j, k1, k2, lwest, maxit, nmin;
    double tol;
    double **a, **b, **g, **q;
    /*  we set up the parameters tol and maxit */
    maxit = 20;
    tol = 0.001;
    /*  before starting computations a data check is made. if the input data
     are invalid, control is immediately repassed to the calling program. */
    *ier = 10;
    if (k <= 0 || k > 5)
        return;
    k1 = k + 1;
    k2 = k1 + 1;
    if (iopt < -1 || iopt > 1)
        return;
    nmin = k1 << 1;
    if (m < k1 || nest < nmin)
        return;
    lwest = m * k1 + nest * (k * 3 + 7);
    if (lwrk < lwest)
        return;
    if (xb > x[0] || xe < x[m - 1] || w[0] <= 0)
        return;
    for (i = 1; i < m; i++)
        if (x[i - 1] >= x[i] || w[i] <= 0)
            return;
    if (iopt >= 0) {
        if (s < 0)
            return;
        if (s == 0 && nest < m + k1)
            return;
        *ier = 0;
    } else {
        if (*n < nmin || *n > nest)
            return;
        j = *n;
        for (i = 0; i < k1; i++) {
            t[i] = xb;
            t[--j] = xe;
        }
        fpchec(x, m, t, *n, k, ier);
        if (*ier)
            return;
    }
    /* we partition the working space and determine the spline approximation. */
    ifp = 0;
    iz = ifp + nest;
    ia = iz + nest;
    ib = ia + nest * k1;
    ig = ib + nest * k2;
    iq = ig + nest * k2;
    a = new double*[nest];
    b = new double*[nest];
    g = new double*[nest];
    q = new double*[nest];
    for (i = 0; i < nest; i++) {
        a[i] = &wrk[ia + i * k1];
        b[i] = &wrk[ib + i * k2];
        g[i] = &wrk[ig + i * k2];
    }
    for (i = 0; i < m; i++)
        q[i] = &wrk[iq + i * k1];
    fpcurf(iopt, x, y, w, m, xb, xe, k, s, nest, tol, maxit, k1, k2, n, t, c,
           fp, &wrk[ifp], &wrk[iz], a, b, g, q, iwrk, ier);
    delete [] q;
    delete [] g;
    delete [] b;
    delete [] a;
    return;
}

void FitPack::splder(const double* t, int n, const double* c, int k, int nu,
                     const double* x, double* y, int m, int* ier)
{
    /*  subroutine splder evaluates in a number of points x(i),i=1,2,...,m
     the derivative of order nu of a spline s(x) of degree k,given in
     its b-spline representation.

     calling sequence:
     call splder(t,n,c,k,nu,x,y,m,wrk,ier)

     input parameters:
     t    : array,length n, which contains the position of the knots.
     n    : integer, giving the total number of knots of s(x).
     c    : array,length n, which contains the b-spline coefficients.
     k    : integer, giving the degree of s(x).
     nu   : integer, specifying the order of the derivative. 0<=nu<=k
     x    : array,length m, which contains the points where the deriv-
     ative of s(x) must be evaluated.
     m    : integer, giving the number of points where the derivative
     of s(x) must be evaluated
     wrk  : real array of dimension n. used as working space.

     output parameters:
     y    : array,length m, giving the value of the derivative of s(x)
     at the different points.
     ier  : error flag
     ier = 0 : normal return
     ier =10 : invalid input data (see restrictions)

     restrictions:
     0 <= nu <= k
     m >= 1
     t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.

     other subroutines required: fpbspl

     references :
     de boor c : on calculating with b-splines, j. approximation theory
     6 (1972) 50-62.
     cox m.g.  : the numerical evaluation of b-splines, j. inst. maths
     applics 10 (1972) 134-149.
     dierckx p. : curve and surface fitting with splines, monographs on
     numerical analysis, oxford university press, 1993.

     author :
     p.dierckx
     dept. computer science, k.u.leuven
     celestijnenlaan 200a, b-3001 heverlee, belgium.
     e-mail : Paul.Dierckx@cs.kuleuven.ac.be

     latest update : march 1987 */
    /*  C++ translation by Werner Stille. */
    int i, j, kk, k1, k2, l, ll, l1, l2, nk1, nk2, nn;
    double ak, arg, fac, sp, tb, te;
    double *wrk;
    double h[6];
    *ier = 10;
    if (nu < 0 || nu > k)
        return;
    if (m < 1)
        return;
    for (i = 1; i < m; ++i)
        if (x[i] < x[i - 1])
            return;
    *ier = 0;
    /*  fetch tb and te, the boundaries of the approximation interval. */
    k1 = k + 1;
    nk1 = n - k1;
    tb = t[k];
    te = t[nk1];
    /*  the derivative of order nu of a spline of degree k is a spline of
     degree k-nu,the b-spline coefficients wrk(i) of which can be found
     using the recurrence scheme of de boor. */
    l = 1;
    kk = k;
    nn = n;
    wrk = new double[n];
    for (i = 0; i < nk1; ++i)
        wrk[i] = c[i];
    if (nu) {
        nk2 = nk1;
        for (j = 0; j < nu; ++j) {
            ak = (double) kk;
            --nk2;
            l1 = l;
            for (i = 1; i <= nk2; ++i) {
                ++l1;
                l2 = l1 + kk;
                fac = t[l2 - 1] - t[l1 - 1];
                if (fac > 0.0)
                    wrk[i - 1] = ak * (wrk[i] - wrk[i - 1]) / fac;
            }
            ++l;
            --kk;
        }
        if (!kk) {
            /*  if nu=k the derivative is a piecewise constant function */
            j = 0;
            for (i = 0; i < m; ++i) {
                arg = x[i];
                while (!(arg < t[l] || l == nk1)) {
                    ++l;
                    ++j;
                }
                y[i] = wrk[j];
            }
            delete [] wrk;
            return;
        }
    }
    l = k1;
    l1 = l + 1;
    k2 = k1 - nu;
    /*  main loop for the different points. */
    for (i = 0; i < m; ++i) {
        arg = x[i];
        if (arg < tb)
            arg = tb;
        if (arg > te)
            arg = te;
        /*  search for knot interval t(l) <= arg < t(l+1) */
        while (!(arg < t[l1 - 1] || l == nk1)) {
            l = l1;
            l1 = l + 1;
        }
        /*  evaluate the non-zero b-splines of degree k-nu at arg. */
        fpbspl(t, kk, arg, l, h);
        /*  find the value of the derivative at x=arg. */
        sp = 0.;
        ll = l - k1;
        for (j = 0; j < k2; ++j) {
            ++ll;
            sp += wrk[ll - 1] * h[j];
        }
        y[i] = sp;
    }
    delete [] wrk;
    return;
}

void FitPack::splev(const double* t, int n, const double* c, int k,
                    const double* x, double* y, int m, int* ier)
{
    /*  subroutine splev evaluates in a number of points x(i),i=1,2,...,m
     a spline s(x) of degree k, given in its b-spline representation.

     calling sequence:
     call splev(t,n,c,k,x,y,m,ier)

     input parameters:
     t    : array,length n, which contains the position of the knots.
     n    : integer, giving the total number of knots of s(x).
     c    : array,length n, which contains the b-spline coefficients.
     k    : integer, giving the degree of s(x).
     x    : array,length m, which contains the points where s(x) must
     be evaluated.
     m    : integer, giving the number of points where s(x) must be
     evaluated.

     output parameter:
     y    : array,length m, giving the value of s(x) at the different
     points.
     ier  : error flag
     ier = 0 : normal return
     ier =10 : invalid input data (see restrictions)

     restrictions:
     m >= 1
     t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1.

     other subroutines required: fpbspl.

     references :
     de boor c  : on calculating with b-splines, j. approximation theory
     6 (1972) 50-62.
     cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
     applics 10 (1972) 134-149.
     dierckx p. : curve and surface fitting with splines, monographs on
     numerical analysis, oxford university press, 1993.

     author :
     p.dierckx
     dept. computer science, k.u.leuven
     celestijnenlaan 200a, b-3001 heverlee, belgium.
     e-mail : Paul.Dierckx@cs.kuleuven.ac.be

     latest update : march 1987 */
    /*  C++ translation by Werner Stille. */
    int i, j, k1, l, ll, nk1;
    double arg, sp, tb, te;
    double h[6];
    *ier = 10;
    if (!m)
        return;
    if (m > 1)
        for (i = 1; i < m; ++i)
            if (x[i] < x[i - 1])
                return;
    *ier = 0;
    /*  fetch tb and te, the boundaries of the approximation interval. */
    k1 = k + 1;
    nk1 = n - k1;
    tb = t[k];
    te = t[nk1];
    l = k1;
    /*  main loop for the different points. */
    for (i = 0; i < m; ++i) {
        /*  fetch a new x-value arg. */
        arg = x[i];
        if (arg < tb)
            arg = tb;
        if (arg > te)
            arg = te;
        /*  search for knot interval t(l) <= arg < t(l+1) */
        while (!(arg < t[l] || l == nk1))
            ++l;
        /*  evaluate the non-zero b-splines at arg. */
        fpbspl(t, k, arg, l, h);
        /*  find the value of s(x) at x=arg. */
        sp = 0.0;
        ll = l - k1;
        for (j = 0; j < k1; ++j)
            sp += c[ll++] * h[j];
        y[i] = sp;
    }
    return;
}

double FitPack::splint(const double* t, int n, const double* c, int k,
                       double a, double b, double* wrk)
{
    /*  function splint calculates the integral of a spline function s(x)
     of degree k, which is given in its normalized b-spline representation

     calling sequence:
     aint = splint(t,n,c,k,a,b,wrk)

     input parameters:
     t    : array,length n,which contains the position of the knots
     of s(x).
     n    : integer, giving the total number of knots of s(x).
     c    : array,length n, containing the b-spline coefficients.
     k    : integer, giving the degree of s(x).
     a,b  : real values, containing the end points of the integration
     interval. s(x) is considered to be identically zero outside
     the interval (t(k+1),t(n-k)).

     output parameter:
     aint : real, containing the integral of s(x) between a and b.
     wrk  : real array, length n.  used as working space
     on output, wrk will contain the integrals of the normalized
     b-splines defined on the set of knots.

     other subroutines required: fpintb.

     references :
     gaffney p.w. : the calculation of indefinite integrals of b-splines
     j. inst. maths applics 17 (1976) 37-41.
     dierckx p. : curve and surface fitting with splines, monographs on
     numerical analysis, oxford university press, 1993.

     author :
     p.dierckx
     dept. computer science, k.u.leuven
     celestijnenlaan 200a, b-3001 heverlee, belgium.
     e-mail : Paul.Dierckx@cs.kuleuven.ac.be

     latest update : march 1987 */
    /*  C++ translation by Werner Stille. */
    int i, nk1;
    double ret_val = 0;
    nk1 = n - k - 1;
    /*  calculate the integrals wrk(i) of the normalized b-splines
     ni,k+1(x), i=1,2,...nk1. */
    fpintb(t, n, wrk, nk1, a, b);
    /*  calculate the integral of s(x). */
    for (i = 0; i < nk1; ++i)
        ret_val += c[i] * wrk[i];
    return ret_val;
}

void FitPack::sproot(double* t, int n, double* c, double* zero, int mest,
                     int* m, int* ier, int nu, double offset)
{
    /*  subroutine sproot finds the zeros of a cubic spline s(x),which is
     given in its normalized b-spline representation.

     calling sequence:
     call sproot(t,n,c,zero,mest,m,ier)

     input parameters:
     t    : real array,length n, containing the knots of s(x).
     n    : integer, containing the number of knots.  n>=8
     c    : real array,length n, containing the b-spline coefficients.
     mest : integer, specifying the dimension of array zero.

     output parameters:
     zero : real array,lenth mest, containing the zeros of s(x).
     m    : integer,giving the number of zeros.
     ier  : error flag:
     ier = 0: normal return.
     ier = 1: the number of zeros exceeds mest.
     ier =10: invalid input data (see restrictions).

     other subroutines required: fpcuro

     restrictions:
     1) n>= 8.
     2) t(4) < t(5) < ... < t(n-4) < t(n-3).
     t(1) <= t(2) <= t(3) <= t(4)
     t(n-3) <= t(n-2) <= t(n-1) <= t(n)

     author :
     p.dierckx
     dept. computer science, k.u.leuven
     celestijnenlaan 200a, b-3001 heverlee, belgium.
     e-mail : Paul.Dierckx@cs.kuleuven.ac.be

     latest update : march 1987 */

    /*  Modified to find the positions of non-zero value of the spline or
     its first or second derivative.

     calling sequence:
     sproot(t, n, c, zero, mest, &m, &ier, nu, offset);

     input parameters:
     nu    : integer, specifying the order of the derivative. 0<=nu<=2
     offset: double, specifying the value of the spline or derivative.

     C++ translation and modifications by Werner Stille. */
    bool z3 = true, nz3 = false;
    bool z, z0, z1, z2, z4, nz0, nz1, nz2, nz4;
    int i, j, l, j1, n4;
    double p0 = 0.0, p1 = 0.0;
    double a0, a1, a2, a3, b0, b1, c1, c2, c3, c4, c5, d4, d5, h1, h2, p2, p3,
    t1, t2, t3, t4, t5, ah, bh, zz;
    double y[3];
    /*  before starting computations a data check is made. if the input data
     are invalid, control is immediately repassed to the calling program. */
    n4 = n - 4;
    *ier = 10;
    if (n < 8)
        return;
    if ((nu < 0) || (nu > 2))
        return;
    j = n;
    for (i = 0; i < 3; ++i) {
        if (t[i] > t[i + 1])
            return;
        if (t[j - 1] < t[j - 2])
            return;
        --j;
    }
    for (i = 3; i < n4; ++i)
        if (t[i] >= t[i + 1])
            return;
    /*  the problem considered reduces to finding the zeros of the cubic
     polynomials pl(x) which define the cubic spline in each knot
     interval t(l)<=x<=t(l+1). a zero of pl(x) is also a zero of s(x) on
     the condition that it belongs to the knot interval.
     the cubic polynomial pl(x) is determined by computing s(t(l)),
     s'(t(l)),s(t(l+1)) and s'(t(l+1)). in fact we only have to compute
     s(t(l+1)) and s'(t(l+1)); because of the continuity conditions of
     splines and their derivatives, the value of s(t(l)) and s'(t(l))
     is already known from the foregoing knot interval. */
    *ier = 0;
    /*  evaluate some constants for the first knot interval */
    h1 = t[3] - t[2];
    h2 = t[4] - t[3];
    t1 = t[3] - t[1];
    t2 = t[4] - t[2];
    t3 = t[5] - t[3];
    t4 = t[4] - t[1];
    t5 = t[5] - t[2];
    /*  calculate a0 = s(t(4)) and ah = s'(t(4)). */
    c1 = c[0];
    c2 = c[1];
    c3 = c[2];
    c4 = (c2 - c1) / t4;
    c5 = (c3 - c2) / t5;
    d4 = (h2 * c1 + t1 * c2) / t4;
    d5 = (t3 * c2 + h1 * c3) / t5;
    a0 = (h2 * d4 + h1 * d5) / t2;
    ah = 3 * (h2 * c4 + h1 * c5) / t2;
    z1 = (ah >= 0);
    nz1 = !z1;
    *m = 0;
    /*  main loop for the different knot intervals. */
    for (l = 3; l < n4; ++l) {
        /*  evaluate some constants for the knot interval t(l) <= x <= t(l+1). */
        h1 = h2;
        h2 = t[l + 2] - t[l + 1];
        t1 = t2;
        t2 = t3;
        t3 = t[l + 3] - t[l + 1];
        t4 = t5;
        t5 = t[l + 3] - t[l];
        /*  find a0 = s(t(l)), ah = s'(t(l)), b0 = s(t(l+1)) and bh = s'(t(l+1)). */
        c1 = c2;
        c2 = c3;
        c3 = c[l];
        c4 = c5;
        c5 = (c3 - c2) / t5;
        d4 = (h2 * c1 + t1 * c2) / t4;
        d5 = (h1 * c3 + t3 * c2) / t5;
        b0 = (h2 * d4 + h1 * d5) / t2;
        bh = 3 * (h2 * c4 + h1 * c5) / t2;
        /*  calculate the coefficients a0,a1,a2 and a3 of the cubic polynomial
         pl(x) = ql(y) = a0+a1*y+a2*y**2+a3*y**3 ; y = (x-t(l))/(t(l+1)-t(l)). */
        a1 = ah * h1;
        b1 = bh * h1;
        a2 = 3 * (b0 - a0) - b1 - 2 * a1;
        a3 = 2 * (a0 - b0) + b1 + a1;
        /*  test whether or not pl(x) could have a zero in the range
         t(l) <= x <= t(l+1). */
        z = true;
        p2 = p3 = 0;
        switch (nu) {
            case 0:
                p0 = a0 - offset;
                p1 = a1;
                p2 = a2;
                p3 = a3;
                z3 = (b1 >= 0);
                nz3 = !z3;
                z = (p0 * (b0 - offset) <= 0);
                break;
            case 1:
                p0 = a1 - offset * h1;
                p1 = a2 + a2;
                p2 = 3 * a3;
                break;
            case 2:
                p0 = a2 + a2 - offset * h1 * h1;
                p1 = 6 * a3;
        }
        if (!z) {
            z0 = (p0 >= 0);
            nz0 = !z0;
            z2 = (p2 >= 0);
            nz2 = !z2;
            z4 = (p3 * 3 + p2 >= 0);
            nz4 = !z4;
            z = ((z0 && ((nz1 && (z3 || (z2 && nz4))) || (nz2 && z3 && z4))) || ((nz0 &&
                                                                                  z1 && (nz3 || (nz2 && z4))) || (z2 && nz3 && nz4)));
        }
        if (z) {
            /*  find the zeros of ql(y). */
            fpcuro(p3, p2, p1, p0, y, &j);
            if (j)
            /*  find which zeros of pl(x) are zeros of s(x). */
                for (i = 0; i < j; ++i) {
                    if (!(y[i] < 0 || y[i] > 1)) {
                        /*  test whether the number of zeros of s(x) exceeds mest. */
                        if (*m >= mest) {
                            *ier = 1;
                            return;
                        }
                        ++(*m);
                        zero[*m - 1] = t[l] + h1 * y[i];
                    }
                }
        }
        a0 = b0;
        ah = bh;
        z1 = z3;
        nz1 = nz3;
    }
    /*  the zeros of s(x) are arranged in increasing order. */
    if (*m < 2)
        return;
    for (i = 1; i < *m; ++i) {
        j = i;
        j1 = j - 1;
        while ((j1 != -1) && (zero[j] < zero[j1])) {
            zz = zero[j];
            zero[j] = zero[j1];
            zero[j1] = zz;
            j = j1;
            j1 = j - 1;
        }
    }
    j = *m;
    *m = 1;
    for (i = 1; i < j; ++i) {
        if (zero[i] != zero[*m - 1]) {
            ++(*m);
            zero[*m - 1] = zero[i];
        }
    }
    return;
}

void FitPack::surfit(int iopt, int m, double* x, double* y, const double* z,
                     const double* w, double xb, double xe, double yb,
                     double ye, int kx, int ky, double s, int nxest, int nyest,
                     int nmax, double eps, int* nx, double* tx, int* ny,
                     double* ty, double* c, double* fp, double* wrk1,
                     int lwrk1, int* ier)
{
    /* given the set of data points (x(i),y(i),z(i)) and the set of positive
     numbers w(i),i=1,...,m, subroutine surfit determines a smooth bivar-
     iate spline approximation s(x,y) of degrees kx and ky on the rect-
     angle xb <= x <= xe, yb <= y <= ye.
     if iopt = -1 surfit calculates the weighted least-squares spline
     according to a given set of knots.
     if iopt >= 0 the total numbers nx and ny of these knots and their
     position tx(j),j=1,...,nx and ty(j),j=1,...,ny are chosen automatic-
     ally by the routine. the smoothness of s(x,y) is then achieved by
     minimalizing the discontinuity jumps in the derivatives of s(x,y)
     across the boundaries of the subpanels (tx(i),tx(i+1))*(ty(j),ty(j+1).
     the amounth of smoothness is determined by the condition that f(p) =
     sum ((w(i)*(z(i)-s(x(i),y(i))))**2) be <= s, with s a given non-neg-
     ative constant, called the smoothing factor.
     the fit is given in the b-spline representation (b-spline coefficients
     c((ny-ky-1)*(i-1)+j),i=1,...,nx-kx-1;j=1,...,ny-ky-1) and can be eval-
     uated by means of subroutine bispev.

     calling sequence:
     call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
     *  nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,iwrk,kwrk,ier)

     parameters:
     iopt  : integer flag. on entry iopt must specify whether a weighted
     least-squares spline (iopt=-1) or a smoothing spline (iopt=0
     or 1) must be determined.
     if iopt=0 the routine will start with an initial set of knots
     tx(i)=xb,tx(i+kx+1)=xe,i=1,...,kx+1;ty(i)=yb,ty(i+ky+1)=ye,i=
     1,...,ky+1. if iopt=1 the routine will continue with the set
     of knots found at the last call of the routine.
     attention: a call with iopt=1 must always be immediately pre-
     ceded by another call with iopt=1 or iopt=0.
     unchanged on exit.
     m     : integer. on entry m must specify the number of data points.
     m >= (kx+1)*(ky+1). unchanged on exit.
     x     : real array of dimension at least (m).
     y     : real array of dimension at least (m).
     z     : real array of dimension at least (m).
     before entry, x(i),y(i),z(i) must be set to the co-ordinates
     of the i-th data point, for i=1,...,m. the order of the data
     points is immaterial. unchanged on exit.
     w     : real array of dimension at least (m). before entry, w(i) must
     be set to the i-th value in the set of weights. the w(i) must
     be strictly positive. unchanged on exit.
     xb,xe : real values. on entry xb,xe,yb and ye must specify the bound-
     yb,ye   aries of the rectangular approximation domain.
     xb<=x(i)<=xe,yb<=y(i)<=ye,i=1,...,m. unchanged on exit.
     kx,ky : integer values. on entry kx and ky must specify the degrees
     of the spline. 1<=kx,ky<=5. it is recommended to use bicubic
     (kx=ky=3) splines. unchanged on exit.
     s     : real. on entry (in case iopt>=0) s must specify the smoothing
     factor. s >=0. unchanged on exit.
     for advice on the choice of s see further comments
     nxest : integer. unchanged on exit.
     nyest : integer. unchanged on exit.
     on entry, nxest and nyest must specify an upper bound for the
     number of knots required in the x- and y-directions respect.
     these numbers will also determine the storage space needed by
     the routine. nxest >= 2*(kx+1), nyest >= 2*(ky+1).
     in most practical situation nxest = kx+1+sqrt(m/2), nyest =
     ky+1+sqrt(m/2) will be sufficient. see also further comments.
     nmax  : integer. on entry nmax must specify the actual dimension of
     the arrays tx and ty. nmax >= nxest, nmax >=nyest.
     unchanged on exit.
     eps   : real.
     on entry, eps must specify a threshold for determining the
     effective rank of an over-determined linear system of equat-
     ions. 0 < eps < 1.  if the number of decimal digits in the
     computer representation of a real number is q, then 10**(-q)
     is a suitable value for eps in most practical applications.
     unchanged on exit.
     nx    : integer.
     unless ier=10 (in case iopt >=0), nx will contain the total
     number of knots with respect to the x-variable, of the spline
     approximation returned. if the computation mode iopt=1 is
     used, the value of nx should be left unchanged between sub-
     sequent calls.
     in case iopt=-1, the value of nx should be specified on entry
     tx    : real array of dimension nmax.
     on succesful exit, this array will contain the knots of the
     spline with respect to the x-variable, i.e. the position of
     the interior knots tx(kx+2),...,tx(nx-kx-1) as well as the
     position of the additional knots tx(1)=...=tx(kx+1)=xb and
     tx(nx-kx)=...=tx(nx)=xe needed for the b-spline representat.
     if the computation mode iopt=1 is used, the values of tx(1),
     ...,tx(nx) should be left unchanged between subsequent calls.
     if the computation mode iopt=-1 is used, the values tx(kx+2),
     ...tx(nx-kx-1) must be supplied by the user, before entry.
     see also the restrictions (ier=10).
     ny    : integer.
     unless ier=10 (in case iopt >=0), ny will contain the total
     number of knots with respect to the y-variable, of the spline
     approximation returned. if the computation mode iopt=1 is
     used, the value of ny should be left unchanged between sub-
     sequent calls.
     in case iopt=-1, the value of ny should be specified on entry
     ty    : real array of dimension nmax.
     on succesful exit, this array will contain the knots of the
     spline with respect to the y-variable, i.e. the position of
     the interior knots ty(ky+2),...,ty(ny-ky-1) as well as the
     position of the additional knots ty(1)=...=ty(ky+1)=yb and
     ty(ny-ky)=...=ty(ny)=ye needed for the b-spline representat.
     if the computation mode iopt=1 is used, the values of ty(1),
     ...,ty(ny) should be left unchanged between subsequent calls.
     if the computation mode iopt=-1 is used, the values ty(ky+2),
     ...ty(ny-ky-1) must be supplied by the user, before entry.
     see also the restrictions (ier=10).
     c     : real array of dimension at least (nxest-kx-1)*(nyest-ky-1).
     on succesful exit, c contains the coefficients of the spline
     approximation s(x,y)
     fp    : real. unless ier=10, fp contains the weighted sum of
     squared residuals of the spline approximation returned.
     wrk1  : real array of dimension (lwrk1). used as workspace.
     if the computation mode iopt=1 is used the value of wrk1(1)
     should be left unchanged between subsequent calls.
     on exit wrk1(2),wrk1(3),...,wrk1(1+(nx-kx-1)*(ny-ky-1)) will
     contain the values d(i)/max(d(i)),i=1,...,(nx-kx-1)*(ny-ky-1)
     with d(i) the i-th diagonal element of the reduced triangular
     matrix for calculating the b-spline coefficients. it includes
     those elements whose square is less than eps,which are treat-
     ed as 0 in the case of presumed rank deficiency (ier<-2).
     lwrk1 : integer. on entry lwrk1 must specify the actual dimension of
     the array wrk1 as declared in the calling (sub)program.
     lwrk1 must not be too small. let
     u = nxest-kx-1, v = nyest-ky-1, km = max(kx,ky)+1,
     ne = max(nxest,nyest), bx = kx*v+ky+1, by = ky*u+kx+1,
     if(bx.le.by) b1 = bx, b2 = b1+v-ky
     if(bx.gt.by) b1 = by, b2 = b1+u-kx  then
     lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
     iwrk  : integer array of dimension (kwrk). used as workspace.
     kwrk  : integer. on entry kwrk must specify the actual dimension of
     the array iwrk as declared in the calling (sub)program.
     kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1).
     ier   : integer. unless the routine detects an error, ier contains a
     non-positive value on exit, i.e.
     ier=0  : normal return. the spline returned has a residual sum of
     squares fp such that abs(fp-s)/s <= tol with tol a relat-
     ive tolerance set to 0.001 by the program.
     ier=-1 : normal return. the spline returned is an interpolating
     spline (fp=0).
     ier=-2 : normal return. the spline returned is the weighted least-
     squares polynomial of degrees kx and ky. in this extreme
     case fp gives the upper bound for the smoothing factor s.
     ier<-2 : warning. the coefficients of the spline returned have been
     computed as the minimal norm least-squares solution of a
     (numerically) rank deficient system. (-ier) gives the rank.
     especially if the rank deficiency which can be computed as
     (nx-kx-1)*(ny-ky-1)+ier, is large the results may be inac-
     curate. they could also seriously depend on the value of
     eps.
     ier=1  : error. the required storage space exceeds the available
     storage space, as specified by the parameters nxest and
     nyest.
     probably causes : nxest or nyest too small. if these param-
     eters are already large, it may also indicate that s is
     too small
     the approximation returned is the weighted least-squares
     spline according to the current set of knots.
     the parameter fp gives the corresponding weighted sum of
     squared residuals (fp>s).
     ier=2  : error. a theoretically impossible result was found during
     the iteration proces for finding a smoothing spline with
     fp = s. probably causes : s too small or badly chosen eps.
     there is an approximation returned but the corresponding
     weighted sum of squared residuals does not satisfy the
     condition abs(fp-s)/s < tol.
     ier=3  : error. the maximal number of iterations maxit (set to 20
     by the program) allowed for finding a smoothing spline
     with fp=s has been reached. probably causes : s too small
     there is an approximation returned but the corresponding
     weighted sum of squared residuals does not satisfy the
     condition abs(fp-s)/s < tol.
     ier=4  : error. no more knots can be added because the number of
     b-spline coefficients (nx-kx-1)*(ny-ky-1) already exceeds
     the number of data points m.
     probably causes : either s or m too small.
     the approximation returned is the weighted least-squares
     spline according to the current set of knots.
     the parameter fp gives the corresponding weighted sum of
     squared residuals (fp>s).
     ier=5  : error. no more knots can be added because the additional
     knot would (quasi) coincide with an old one.
     probably causes : s too small or too large a weight to an
     inaccurate data point.
     the approximation returned is the weighted least-squares
     spline according to the current set of knots.
     the parameter fp gives the corresponding weighted sum of
     squared residuals (fp>s).
     ier=10 : error. on entry, the input data are controlled on validity
     the following restrictions must be satisfied.
     -1<=iopt<=1, 1<=kx,ky<=5, m>=(kx+1)*(ky+1), nxest>=2*kx+2,
     nyest>=2*ky+2, 0<eps<1, nmax>=nxest, nmax>=nyest,
     xb<=x(i)<=xe, yb<=y(i)<=ye, w(i)>0, i=1,...,m
     lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
     kwrk >= m+(nxest-2*kx-1)*(nyest-2*ky-1)
     if iopt=-1: 2*kx+2<=nx<=nxest
     xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
     2*ky+2<=ny<=nyest
     yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
     if iopt>=0: s>=0
     if one of these conditions is found to be violated,control
     is immediately repassed to the calling program. in that
     case there is no approximation returned.
     
     further comments:
     by means of the parameter s, the user can control the tradeoff
     between closeness of fit and smoothness of fit of the approximation.
     if s is too large, the spline will be too smooth and signal will be
     lost ; if s is too small the spline will pick up too much noise. in
     the extreme cases the program will return an interpolating spline if
     s=0 and the weighted least-squares polynomial (degrees kx,ky)if s is
     very large. between these extremes, a properly chosen s will result
     in a good compromise between closeness of fit and smoothness of fit.
     to decide whether an approximation, corresponding to a certain s is
     satisfactory the user is highly recommended to inspect the fits
     graphically.
     recommended values for s depend on the weights w(i). if these are
     taken as 1/d(i) with d(i) an estimate of the standard deviation of
     z(i), a good s-value should be found in the range (m-sqrt(2*m),m+
     sqrt(2*m)). if nothing is known about the statistical error in z(i)
     each w(i) can be set equal to one and s determined by trial and
     error, taking account of the comments above. the best is then to
     start with a very large value of s ( to determine the least-squares
     polynomial and the corresponding upper bound fp0 for s) and then to
     progressively decrease the value of s ( say by a factor 10 in the
     beginning, i.e. s=fp0/10, fp0/100,...and more carefully as the
     approximation shows more detail) to obtain closer fits.
     to choose s very small is strongly discouraged. this considerably
     increases computation time and memory requirements. it may also
     cause rank-deficiency (ier<-2) and endager numerical stability.
     to economize the search for a good s-value the program provides with
     different modes of computation. at the first call of the routine, or
     whenever he wants to restart with the initial set of knots the user
     must set iopt=0.
     if iopt=1 the program will continue with the set of knots found at
     the last call of the routine. this will save a lot of computation
     time if surfit is called repeatedly for different values of s.
     the number of knots of the spline returned and their location will
     depend on the value of s and on the complexity of the shape of the
     function underlying the data. if the computation mode iopt=1
     is used, the knots returned may also depend on the s-values at
     previous calls (if these were smaller). therefore, if after a number
     of trials with different s-values and iopt=1, the user can finally
     accept a fit as satisfactory, it may be worthwhile for him to call
     surfit once more with the selected value for s but now with iopt=0.
     indeed, surfit may then return an approximation of the same quality
     of fit but with fewer knots and therefore better if data reduction
     is also an important objective for the user.
     the number of knots may also depend on the upper bounds nxest and
     nyest. indeed, if at a certain stage in surfit the number of knots
     in one direction (say nx) has reached the value of its upper bound
     (nxest), then from that moment on all subsequent knots are added
     in the other (y) direction. this may indicate that the value of
     nxest is too small. on the other hand, it gives the user the option
     of limiting the number of knots the routine locates in any direction
     for example, by setting nxest=2*kx+2 (the lowest allowable value for
     nxest), the user can indicate that he wants an approximation which
     is a simple polynomial of degree kx in the variable x.
     
     other subroutines required:
     fpback,fpbspl,fpsurf,fpdisc,fpgivs,fprank,fprati,fprota,fporde
     
     references:
     dierckx p. : an algorithm for surface fitting with spline functions
     ima j. numer. anal. 1 (1981) 267-283.
     dierckx p. : an algorithm for surface fitting with spline functions
     report tw50, dept. computer science,k.u.leuven, 1980.
     dierckx p. : curve and surface fitting with splines, monographs on
     numerical analysis, oxford university press, 1993.
     
     author:
     p.dierckx
     dept. computer science, k.u. leuven
     celestijnenlaan 200a, b-3001 heverlee, belgium.
     e-mail : Paul.Dierckx@cs.kuleuven.ac.be
     
     creation date : may 1979
     latest update : march 1987 */
    /*  C++ translation by Werner Stille. */
    
    /*  we set up the parameters tol and maxit. */
    int i, ib1, ib3, jb1, kmax, km1, km2, kx1, ky1, la, lbx, lby, lco, lf, lff,
    lfp, lh, lq, lsx, lsy, lwest, maxit, ncest, nest, nek, nminx, nminy, nmx,
    nmy, nreg, nrint, nxk, nyk;
    double tol = 0.001;
    maxit = 20;
    /*  before starting computations a data check is made. if the input data
     are invalid,control is immediately repassed to the calling program. */
    *ier = 10;
    if (eps <= 0 || eps >= 1)
        return;
    if (kx <= 0 || kx > 5)
        return;
    kx1 = kx + 1;
    if (ky <= 0 || ky > 5)
        return;
    ky1 = ky + 1;
    kmax = qMax(kx, ky);
    km1 = kmax + 1;
    km2 = km1 + 1;
    if (iopt < -1 || iopt > 1)
        return;
    if (m < kx1 * ky1)
        return;
    nminx = kx1 << 1;
    if (nxest < nminx || nxest > nmax)
        return;
    nminy = ky1 << 1;
    if (nyest < nminy || nyest > nmax)
        return;
    nest = qMax(nxest, nyest);
    nxk = nxest - kx1;
    nyk = nyest - ky1;
    ncest = nxk * nyk;
    nmx = nxest - nminx + 1;
    nmy = nyest - nminy + 1;
    nrint = nmx + nmy;
    nreg = nmx * nmy;
    ib1 = kx * nyk + ky1;
    jb1 = ky * nxk + kx1;
    ib3 = kx1 * nyk + 1;
    if (ib1 > jb1) {
        ib1 = jb1;
        ib3 = ky1 * nxk + 1;
    }
    lwest = ncest * (ib1 + 2 + ib3) + 2 * (nrint + nest * km2 + m * km1) + ib3;
    if (lwrk1 < lwest)
        return;
    if (xb >= xe || yb >= ye)
        return;
    for (i = 0; i < m; ++i) {
        if (w[i] <= 0)
            return;
        if (x[i] < xb || x[i] > xe)
            return;
        if (y[i] < yb || y[i] > ye)
            return;
    }
    if (iopt >= 0) {
        if (s < 0)
            return;
    } else {
        if (*nx < nminx || *nx > nxest)
            return;
        nxk = *nx - kx1;
        tx[kx1 - 1] = xb;
        tx[nxk] = xe;
        for (i = kx1; i <= nxk; ++i)
            if (tx[i] <= tx[i - 1])
                return;
        if (*ny < nminy || *ny > nyest)
            return;
        nyk = *ny - ky1;
        ty[ky1 - 1] = yb;
        ty[nyk] = ye;
        for (i = ky1; i <= nyk; ++i)
            if (ty[i] <= ty[i - 1])
                return;
    }
    *ier = 0;
    /*  we partition the working space and determine the spline approximation */
    lq = 1;
    la = lq + ncest * ib3;
    lf = la + ncest * ib1;
    lff = lf + ncest;
    lfp = lff + ncest;
    lco = lfp + nrint;
    lh = lco + nrint;
    lbx = lh + ib3;
    nek = nest * km2;
    lby = lbx + nek;
    lsx = lby + nek;
    lsy = lsx + m * km1;
    {
        double** a = new double*[ncest];
        double** q = new double*[ncest];
        double** bx = new double*[ncest];
        double** by= new double*[ncest];
        double** sx = new double*[m];
        double** sy= new double*[m];
        for (i = 0; i < ncest; ++i) {
            a[i] = &wrk1[la + i * ib1];
            q[i] = &wrk1[lq + i * ib3];
            bx[i] = &wrk1[lbx + i * km2];
            by[i] = &wrk1[lby + i * km2];
        }
        for (i = 0; i < m; ++i) {
            sx[i] = &wrk1[lsx + i * km1];
            sy[i] = &wrk1[lsy + i * km1];
        }
        fpsurf(iopt, m, x, y, z, w, xb, xe, yb, ye, kx, ky, s, nxest, nyest, eps,
               tol, maxit, nreg, nx, tx, ny, ty, c, fp, wrk1, &wrk1[lfp],
               &wrk1[lco], &wrk1[lf], &wrk1[lff], a, q, bx, by, sx, sy, &wrk1[lh],
               ier);
        delete [] sy;
        delete [] sx;
        delete [] by;
        delete [] bx;
        delete [] q;
        delete [] a;
    }
    return;
}
