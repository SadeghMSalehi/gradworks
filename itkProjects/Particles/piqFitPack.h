/***************************************************************************
                          fitpack.h  -  description
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

#ifndef FITPACK_H
#define FITPACK_H

/**
  * Spline class. Provides functions for curve and surface fitting with
  * splines. The functions base on the FORTRAN program collection FITPACK by
  * P. Diercks. C++ translation by Werner Stille. The mathematical fundamentals
  * of FITPACK are described in the book <pre>
  Curve and Surface Fitting with Splines, by P. Dierckx,
  Monographs on Numerical Analysis, Oxford University Press, 1993,
  (ISBN 0-19-853441-8, 286 pages, 56 line figures) </pre>
  *
  */

class FitPack {
public:
  FitPack();
  ~FitPack();
  static void bispev(const double* tx, int nx, const double* ty, int ny,
                     const double* c, int kx, int ky, const double* x, int mx,
                     const double* y, int my, double *z, int* ier);
  static void curfit(int iopt, int m, const double* x, const double* y,
                     const double* w, double xb, double xe, int k, double s,
                     int nest, int* n, double* t, double* c, double* fp,
                     double* wrk, int lwrk, int* iwrk, int* ier);
  static void splder(const double* t, int n, const double* c, int k, int nu,
                     const double* x, double* y, int m, int* ier);
  static void splev(const double* t, int n, const double* c, int k,
                    const double* x, double* y, int m, int *ier);
  static double splint(const double* t, int n, const double* c, int k,
                       double a, double b, double* wrk);
  static void sproot(double* t, int n, double* c, double* zero, int mest,
                     int* m, int* ier, int nu = 0, double offset = 0.0);
  static void surfit(int iopt, int m, double* x, double* y, const double* z,
                     const double* w, double xb, double xe, double yb,
                     double ye, int kx, int ky, double s, int nxest, int nyest,
                     int nmax, double eps, int* nx, double* tx, int* ny,
                     double* ty, double* c, double* fp, double* wrk1,
                     int lwrk1, int* ier);
};

#endif
