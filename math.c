/*
 * Power Crust software, by Nina Amenta, Sunghee Choi and Ravi Krishna Kolluri.
 * Copyright (c) 2000 by the University of Texas
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee under the GNU Public License is hereby granted, 
 * provided that this entire notice  is included in all copies of any software 
 * which is or includes a copy or modification of this software and in all copies 
 * of the supporting documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

#include <stdio.h>
#include <math.h>
#include <float.h> /* used for DBL_MAX macro definition */
#include "hull.h" /* sunghee */

void normalize(double a[3])
{
    double t;
  
    t =SQ(a[0])+SQ(a[1])+SQ(a[2]);
    t = sqrt(t);
    a[0]=a[0]/t;
    a[2]=a[2]/t;
    a[1]=a[1]/t;
}

double sqdist(double a[3], double b[3]) 
{
	/* returns the squared distance between a and b */ 
	return SQ(a[0]-b[0])+SQ(a[1]-b[1])+SQ(a[2]-b[2]);
}

void dir_and_dist(double a[3], double b[3], double dir[3], double* dist) {
    int k;

    for (k=0; k<3; k++) dir[k] = b[k] - a[k];
    *dist = sqrt( SQ(dir[0])+SQ(dir[1])+SQ(dir[2]));
    for (k=0; k<3; k++) dir[k] = dir[k] / (*dist);
}



void crossabc(double a[3], double b[3], double c[3], double n[3])
{
    double t;

    n[0] = (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]);
    n[1] = (b[2]-a[2])*(c[0]-a[0]) - (b[0]-a[0])*(c[2]-a[2]);
    n[2] = (a[0]-b[0])*(a[1]-c[1]) - (a[1]-b[1])*(a[0]-c[0]);
    t =SQ(n[0])+SQ(n[1])+SQ(n[2]);
    t = sqrt(t);
    n[0]=n[0]/t;
    n[2]=n[2]/t;
    n[1]=n[1]/t;
    /* normalized */
}

double dotabac(double a[3], double b[3], double c[3])
{
    return (b[0]-a[0])*(c[0]-a[0])+(b[1]-a[1])*(c[1]-a[1])+(b[2]-a[2])*(c[2]-a[2]);
}

double dotabc(double a[3], double b[3], double c[3])
{
    return (b[0]-a[0])*c[0]+(b[1]-a[1])*c[1]+(b[2]-a[2])*c[2];
}

/*****************************************************************************/
/*                                                                           */
/*  tetcircumcenter()   Find the circumcenter of a tetrahedron.              */
/*                                                                           */
/*  The result is returned both in terms of xyz coordinates and xi-eta-zeta  */
/*  coordinates, relative to the tetrahedron's point `a' (that is, `a' is    */
/*  the origin of both coordinate systems).  Hence, the xyz coordinates      */
/*  returned are NOT absolute; one must add the coordinates of `a' to        */
/*  find the absolute coordinates of the circumcircle.  However, this means  */
/*  that the result is frequently more accurate than would be possible if    */
/*  absolute coordinates were returned, due to limited floating-point        */
/*  precision.  In general, the circumradius can be computed much more       */
/*  accurately.                                                              */
/*                                                                           */
/*  The xi-eta-zeta coordinate system is defined in terms of the             */
/*  tetrahedron.  Point `a' is the origin of the coordinate system.          */
/*  The edge `ab' extends one unit along the xi axis.  The edge `ac'         */
/*  extends one unit along the eta axis.  The edge `ad' extends one unit     */
/*  along the zeta axis.  These coordinate values are useful for linear      */
/*  interpolation.                                                           */
/*                                                                           */
/*  If `xi' is NULL on input, the xi-eta-zeta coordinates will not be        */
/*  computed.                                                                */
/*                                                                           */
/*****************************************************************************/

void tetcircumcenter(a, b, c, d, circumcenter, cond)
    double a[3];
    double b[3];
    double c[3];
    double d[3];
    double circumcenter[3];
    double *cond; 
{
    double xba, yba, zba, xca, yca, zca, xda, yda, zda;
    double balength, calength, dalength;
    double xcrosscd, ycrosscd, zcrosscd;
    double xcrossdb, ycrossdb, zcrossdb;
    double xcrossbc, ycrossbc, zcrossbc;
    double denominator;
    double xcirca, ycirca, zcirca;

    /* Use coordinates relative to point `a' of the tetrahedron. */
    xba = b[0] - a[0];
    yba = b[1] - a[1];
    zba = b[2] - a[2];
    xca = c[0] - a[0];
    yca = c[1] - a[1];
    zca = c[2] - a[2];
    xda = d[0] - a[0];
    yda = d[1] - a[1];
    zda = d[2] - a[2];
    /* Squares of lengths of the edges incident to `a'. */
    balength = xba * xba + yba * yba + zba * zba;
    calength = xca * xca + yca * yca + zca * zca;
    dalength = xda * xda + yda * yda + zda * zda;
    /* Cross products of these edges. */
    xcrosscd = yca * zda - yda * zca;
    ycrosscd = zca * xda - zda * xca;
    zcrosscd = xca * yda - xda * yca;
    xcrossdb = yda * zba - yba * zda;
    ycrossdb = zda * xba - zba * xda;
    zcrossdb = xda * yba - xba * yda;
    xcrossbc = yba * zca - yca * zba;
    ycrossbc = zba * xca - zca * xba;
    zcrossbc = xba * yca - xca * yba;

    /* Calculate the denominator of the formulae. */
#ifdef EXACT
    /* Use orient3d() from http://www.cs.cmu.edu/~quake/robust.html     */
        /*   to ensure a correctly signed (and reasonably accurate) result, */
        /*   avoiding any possibility of division by zero.                  */
        *cond = orient3d(b,c,d,a);
        denominator = 0.5 / (*cond);
#else
        /* Take your chances with floating-point roundoff. */
        denominator = 0.5 / (xba * xcrosscd + yba * ycrosscd + zba * zcrosscd);
#endif

        /* Calculate offset (from `a') of circumcenter. */
        xcirca = (balength * xcrosscd + calength * xcrossdb + dalength * xcrossbc) *
            denominator;
        ycirca = (balength * ycrosscd + calength * ycrossdb + dalength * ycrossbc) *
            denominator;
        zcirca = (balength * zcrosscd + calength * zcrossdb + dalength * zcrossbc) *
            denominator;
        circumcenter[0] = xcirca;
        circumcenter[1] = ycirca;
        circumcenter[2] = zcirca;  
}


/*****************************************************************************/
/*                                                                           */
/*  tricircumcenter3d()   Find the circumcenter of a triangle in 3D.         */
/*                                                                           */
/*  The result is returned both in terms of xyz coordinates and xi-eta       */
/*  coordinates, relative to the triangle's point `a' (that is, `a' is       */
/*  the origin of both coordinate systems).  Hence, the xyz coordinates      */
/*  returned are NOT absolute; one must add the coordinates of `a' to        */
/*  find the absolute coordinates of the circumcircle.  However, this means  */
/*  that the result is frequently more accurate than would be possible if    */
/*  absolute coordinates were returned, due to limited floating-point        */
/*  precision.  In general, the circumradius can be computed much more       */
/*  accurately.                                                              */
/*                                                                           */
/*  The xi-eta coordinate system is defined in terms of the triangle.        */
/*  Point `a' is the origin of the coordinate system.  The edge `ab' extends */
/*  one unit along the xi axis.  The edge `ac' extends one unit along the    */
/*  eta axis.  These coordinate values are useful for linear interpolation.  */
/*                                                                           */
/*  If `xi' is NULL on input, the xi-eta coordinates will not be computed.   */
/*                                                                           */
/*****************************************************************************/

void tricircumcenter3d(a, b, c, circumcenter,cond)
    double a[3];
    double b[3];
    double c[3];
    double circumcenter[3];
    double *cond;
{
    double xba, yba, zba, xca, yca, zca;
    double balength, calength;
    double xcrossbc, ycrossbc, zcrossbc;
    double denominator;
    double xcirca, ycirca, zcirca;
    double ta[2],tb[2],tc[2];

    /* Use coordinates relative to point `a' of the triangle. */
    xba = b[0] - a[0];
    yba = b[1] - a[1];
    zba = b[2] - a[2];
    xca = c[0] - a[0];
    yca = c[1] - a[1];
    zca = c[2] - a[2];
    /* Squares of lengths of the edges incident to `a'. */
    balength = xba * xba + yba * yba + zba * zba;
    calength = xca * xca + yca * yca + zca * zca;
  
    /* Cross product of these edges. */
#ifdef EXACT
    /* Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html     */
        /*   to ensure a correctly signed (and reasonably accurate) result, */
        /*   avoiding any possibility of division by zero.                  */
        ta[0]=b[1];ta[1]=b[2];tb[0]=c[1];tb[1]=c[2];tc[0]=a[1];tc[1]=a[2];
        xcrossbc = orient2d(ta,tb,tc);
        ta[0]=b[2];ta[1]=b[0];tb[0]=c[2];tb[1]=c[0];tc[0]=a[2];tc[1]=a[0];
        ycrossbc = orient2d(ta,tb,tc);
        ta[0]=b[0];ta[1]=b[1];tb[0]=c[0];tb[1]=c[1];tc[0]=a[0];tc[1]=a[1];
        zcrossbc = orient2d(ta,tb,tc);
#else
        /* Take your chances with floating-point roundoff. */
        xcrossbc = yba * zca - yca * zba;
        ycrossbc = zba * xca - zca * xba;
        zcrossbc = xba * yca - xca * yba;
#endif

        /* Calculate the denominator of the formulae. */
        denominator = 0.5 / (xcrossbc * xcrossbc + ycrossbc * ycrossbc +
                             zcrossbc * zcrossbc);

        /* Calculate offset (from `a') of circumcenter. */
        xcirca = ((balength * yca - calength * yba) * zcrossbc -
                  (balength * zca - calength * zba) * ycrossbc) * denominator;
        ycirca = ((balength * zca - calength * zba) * xcrossbc -
                  (balength * xca - calength * xba) * zcrossbc) * denominator;
        zcirca = ((balength * xca - calength * xba) * ycrossbc -
                  (balength * yca - calength * yba) * xcrossbc) * denominator;
        circumcenter[0] = xcirca;
        circumcenter[1] = ycirca;
        circumcenter[2] = zcirca;

}

void tetorthocenter(a, b, c, d, orthocenter, cnum)
    double a[4];
    double b[4];
    double c[4];
    double d[4];
    double orthocenter[3];
    double *cnum;
{
    double xba, yba, zba, xca, yca, zca, xda, yda, zda, wba, wca, wda;
    double balength, calength, dalength;
    double xcrosscd, ycrosscd, zcrosscd;
    double xcrossdb, ycrossdb, zcrossdb;
    double xcrossbc, ycrossbc, zcrossbc;
    double denominator;
    double xcirca, ycirca, zcirca;
    double wa,wb,wc,wd;

    wa = a[0]*a[0] + a[1]*a[1] + a[2]*a[2] - a[3];
    wb = b[0]*b[0] + b[1]*b[1] + b[2]*b[2] - b[3];
    wc = c[0]*c[0] + c[1]*c[1] + c[2]*c[2] - c[3];
    wd = d[0]*d[0] + d[1]*d[1] + d[2]*d[2] - d[3];
    /* Use coordinates relative to point `a' of the tetrahedron. */
    xba = b[0] - a[0];
    yba = b[1] - a[1];
    zba = b[2] - a[2];
    wba = wb - wa;
    xca = c[0] - a[0];
    yca = c[1] - a[1];
    zca = c[2] - a[2];
    wca = wc - wa;
    xda = d[0] - a[0];
    yda = d[1] - a[1];  
    zda = d[2] - a[2];
    wda = wd - wa;
  
    /* Squares of lengths of the edges incident to `a'. */
    balength = xba * xba + yba * yba + zba * zba - wba;
    calength = xca * xca + yca * yca + zca * zca - wca;
    dalength = xda * xda + yda * yda + zda * zda - wda;
    /* Cross products of these edges. */
    xcrosscd = yca * zda - yda * zca;
    ycrosscd = zca * xda - zda * xca;
    zcrosscd = xca * yda - xda * yca;
    xcrossdb = yda * zba - yba * zda;
    ycrossdb = zda * xba - zba * xda;
    zcrossdb = xda * yba - xba * yda;
    xcrossbc = yba * zca - yca * zba;
    ycrossbc = zba * xca - zca * xba;
    zcrossbc = xba * yca - xca * yba;
  
    /* Calculate the denominator of the formulae. */
#ifdef EXACT 
    /* Use orient3d() from http://www.cs.cmu.edu/~quake/robust.html     */
        /*   to ensure a correctly signed (and reasonably accurate) result, */
        /*   avoiding any possibility of division by zero.                  */
        *cnum = orient3d(b, c, d, a);
        denominator = 0.5 / (*cnum);
#else
        /* Take your chances with floating-point roundoff. */
        denominator = 0.5 / (xba * xcrosscd + yba * ycrosscd + zba * zcrosscd);

#endif
  
        /* Calculate offset (from `a') of circumcenter. */
        xcirca = (balength * xcrosscd + calength * xcrossdb + dalength * xcrossbc) *
            denominator;
        ycirca = (balength * ycrosscd + calength * ycrossdb + dalength * ycrossbc) *
            denominator;
        zcirca = (balength * zcrosscd + calength * zcrossdb + dalength * zcrossbc) *
            denominator;
        orthocenter[0] = xcirca;
        orthocenter[1] = ycirca;
        orthocenter[2] = zcirca;
} 
  
/*****************************************************************************/
/*                                                                           */
/*  tricircumcenter()   Find the circumcenter of a triangle.                 */
/*                                                                           */
/*  The result is returned both in terms of x-y coordinates and xi-eta       */
/*  coordinates, relative to the triangle's point `a' (that is, `a' is       */
/*  the origin of both coordinate systems).  Hence, the x-y coordinates      */
/*  returned are NOT absolute; one must add the coordinates of `a' to        */
/*  find the absolute coordinates of the circumcircle.  However, this means  */
/*  that the result is frequently more accurate than would be possible if    */
/*  absolute coordinates were returned, due to limited floating-point        */
/*  precision.  In general, the circumradius can be computed much more       */
/*  accurately.                                                              */
/*                                                                           */
/*  The xi-eta coordinate system is defined in terms of the triangle.        */
/*  Point `a' is the origin of the coordinate system.  The edge `ab' extends */
/*  one unit along the xi axis.  The edge `ac' extends one unit along the    */
/*  eta axis.  These coordinate values are useful for linear interpolation.  */
/*                                                                           */
/*****************************************************************************/

void triorthocenter(a, b, c, orthocenter, cnum)
    double a[3];
    double b[3];
    double c[3];
    double orthocenter[2];
    double *cnum;
{
    double xba, yba, wba, xca, yca, wca;
    double balength, calength;
    double denominator;
    double xcirca, ycirca;

    /* Use coordinates relative to point `a' of the triangle. */
    xba = b[0] - a[0];
    yba = b[1] - a[1];
    wba = b[2] - a[2];
    xca = c[0] - a[0];
    yca = c[1] - a[1];
    wca = b[2] - a[2];
    /* Squares of lengths of the edges incident to `a'. */
    balength = xba * xba + yba * yba - wba;
    calength = xca * xca + yca * yca - wca;

    /* Calculate the denominator of the formulae. */
#ifdef EXACT
    /* Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html     */
        /*   to ensure a correctly signed (and reasonably accurate) result, */
        /*   avoiding any possibility of division by zero.                  */
        *cnum = orient2d(b, c, a);
        denominator = 0.5 / (*cnum);
#else
        /* Take your chances with floating-point roundoff. */
        denominator = 0.5 / (xba * yca - yba * xca);
#endif

        /* Calculate offset (from `a') of circumcenter. */
        xcirca = (yca * balength - yba * calength) * denominator;  
        ycirca = (xba * calength - xca * balength) * denominator;  
        orthocenter[0] = xcirca;
        orthocenter[1] = ycirca;

}
