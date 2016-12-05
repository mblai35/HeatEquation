/*
 * =====================================================================================
 *
 *       Filename:  HeatEquation.h
 *
 *    Description:  For heat equation solvers
 *
 *        Version:  1.0
 *        Created:  12/04/2016 19:47:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xiukun Hu (xhu), xiukun.hu@outlook.com
 *   Organization:  University of Wyoming, Math Dept.
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef	    L
#define	    L	    3.0         /* Length of muffin	3	inches */
#endif

#ifndef	    H
#define	    H	    .5          /* Height of muffin	.5	inches */
#endif

#ifndef	    T
#define	    T	    70.0        /* Total time		70     minutes */
#endif

#ifndef	    ALPHA
#define	    ALPHA   1.0         /* Thermal diffusivity  1   inch^2/min */
#endif

#if defined(MALLORY)
    static const double par_cnt[4] = {90.49, -.06361, 80.13, -.001023};	/* parameters for boundary model */
    static const double par_bnd[4] = {73.97, -.08249, 81.56, -.001303};	/* parameters for central  model */
    static const double t0 = -10.25;
#elif defined(GEETA)
    static const double par_cnt[4] = {81.09, -.09036, 92.93, -.002168};	/* parameters for boundary model */
    static const double par_bnd[4] = {80.35, -.1156, 93.69, -.002442};	/* parameters for central  model */
    static const double t0 = 0;
#else
    static const double par_cnt[4] = {85.6,  -.09233, 90.35, -.002332}; /* parameters for boundary model */
    static const double par_bnd[4] = {63.76,  -.2109, 96.21, -.003575}; /* parameters for central  model */
    static const double t0 = -1.97;
#endif

double BoundaryModel( const double t );
double CenterModel ( const double t );
void GetPar_1D (double *dx, double *dt, double *theta);
void GetPar_2D (double *dx, double *dy, double *dt, double *theta);
double InitialModel ( const double x, const double bnd0, const double cnt0 );
void   CalcCoef_1D ( const int nx, const double *u, const double mu, const double theta,
       	double * const e, double * const f, const double a, const double b, const double k );
void   CalcCoef_2D ( const int nx, const int ny, const double *restrict *restrict u,
	const double kx, const double ky, const double ax, const double ay,
	const double bx, const double by,
	double *restrict *restrict e, double *restrict *restrict f, const char mode );
