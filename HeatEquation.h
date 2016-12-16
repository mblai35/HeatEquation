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

#ifndef	    HEAT_L
#define	    HEAT_L	    3.0         /* Length of muffin	3	inches */
#endif

#ifndef	    HEAT_H
#define	    HEAT_H	    .5          /* Height of muffin	.5	inches */
#endif

#if defined(MALLORY)
    static const double par_cnt[4] = {90.49, -.06361, 80.13, -.001023};	/* parameters for boundary model */
    static const double par_bnd[4] = {73.97, -.08249, 81.56, -.001303};	/* parameters for central  model */
    static const double t0 = -10.25;
#define	    HEAT_T	    100.0
#elif defined(GEETA)
    static const double par_cnt[4] = {81.09, -.09036, 92.93, -.002168};	/* parameters for boundary model */
    static const double par_bnd[4] = {80.35, -.1156, 93.69, -.002442};	/* parameters for central  model */
    static const double t0 = 0;
#define	    HEAT_T	    80.0
#else
    static const double par_cnt[4] = {85.6,  -.09233, 90.35, -.002332}; /* parameters for boundary model */
    static const double par_bnd[4] = {63.76,  -.2109, 96.21, -.003575}; /* parameters for central  model */
    static const double t0 = -1.97;
#endif

#ifndef	    HEAT_T
#define	    HEAT_T	    70.0        /* Total time		70     minutes */
#endif


int heat_info_write_1D ( int nx, double dx, int nt, double dt );
int heat_info_write_2D ( int nx, double dx, int ny, double dy, int nt, double dt, double t0 );
double BoundaryModel ( const double t );
double CenterModel ( const double t );
void GetPar_1D (double *dx, double *dt, double *theta);
void GetPar_2D (double *dx, double *dy, double *dt, double *theta);
double InitialModel ( const double x, const double bnd0, const double cnt0 );
void   CalcCoef_1D ( const int nx, const double *u, const double mu, const double theta,
       	double * const e, double * const f, const double a, const double b, const double k );
void   CalcCoef_2D ( const int nx, const int ny, const double *__restrict__ *__restrict__ u,
	const double kx, const double ky, const double ax, const double ay,
	const double bx, const double by,
	double *__restrict__ *__restrict__ e, double *__restrict__ *__restrict__ f, const char mode );
