/*
 * =====================================================================================
 *
 *       Filename:  HeatEquation.c
 *
 *    Description:  For heat equation solvers
 *
 *        Version:  1.0
 *        Created:  12/04/2016 19:33:41
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xiukun Hu (xhu), xiukun.hu@outlook.com
 *   Organization:  University of Wyoming, Math Dept.
 *
 * =====================================================================================
 */

#include "HeatEquation.h"

char file_info_1D[] = { "theta1D.txt" };
char file_info_2D[] = { "theta2D.txt" };


int heat_info_write_1D ( int nx, double dx, int nt, double dt ) {
    FILE * fp;
    if ( (fp = fopen(file_info_1D,"w")) == NULL){
	fprintf( stderr, "heat_info_write_1D: Could not open %s\n", file_info_1D );
	return 1;
    }

    rewind( fp );
    fprintf( fp, "%d %e\n", nx, dx );
    fprintf( fp, "%d %e\n", nt, dt );

    fclose( fp );
    return 0;
}


int heat_info_write_2D ( int nx, double dx, int ny, double dy, int nt, double dt, double t0 ) {
    FILE * fp;
    if ( (fp = fopen(file_info_2D,"w")) == NULL){
	fprintf( stderr, "heat_info_write_2D: Could not open %s\n", file_info_2D );
	return 1;
    }

    rewind( fp );
    fprintf( fp, "%d %e\n", nx, dx );
    fprintf( fp, "%d %e\n", ny, dy );
    fprintf( fp, "%d %e %e\n", nt, dt, t0 );

    fclose( fp );
    return 0;
    
}

double BoundaryModel ( const double t ) {
    return par_bnd[0] * exp( par_bnd[1] * t ) + par_bnd[2] * exp( par_bnd[3] * t );
}


double CenterModel ( const double t ) {
    return par_cnt[0] * exp( par_cnt[1] * t ) + par_cnt[2] * exp( par_cnt[3] * t );
}


void GetPar_1D (double *dx, double *dt, double *theta) {
	/* Input dx and make sure it can devide total length */
	while (1) {
	    printf("Please input dx: ");
	    scanf("%lf", dx);
	    fflush(stdin);
	    if ( remainder(HEAT_L, *dx) < 1e-6 ) {
		printf("dx = %lf\n", *dx);
		break;
	    }
	    else printf("Total length cannot be devided by dx!\n\n");
	}

	/* Input dt and make sure it can devide total time */
	while (1) {
	    printf("Please input dt: ");
	    scanf("%lf", dt);
	    fflush(stdin);
	    if ( remainder(HEAT_T, *dt) < 1e-8 ) {
		printf("dt = %lf\n", *dt);
		break;
	    }
	    else printf("Total time cannot be devided by dt!\n\n");
	}

	/* Input theta and make sure 0 <= theta <= 1 */
	while (1){
	    printf("Please input theta (type '-1' to make theta 0.5-∆x^2/(12∆t)): ");
	    scanf("%lf", theta);
	    if ( 0 <= *theta && *theta <= 1 ) {
		printf("theta = %lf\n", *theta);
		break;
	    }
	    else if ( *theta == -1 ) {
		*theta = .5 - *dx * *dx / (12.0 * *dt);
		printf("theta = %lf\n", *theta);
		break;
	    }
	    else printf("theta must be in [0,1] or -1!\n\n");
	}
}


void GetPar_2D (double *dx, double *dy, double *dt, double *theta) {
    /* Input dx and make sure it can devide total length */
    while (1) {
	printf("Please input dy: ");
	scanf("%lf", dy);
	fflush(stdin);
	if ( remainder(HEAT_H, *dy) < 1e-6 ) {
	    printf("dy = %lf\n", *dy);
	    break;
	}
	else printf("\nTotal height cannot be devided by dy!\n\n");
    }
    GetPar_1D (dx, dt, theta);

}


double InitialModel ( const double x, const double bnd0, const double cnt0 ) {
    return (bnd0 - cnt0) * 4.0 / (HEAT_L*HEAT_L) 
	* (x - HEAT_L/2.0) * (x - HEAT_L/2.0) + cnt0;
}


void   CalcCoef_1D ( const int nx, const double *u, const double mu, const double theta,
	double * const e, double * const f, const double a, const double b, const double k ) {
    int i;
    double d;
    for ( i = 1 ; i < nx-1 ; i++ ) {
	d = u[i-1] * k + u[i] * (1 - 2*k) + u[i+1] * k;
	e[i] = a / (b - a * e[i-1]);
	f[i] = (d + a * f[i-1]) / (b - a * e[i-1]);
    }
}


void   CalcCoef_2D ( const int nx, const int ny, const double *restrict *restrict u,
	const double kx, const double ky, const double ax, const double ay,
	const double bx, const double by,
	double *restrict *restrict e, double *restrict *restrict f, const char mode ) {
    int ix, iy;
    double d;
    if ( mode == 'x' ) 
	for ( iy = 1 ; iy < ny-1 ; iy++ ) {
	    for ( ix = 1 ; ix < nx -1  ; ix++ ) {
		d = u[iy-1][ix] * ky + u[iy][ix] * (1 - 2*ky) + u[iy+1][ix]*ky;
		e[iy][ix] = ax / (bx - ax * e[iy][ix-1]);
		f[iy][ix] = (d + ax * f[iy][ix-1]) / (bx - ax * e[iy][ix-1]);
	    }
	}
    else if ( mode == 'y' ) 
	for ( iy = 1 ; iy < ny-1 ; iy++ ) 
	    for ( ix = 1 ; ix < nx-1 ; ix++ ) {
		d = u[iy][ix-1] * kx + u[iy][ix] * (1 - 2*kx) + u[iy][ix+1]*kx;
		e[iy][ix] = ay / (by - ay * e[iy-1][ix]);
		f[iy][ix] = (d + ay * f[iy-1][ix]) / (by - ay * e[iy-1][ix]);
	    }
}
