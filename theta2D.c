/*
 * =====================================================================================
 *
 *       Filename:  theta2D.c
 *
 *    Description:  Heat Equation solver with theta method in 2D
 *
 *        Version:  1.0
 *        Created:  10/16/2016 14:38:29
 *       Revision:  12/04/2016 19:26:41
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
#include "HeatEquation.h"

double** d2malloc( const int nx, const int ny, const char* s );

int main( int argc, char *argv[] ) {

    /* Local Declaration */
    double  dx;                 /* ∆x */
    double  dy;                 /* ∆y */
    double  dt;                 /* ∆t */

    double  mux, muy;           /* µ and k */
    double  theta;              /* theta for theta method */

    double  bnd0;               /* Initial boundary value and center value */
    double  ax, ay, bx, by,
	    kx, ky;             /* Parameters for Thomas Algorithm */

    double  **u;                /* matrix of temperature */
    double  **u1;               /* matrix of temp for next time step */
    double  **temppb;           /* temporary 2D double pointer */

    double  **e, **f;		/* parameters for Thomas Algorithm */

    FILE    *fp;                /* file pointer for output */
    int	    i,j;                /* loop index */
    int	    it;                 /* loop index */
    int	    ix,iy;              /* loop index */
    int	    nx,ny,nt;


    while(1) {
    /* Get dx, dy, dt and theta */
    GetPar_2D ( &dx, &dy, &dt, &theta);

    /* Test stability */
    mux = dt / pow(dx,2) * ALPHA;
    kx  = mux * (1 - theta);
    muy = dt / pow(dy,2) * ALPHA;
    ky  = muy * (1 - theta);
    if ( kx > .5 || ky > .5 || theta < 0 ) 
	printf("\nNot stable or theta < 0.\n\n");   
    else
	break;
    }


    /* Initialize parameters */
    bnd0 = BoundaryModel(t0);

    ax = mux * theta, ay = muy * theta;
    bx = 2 * ax + 1,  by = 2 * ay + 1; 


    /* Initialize x, t and u array */
    nx = (double)L/dx + 1;
    ny = (double)H/dy + 1;
    nt = ((double)T-t0)/dt + 1;


    if ( (u = d2malloc(nx,ny,"u")) == NULL ) {
	return -2;
    }

    if ( (u1 = d2malloc(nx,ny,"u1")) == NULL ) {
	return -2;
    }

    if ( (e = d2malloc(nx,ny,"e")) == NULL ) {
	return -2;
    }

    if ( (f = d2malloc(nx,ny,"f")) == NULL ) {
	return -2;
    }


    for ( ix = 0 ; ix < nx ; ix++ ) 
	for ( iy = 0 ; iy < ny ; iy++ ) 
	    u[iy][ix] = bnd0;


    /* Open file for output */
    if ( (fp = fopen("Theta2D.csv","w+")) == NULL ) {
	perror("Error opening file\n");
	return -1;
    }


    /* Print initial condition to Theta1D.txt */
    fprintf(fp, "Time (min)");
    for ( ix = 0 ; ix < nx ; ix++ ) 
	fprintf(fp, ", x = %4.2f ", ix*dx);
    fprintf(fp, "\n%-6.3f     ", t0);
    for ( iy = 0 ; iy < ny ; iy++ ) { 
	for ( ix = 0 ; ix < nx ; ix++ ) 
	    fprintf(fp, ",%9.5f ", u[iy][ix]);
	fprintf(fp, "\n");
    }
    fprintf(fp, "\n");


    /* Solve for u using Thomas Algorithm */
    for ( iy = 0; iy < ny; iy++ )
	    e[iy][0] = 0;
    for ( ix = 1 ; ix < nx-1 ; ix++ ) 
	    e[0][ix] = 0;
    for ( it = 1 ; it < nt ; it++ ) {
	fprintf(fp, "%-6.3f     ", t0+dt*it);

	/* Calculate e and f */
	for ( iy = 0 ; iy < ny ; iy++ ) {
	    f[iy][0] = (u1[iy][0] = (u1[iy][nx-1] = BoundaryModel(t0+dt*(it-.5))));
	}
	CalcCoef_2D(nx, ny, (const double **)u, kx, ky, ax, ay, bx, by, e, f, 'x');

	/* Calculate u for t[it+1/2] */
	for ( iy = 1 ; iy < ny-1 ; iy ++ ) 
	    for ( ix = nx-2 ; 0 < ix ; ix-- ) 
		u1[iy][ix] = f[iy][ix] + e[iy][ix] * u1[iy][ix+1];

	/* Swap u and u1 */
	temppb = u1;
	u1 = u;
	u = temppb;

	for ( ix = 0 ; ix < nx ; ix++ ) 
	    f[0][ix] = (u1[0][ix] = (u1[ny-1][ix] = BoundaryModel(t0+dt*it)));
	CalcCoef_2D(nx, ny, (const double **)u, kx, ky, ax, ay, bx, by, e, f, 'y');

	for ( iy = ny-2 ; 0 < iy ; iy-- ) 
	    for ( ix = 0 ; ix < nx-1 ; ix++ ) 
		u1[iy][ix] = f[iy][ix] + e[iy][ix] * u1[iy][ix-1];

	/* Output result */
	for ( iy = 0 ; iy < ny ; iy++ ) {
	    for ( ix = 0 ; ix < nx ; ix++ ) {
		fprintf(fp, ",%9.5f ", u1[iy][ix]);
	    }
	    fprintf(fp,"\n");
	}
	fprintf(fp, "\n");
    }

    fclose(fp);
    free(e);
    free(f);
    for ( iy = 0 ; iy < ny ; iy++ ) {
	free(u[iy]);
	free(u1[iy]);
    }
    free(u);
    free(u1);
    return 0;
}


double** d2malloc( const int nx, const int ny, const char* s ){
    int	    i,j;
    double** pb;
    if ( (pb = (double**) malloc( ny * sizeof(double *) )) == NULL ) {
	fprintf(stderr,"memory allocation for %s" ,s);
	return NULL;
    }

    for ( i = 0 ; i < ny ; i++ ) {
	if ( (pb[i] = (double*) malloc( nx * sizeof(double) )) == NULL ) {
	    fprintf(stderr,"memory allocation for %s[%d]",s,i);
	    for ( j = i-1 ; j >=0 ; j-- ) {
		free(pb[j]);
	    }
	    free(pb);
	    return NULL;
	}
    }
    return pb;
}
