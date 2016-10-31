/*
 * =====================================================================================
 *
 *       Filename:  theta1D.c
 *
 *    Description:  Heat Equation solver with theta method in 1D
 *
 *        Version:  1.0
 *        Created:  10/16/2016 14:38:29
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

#define	    L	    3.0         /* Length of muffin	3	inches */
#define	    H	    .5          /* Height of muffin	.5	inches */
#define	    T	    70.0        /* Total time		70     minutes */
#define	    ALPHA   1.0         /* Thermal diffusivity  1   inch^2/min */

#ifndef FILENAME
#define FILENAME "theta2D.txt"
#endif

#ifndef	    DX
#define DX	.01
#endif

#ifndef	    DY
#define DY	.01
#endif

#ifndef	    DT
#define DT	.05
#endif

#ifndef	    THETA
static const double THETA = 1;
#endif

#if defined(MALLORY)
static const double par_cnt[4] = {90.49, -.06361, 80.13, -.001023};	/* parameters for boundary model */
static const double par_bnd[4] = {73.97, -.08249, 81.56, -.001303};	/* parameters for central  model */
  static const double t0 = -10.16
#elif defined(GEETA)
static const double par_cnt[4] = {81.09, -.09036, 92.93, -.002168};	/* parameters for boundary model */
static const double par_bnd[4] = {80.35, -.1156, 93.69, -.002442};	/* parameters for central  model */
  static const double t0 = 0;
#else
static const double par_cnt[4] = {85.6,  -.09233, 90.35, -.002332}; /* parameters for boundary model */
static const double par_bnd[4] = {63.76,  -.2109, 96.21, -.003575}; /* parameters for central  model */
static const double t0 = -1.97;
#endif


double** d2malloc( const int nx, const int ny, const char* s );
double BoundaryModel( const double t );
double CenterModel  ( const double t );
void   CalcPar ( const int nx, const int ny, const double *restrict *restrict u,
                const double kx, const double ky, const double ax, const double ay,
                const double bx, const double by,
                double *restrict *restrict e, double *restrict *restrict f, const char mode );


int main( int argc, char *argv[] ) {
    
    /* Local Declaration */
    
    double  **u;                /* matrix of temperature */
    double  **u1;               /* matrix of temp for next time step */
    
    double  **e, **f;		/* parameters for Thomas Algorithm */
    
    FILE    *fp;                /* file pointer for output */
    int	    it;                 /* loop index */
    int	    ix,iy;              /* loop index */
    
    
    /* Test stability */
    const double mux = DT / pow(DX,2) * ALPHA;
    const double kx  = mux * (1 - THETA);
    const double muy = DT / pow(DY,2) * ALPHA;
    const double ky  = muy * (1 - THETA);
    if ( kx > .5 || ky > .5 || THETA < 0 ) {
        printf("\nNot stable or THETA < 0.\n\n");
        return 1;
    }
    
    
    /* Initialize parameters */
    const double bnd0 = BoundaryModel(t0);
    
    const double ax = mux * THETA, ay = muy * THETA;
    const double bx = 2 * ax + 1,  by = 2 * ay + 1;
    
    
    /* Initialize x, t and u array */
    const int nx = L/DX + 1;
    const int ny = H/DY + 1;
    const int nt = (T-t0)/DT + 1;
    
    
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
    

    for ( iy = 0 ; iy < ny ; iy++ )
	for ( ix = 0 ; ix < nx ; ix++ )
	    u[iy][ix] = bnd0;


    /* Open file for output */
    if ( (fp = fopen(FILENAME,"w")) == NULL ) {
	perror("Error opening file!!\n");
	return -1;
    }


    /* Print initial condition to Theta1D.txt */
    fprintf(fp, "\n%-6.2f(min)", t0);
    for ( ix = 0 ; ix < nx ; ix++ )
	fprintf(fp, "| x = %4.2f ", ix*DX);
    fprintf(fp, "\n");
    for ( iy = 0 ; iy < ny ; iy++ ) {
	fprintf(fp, "y = %-7.3f", iy*DY);
	for ( ix = 0 ; ix < nx ; ix++ )
	    fprintf(fp, "|%9.5f ", u[iy][ix]);
	fprintf(fp, "\n");
    }
    fprintf(fp, "\n");


    /* Solve for u using Thomas Algorithm */
    for ( iy = 0; iy < ny; iy++ )
	e[iy][0] = 0;
    for ( ix = 1 ; ix < nx-1 ; ix++ )
	e[0][ix] = 0;

    for ( it = 1 ; it < nt ; it++ ) {
	fprintf(fp, "%-6.2f(min)", t0+DT*it);
	for ( ix = 0 ; ix < nx ; ix++ )
	    fprintf(fp, "| x = %4.2f ", ix*DX);
	fprintf(fp, "\n");

	/* Calculate e and f */
	for ( iy = 0 ; iy < ny ; iy++ ) {
	    f[iy][0] = (u1[iy][0] = (u1[iy][nx-1] = BoundaryModel(t0+DT*(it-.5))));
	}
	CalcPar(nx, ny, (const double **)u, kx, ky, ax, ay, bx, by, e, f, 'x');

	/* Calculate u for t[it+1/2] */
	for ( iy = 1 ; iy < ny-1 ; iy ++ )
	    for ( ix = nx-2 ; 0 < ix ; ix-- )
		u1[iy][ix] = f[iy][ix] + e[iy][ix] * u1[iy][ix+1];

	/* Copy u1 into u */
	for ( iy = 0 ; iy < ny ; iy++ ) {
	    memcpy(u[iy], u1[iy], nx*sizeof(double));
	}

	for ( ix = 0 ; ix < nx ; ix++ )
	    f[0][ix] = (u1[0][ix] = (u1[ny-1][ix] = BoundaryModel(t0+DT*it)));
	CalcPar(nx, ny, (const double **)u, kx, ky, ax, ay, bx, by, e, f, 'y');

	for ( iy = ny-2 ; 0 < iy ; iy-- )
	    for ( ix = 0 ; ix < nx-1 ; ix++ )
		u1[iy][ix] = f[iy][ix] + e[iy][ix] * u1[iy+1][ix];

	/* Output result */
	for ( iy = 0 ; iy < ny ; iy++ ) {
	    fprintf(fp, "y = %-7.3f", iy*DY);
	    for ( ix = 0 ; ix < nx ; ix++ ) {
		fprintf(fp, "|%9.5f ", u1[iy][ix]);
		u[iy][ix] = u1[iy][ix];
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


double BoundaryModel( const double t ) {
    return par_bnd[0] * exp( par_bnd[1] * t ) + par_bnd[2] * exp( par_bnd[3] * t );
}


double CenterModel ( const double t ) {
    return par_cnt[0] * exp( par_cnt[1] * t ) + par_cnt[2] * exp( par_cnt[3] * t );
}


void   CalcPar ( const int nx, const int ny, const double *restrict *restrict u,
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
    else {
	fprintf(stderr, "mode should be 'x' or 'y'");
    }
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
