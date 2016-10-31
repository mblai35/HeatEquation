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

#ifndef	    DX
    #define DX	.01
#endif

#ifndef	    DY
    #define DY	.01
#endif

#ifndef	    DT	
    #define DT	.00005
#endif

#ifndef	    THETA
    static const double THETA = .5 - DX * DX / (12.0 * DT);
#endif

#if defined(MALLORY)
    static const double par_cnt[4] = {90.49, -.06361, 80.13, -.001023};	/* parameters for boundary model */
    static const double par_bnd[4] = {73.97, -.08249, 81.56, -.001303};	/* parameters for central  model */
//  static const double t0 =
#elif defined(GEETA)
    static const double par_cnt[4] = {81.09, -.09036, 92.93, -.002168};	/* parameters for boundary model */
    static const double par_bnd[4] = {80.35, -.1156, 93.69, -.002442};	/* parameters for central  model */
//  static const double t0 =
#else
    static const double par_cnt[4] = {85.6,  -.09233, 90.35, -.002332}; /* parameters for boundary model */
    static const double par_bnd[4] = {63.76,  -.2109, 96.21, -.003575}; /* parameters for central  model */
    static const double t0 = -1.97;
#endif


double** d2malloc( const int nx, const int ny, const char* s );
double BoundaryModel( const double t );
double CenterModel  ( const double t );
void   CalcPar	    ( const int nx, const double **u, const double mu, const double THETA, double **e, double **f );
void   CalcU1	    ( const int nx, const double **e, const double **f, double *u1 );


int main( int argc, char *argv[] ) {

    /* Local Declaration */

    double  *x;                 /* array of x    points */
    double  *y;                 /* array of y	 points */
    double  *t;                 /* array of time points */
    double  **u;                /* matrix of temperature */
    double  **u1;               /* matrix of temp for next time step */

    double  **e, **f;		/* parameters for Thomas Algorithm */

    FILE    *fp;                /* file pointer for output */
    int	    i,j;                  /* loop index */
    int	    it;                 /* loop index */
    int	    ix,iy;              /* loop index */


    /* Check for arguments */
    if ( argc != 1 ) {
	if ( argc != 6 && argc != 11 ) {
	    perror("Wrong number of arguments.");
	    return 1;
	}
	i = 1;
	while ( argc != i ) {
	    if ( strcmp(argv[i], "c") ) 
		for ( i++ ; i < 6 ; i++ )
		    sscanf(argv[i], "%lf", &par_cnt[i]);
	    else if ( strcmp(argv[i], "b") ) 
		for ( i++ ; i < 6 ; i++ ) 
		    sscanf(argv[i], "%lf", &par_bnd[i]);
	    else {
		printf("Illegal argunent %d, must be 'c' or 'b'\n", i);
		return 2;
	    }
	}
    }


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

    if ( (t = (double *) malloc( nt * sizeof(double) )) == NULL ) {
	return -2;
    }

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

    t[0] = t0;
    for ( i = 1 ; i < nt ; i++ ) {
	t[i] = t[0] + DT;
    }

    for ( ix = 0 ; ix < nx ; ix++ ) 
	for ( iy = 0 ; iy < ny ; iy++ ) 
	    u[ix][iy] = bnd0;


    /* Open file for output */
    if ( (fp = fopen("Theta2D.txt","w+")) == NULL ) {
	perror("Error opening file\n");
	return -1;
    }


    /* Print initial condition to Theta1D.txt */
    fprintf(fp, "Time	   ");
    for ( ix = 0 ; ix < nx ; ix++ ) 
	fprintf(fp, "| x = %4.2f ", ix*DX);
    fprintf(fp, "\n%-6.3f(min)", t[0]);
    for ( iy = 0 ; iy < ny ; iy++ ) { 
	for ( ix = 0 ; ix < nx ; ix++ ) 
	    fprintf(fp, "|%9.5f ", u[ix]);
	fprintf(fp, "\n");
    }
    fprintf(fp, "\n");


    /* Solve for u using Thomas Algorithm */
    for ( it = 1 ; it < nt ; it++ ) {
	fprintf(fp, "%-6.3f(min)", t[it]);

	/* Calculate e and f */
	e[0] = 0;
	f[0] = (u1[0] = (u1[nx-1] = BoundaryModel(t[it])));
	CalcPar(nx, u, mu, THETA, e, f);

	/* Calculate u for t[it] */
	for ( ix = nx-2 ; 0 < ix ; ix-- ) 
	    u1[ix] = f[ix] + e[ix] * u1[ix+1];

	/* Output result */
	for ( ix = 0 ; ix < nx ; ix++ ) {
	    fprintf(fp, "|%9.5f ", u1[ix]);
	    u[ix] = u1[ix];
	}
	fprintf(fp, "\n");
    }

    fclose(fp);
    free(x);
    free(t);
    free(e);
    free(f);
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


void   CalcPar ( const int nx, const double **u, const double mu, const double THETA, double **e, double **f, const char mode ) {
    int i;
    double d;
    for ( i = 1 ; i < nx-1 ; i++ ) {
	d = u[i-1] * k + u[i] * (1 - 2*k) + u[i+1] * k;
	e[i] = a / (b - a * e[i-1]);
	f[i] = (d + a * f[i-1]) / (b - a * e[i-1]);
    }
}


double** d2malloc( const int nx, const int ny, const char* s ){
    int	     i,j;
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
