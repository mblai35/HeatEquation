/*
 * =====================================================================================
 *
 *       Filename:  theta1D.c
 *
 *    Description:  Heat Equation solver for all theta in 1D. When theta equals zero,
 *		    explicit solver will be used. Otherwise, Thomas algorithm is applied.
 *		    The output will be a csv file named #theta1D.csv#.
 *
 *        Version:  1.0
 *        Created:  10/16/2016 14:38:29
 *       Revision:  12/04/2016 17:36:50
 *       Compiler:  gcc
 *
 *         Author:  Xiukun Hu, Mallory Lai, Geeta Monpara
 *   Organization:  University of Wyoming.
 *
 * =====================================================================================
 */


#include "HeatEquation.h"
#ifndef	    HEAT_ALPHA
#define	    HEAT_ALPHA	    .2	/* Thermal diffusivity in inch/min */
#endif

int ExplicitSolution( const double dx, const double dt, const double bnd0, const double cnt0);


int main( int argc, char *argv[] ) {

    /* Local Declaration */
    double  dx;                 /* ∆x */
    double  dt;                 /* ∆t */

    double  mu;                 /* µ = alpha * ∆t/(∆x^2) */
    double  theta;              /* theta for theta method */

    double  bnd0, cnt0;         /* Initial boundary value and center value */
    double  a,b,k;              /* Parameters for Thomas Algorithm */

    double  *x;                 /* array of space points */
    double  *t;                 /* array of time  points */
    double  *u;                 /* array of solution */
    double  *u1;                /* array of solution for next time step */

    double  *e, *f;		/* parameters for Thomas Algorithm */

    FILE    *fp;                /* file pointer for output */

    int	    nx;                 /* number of space points */
    int	    nt;                 /* number of time  points */

    int	    i;                  /* loop index */
    int	    it;                 /* loop index */
    int	    ix;                 /* loop index */


    /* Input dx, dt and theta and ensure stability */
    while (1) {

	/* Get dx, dt and theta */
	GetPar_1D ( &dx, &dt, &theta );

	/* Test stability */
	mu = dt / (dx*dx) * HEAT_ALPHA;
	k  = mu * (1 - theta);
	if ( k > .5 || theta < 0 ) printf("\nNot stable or theta < 0.\n\n");
	else break;
    }


    /* Initialize parameters */
    bnd0 = BoundaryModel(0);
    cnt0 = CenterModel(0);

    a = mu * theta;
    b = 2 * a + 1; 


    /* If theta == 0, then use explicit solver */
    if ( theta == 0 ) {
	i = ExplicitSolution( dx, dt, bnd0, cnt0 );
	return i;
    }

    /* Initialize x, t and u array */
    nx = (double)HEAT_L/dx + 1;
    nt = (double)HEAT_T/dt + 1;

    if ( (x = (double *) malloc( nx * sizeof(double) )) == NULL ) {
	perror("memory allocation for x");
	return -2;
    }
    if ( (u = (double *) malloc( nx * sizeof(double) )) == NULL ) {
	perror("memory allocation for u");
	return -2;
    }
    if ( (u1= (double *) malloc( nx * sizeof(double) )) == NULL ) {
	perror("memory allocation for u1");
	return -2;
    }
    if ( (t = (double *) malloc( nt * sizeof(double) )) == NULL ) {
	perror("memory allocation for t");
	return -2;
    }
    if ( (e = (double *) malloc( nx * sizeof(double) )) == NULL ) {
	perror("memory allocation for u");
	return -2;
    }
    if ( (f = (double *) malloc( nx * sizeof(double) )) == NULL ) {
	perror("memory allocation for u");
	return -2;
    }

    for ( i = 0 ; i < nx ; i++ ) {
	x[i] = i * dx;
	u[i] = InitialModel( x[i], bnd0, cnt0 );
    }

    for ( i = 0 ; i < nt ; i++ ) {
	t[i] = i * dt;
    }


    /* Write info */
    heat_info_write_1D( nx, dx, nt, dt );


    /* Open file for output */
    if ( (fp = fopen("theta1D.bin","w")) == NULL ) {
	perror("Error opening file\n");
	return -1;
    }


    /* Print initial condition to Theta1D.bin */
    rewind( fp );
    fwrite( u, sizeof(double), nx, fp );


    /* Solve for u using Thomas Algorithm */
    for ( it = 1 ; it < nt ; it++ ) {

	/* Calculate e and f */
	e[0] = 0;
	f[0] = (u1[0] = (u1[nx-1] = BoundaryModel(t[it])));
	CalcCoef_1D(nx, u, mu, theta, e, f, a, b, k);

	/* Calculate u for t[it] */
	for ( ix = nx-2 ; 0 < ix ; ix-- ) 
	    u1[ix] = f[ix] + e[ix] * u1[ix+1];

	/* u1 -> u */
	for ( ix = 0 ; ix < nx ; ix++ ) 
	    u[ix] = u1[ix];

	/* Write u */
	fwrite( u, sizeof(double), nx, fp );
	    
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


int ExplicitSolution( const double dx, const double dt, const double bnd0, const double cnt0 ) {

    // Numerical parameters
    const int rows = ((double)HEAT_T/dt) + 1; //number of rows in grid
    const int cols = ((double)HEAT_L/dx) + 1; //number of columns in grid

    // Counter variables
    int i;
    int j;

    // Create temperature grid with dimensions matching the number of columns and rows
    double ** TemperatureGrid;
    if ( (TemperatureGrid = (double **) malloc( rows * sizeof(double *) )) == NULL ) {
	perror("memory allocation for TemperatureGrid");
	return -1;
    }

     for ( i = 0 ; i < rows ; i++ ) {
	 if ( (TemperatureGrid[i] = (double *) malloc( cols * sizeof(double ) )) == NULL ) {
	     perror("memory allocation for TemperatureGrid");
	     for ( j = i ; j >= 0 ; j-- ) 
		 free(TemperatureGrid[j]);
	     free(TemperatureGrid);
	     return -1;
	 }
     }

    // Find the appropriate row increments to plug into exponential function for boundary conditions.
    double rowIncrements[rows];


    // Row increments
    for(i = 0; i < rows; i++){
	rowIncrements[i] = i * dt;
    }

    // Initialize boundary conditions of grid
    // Boundary; Populate first column of TemperatureGrid using exponential line of best fit.
    for(i = 0; i < rows; i++){
	TemperatureGrid[i][cols-1] = TemperatureGrid[i][0] = BoundaryModel(rowIncrements[i]);
    }

    // Top boundary; Linear interpolation of first row.
    for(i = 1; i < cols; i++){
	TemperatureGrid[0][i] = InitialModel(i*dx, bnd0, cnt0);
    }


    // Calculation; explicit solution.
    for(i = 1; i < rows; i++){
	for(j = 1; j < (cols - 1); j++){
	    TemperatureGrid[i][j] = TemperatureGrid[i-1][j] + 
		HEAT_ALPHA*dt*((TemperatureGrid[i-1][j-1] - 2*TemperatureGrid[i-1][j] + TemperatureGrid[i-1][j+1])
			/(dx*dx));
	}
    }


    // Write TemperatureGrid to csv file.
    heat_info_write_1D( cols, dx, rows, dt );
    FILE * grid=fopen("theta1D.bin", "w");
    for ( i = 0 ; i < rows ; i++ ) 
	fwrite( TemperatureGrid[i], sizeof(double), cols, grid );


    free(TemperatureGrid);
    fclose(grid);


    return 0;
}
