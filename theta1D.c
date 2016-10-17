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

#define	    L	    3.0         /* Length of muffin	3	inches */
#define	    T	    70.0        /* Total time		70     minutes */
#define	    ALPHA   1.0         /* Thermal diffusivity  1   inch^2/min */


double par_cnt[4] = {85.6,  -.09233, 90.35, -.002332}; /* parameters for boundary model */
double par_bnd[4] = {63.76,  -.2109, 96.21, -.003575}; /* parameters for central  model */
double bnd0, cnt0;              /* Initial boundary value and center value */
double a,b,k;                   /* Parameters for Thomas Algorithm */

double BoundaryModel( const double t );
double CenterModel  ( const double t );
double InitialModel ( const double x );
void   CalcPar	    ( const int nx, const double *u, const double mu, const double theta, double *e, double *f );
void   CalcU1	    ( const int nx, const double *e, const double *f, double *u1 );


int main( int argc, char *argv[] ) {

    /* Local Declaration */
    double  dx;                 /* ∆x */
    double  dt;                 /* ∆t */

    double  mu;                 /* µ = alpha * ∆t/(∆x^2) */
    double  theta;              /* theta for theta method */

    double  *x;                 /* array of space points */
    double  *t;                 /* array of time  points */
    double  *u;                 /* array of solution */
    double  *u1;                /* array of solution for next time step */

    double  *e, *f;		/* parameters for Thomas Algorithm */

    FILE    *fp;                 /* file pointer for output */

    int	    nx;                 /* number of space points */
    int	    nt;                 /* number of time  points */

    int	    i;                  /* loop index */
    int	    it;                 /* loop index */
    int	    ix;                 /* loop index */


    /* Check for arguments */
    if ( argc != 1 ) {
	if ( argc != 6 && argc != 11 ) {
	    perror("Wrong number of arguments.");
	    return 1;
	}
	i = 1;
	while ( argc != i ) {
	    if ( argv[i][0] == 'c' ) 
		for ( i++ ; i < 6 ; i++ )
		    sscanf(argv[i], "%lf", &par_cnt[i]);
	    else if ( argv[i][0] == 'b' ) 
		for ( i++ ; i < 6 ; i++ ) 
		    sscanf(argv[i], "%lf", &par_bnd[i]);
	    else {
		perror("Illegal argunent 1 or 5, must be 'c' or 'b'\n");
		return 2;
	    }
	}
    }


    /* Input dx, dt and theta and ensure stability */
    while (1) {

	/* Input dx and make sure it can devide total length */
	while (1) {
	    printf("Please input dx: ");
	    scanf("%lf", &dx);
	    fflush(stdin);
	    if ( remainder(L, dx) < 1e-6 ) {
		printf("\ndx = %lf\n", dx);
		break;
	    }
	    else printf("\nTotal length cannot be devided by dx!");
	}

	/* Input dt and make sure it can devide total time */
	while (1) {
	    printf("\nPlease input dt: ");
	    scanf("%lf", &dt);
	    fflush(stdin);
	    if ( remainder(T, dt) < 1e-8 ) {
		printf("\ndt = %lf\n", dt);
		break;
	    }
	    else printf("\nTotal time cannot be devided by dt!\n");
	}

	/* Input theta and make sure 0 <= theta <= 1 */
	while (1){
	    printf("Please input theta (type '-1' to make theta 0.5-∆x^2/(12∆t)): ");
	    scanf("%lf", &theta);
	    if ( 0 <= theta && theta <= 1 ) {
		printf("\ntheta = %lf\n", theta);
		break;
	    }
	    else if ( theta == -1 ) {
		theta = .5 - dx * dx / (12.0 * dt);
		printf("\ntheta = %lf\n", theta);
		break;
	    }
	    else printf("\ntheta must be in [0,1] or -1!\n");
	}

	/* Test stability */
	mu = dt / pow(dx,2) * ALPHA;
	k  = mu * (1 - theta);
	if ( k > .5 || theta < 0 ) printf("\nNot stable or theta < 0.\n\n");
	else break;
    }


    /* Initialize parameters */
    bnd0 = BoundaryModel(0);
    cnt0 = CenterModel(0);

    a = mu * theta;
    b = 2 * a + 1; 


    /* Initialize x, t and u array */
    nx = L/dx + 1;
    nt = T/dt + 1;

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
	u[i] = InitialModel( x[i] );
    }

    for ( i = 0 ; i < nt ; i++ ) {
	t[i] = i * dt;
    }


    /* Open file for output */
    if ( (fp = fopen("Theta1D.txt","w+")) == NULL ) {
	perror("Error opening file\n");
	return -1;
    }


    /* Print initial condition to Theta1D.txt */
    fprintf(fp, "Time	  ");
    for ( ix = 0 ; ix < nx ; ix++ ) 
	fprintf(fp, "| x = %4.2f ", x[ix]);
    fprintf(fp, "\n%-5.2f(min)", t[0]);
    for ( ix = 0 ; ix < nx ; ix++ ) 
	fprintf(fp, "|%9.5f ", u[ix]);
    fprintf(fp, "\n");


    /* Solve for u using Thomas Algorithm */
    for ( it = 1 ; it < nt ; it++ ) {
	fprintf(fp, "%-5.2f(min)", t[it]);

	/* Calculate e and f */
	e[0] = 0;
	f[0] = (u1[0] = (u1[nx-1] = BoundaryModel(t[it])));
	CalcPar(nx, u, mu, theta, e, f);

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


double InitialModel ( const double x ) {
    return (bnd0 - cnt0) * 4.0 / (L*L) * (x - L/2.0) * (x - L/2.0) + cnt0;
}


void   CalcPar ( const int nx, const double *u, const double mu, const double theta, double *e, double *f ) {
    int i;
    double d;
    for ( i = 1 ; i < nx-1 ; i++ ) {
	d = u[i-1] * k + u[i] * (1 - 2*k) + u[i+1] * k;
	e[i] = a / (b - a * e[i-1]);
	f[i] = (d + a * f[i-1]) / (b - a * e[i-1]);
    }
}
