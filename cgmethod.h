#ifndef CGMETHOD_H
#define CGMETHOD_H

#include "structs.h"

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
            double (*func)(double,networkinfo info),networkinfo info);
void linmin(double p[],double xi[], int n,double *fret,double (*func)(double [],networkinfo info),networkinfo info);
void dlinmin(double p[], double xi[], int n, double *fret, double (*func)(double [],networkinfo),
             void (*dfunc)(double [], double [],networkinfo),networkinfo info);
double f1dim(double alpha,networkinfo info);
double df1dim(double x,networkinfo info);
double dbrent(double ax, double bx, double cx, double (*f)(double,networkinfo info),
              double (*df)(double,networkinfo info),double tol, double *xmin,networkinfo info);
double brent(double ax,double bx,double cx,double (*f)(double,networkinfo info),double tol,
            double *xmin,networkinfo info);
double *dvector( int size );
void frprmn(double p[],int n, double ftol, int *iter, double *fret,
            double (*func)(double [],networkinfo info), void (*dfunc)(double [],double [],networkinfo info),networkinfo info);






#endif