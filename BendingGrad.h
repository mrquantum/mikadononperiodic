#ifndef BENDINGGRAD_H
#define BENDINGGRAD_H

#include<iostream>
#include<vector>
#include<algorithm>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
//#include "makemikadonetwork.h"
#include "structs.h"


using namespace std;

Eigen::VectorXd BendingGrad(const vector<vector<int>> &springpairs,
                     const vector<spring> &springlist,
                     const Eigen::VectorXd &XY,
                     double kappa,
                     double g11,
                     double g12,
                     double g22);

void BendinggradN(double *p,double *xi,networkinfo params);


double l_ij(int springnr, Eigen::VectorXd &XY,vector<spring> &springlist,double g11, double g12,double g22);

double L_ij(double x1,double y1,double x2,double y2);
vector<gradstruct> gradL(int one, int two,double x1,double y1,double x2,double y2,double Lij);
double thsq(double x1,double y1, double x2, double y2, double x3, double y3, double L12, double L23);

vector<gradstruct> bendgradpart1(double f, double L12, double L23,vector<gradstruct> &gradL12,vector<gradstruct> &gradL23);

vector<gradstruct> bendgradpart2(double L12, double L23, 
                             double x1,double y1, double x2, double y2, double x3, double y3,
                             int i1,int i2,int i3,
                             vector<gradstruct> &gradL12, vector<gradstruct> &gradL23);

void physbendinggradient(double *p, double *g,networkinfo &info);



#endif