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



double l_ij(int springnr, Eigen::VectorXd &XY,vector<spring> &springlist,double g11, double g12,double g22);

double L_ij(double x1,double y1,double x2,double y2);
vector<gradstruct> gradL(int one, int two,double x1,double y1,double x2,double y2,double Lij);
double thsq(double x1,double y1, double x2, double y2, double x3, double y3, double L12, double L23);
void physbendinggradient(double *p, double *g,networkinfo &info);



#endif