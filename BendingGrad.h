#ifndef BENDINGGRAD_H
#define BENDINGGRAD_H

#include<iostream>
#include<vector>
#include<algorithm>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "makemikadonetwork.h"

using namespace std;

Eigen::VectorXd BendingGrad(const vector<vector<int>> &springpairs,
                     const vector<spring> &springlist,
                     const Eigen::VectorXd &XY,
                     double kappa,
                     double g11,
                     double g12,
                     double g22);



double l_ij(int springnr, Eigen::VectorXd &XY,vector<spring> &springlist,double g11, double g12,double g22);


#endif