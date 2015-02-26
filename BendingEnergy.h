#ifndef BENDINGENERGY_H
#define BENDINGENERGY_H

#include<iostream>
#include<vector>
#include<algorithm>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include "makemikadonetwork.h"



double Ebending(const vector<vector<int>> &springpairs,
                const vector<spring> &springlist,
                const Eigen::VectorXd &XY,
                const double kappa,
                const double g11,
                const double g12,
                const double g22);
#endif