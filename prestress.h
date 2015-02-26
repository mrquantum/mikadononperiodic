#ifndef PRESTRESS_H
#define PRESTRESS_H

#include<vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include "minimizers.h"
#include<iostream>





double prestress(const std::vector<spring> &springlist,
                 const std::vector<std::vector<int>> &springpairs,
                 Eigen::VectorXd &XY,
                 double kappa,int Nit);


#endif // MAKEMIKADONETWORK_H