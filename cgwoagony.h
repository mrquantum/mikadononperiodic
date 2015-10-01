#ifndef CWOAGONY_H
#define CWOAGONY_H

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"
#include <vector>
#include "structs.h"



using namespace std;


void CGAGONY(Eigen::VectorXd &XY,
             vector<spring> &springlist, 
             vector<vector<int>> &springpairs, 
             int bendingon, 
             double kappa, 
             double g11, double g12, double g22);

void CGAGONY2(Eigen::VectorXd &XY,
             vector<spring> &springlist, 
             vector<vector<int>> &springpairs, 
             int bendingon, 
             double kappa, 
             double g11, double g12, double g22);
             
             
             
#endif