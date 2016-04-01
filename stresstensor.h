#ifndef STRESSTENSOR_H
#define STRESSTENSOR_H

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"

stresstensor StressTensor(std::vector<spring> &springlist, Eigen::VectorXd &XY,double strain);





#endif
