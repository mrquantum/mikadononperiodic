#ifndef EMERGYANDGRADIENTS_H
#define EMERGYANDGRADIENTS_H

#include<vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>

double Energynetwork(const std::vector<spring> &springlist, const Eigen::VectorXd &XY,double k,double L);
Eigen::VectorXd Gradient(const std::vector<spring> &springlist,const Eigen::VectorXd &XY,double k, double L);
double dEda(const Eigen::VectorXd &XY,const Eigen::VectorXd &s0,const std::vector<spring> &springlist,double k,double L); 



#endif // MAKEMIKADONETWORK_H