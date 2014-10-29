#ifndef EMERGYANDGRADIENTS_H
#define EMERGYANDGRADIENTS_H

#include<vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>

struct anchor{
  int label; 
  double xpos,ypos;
};

double distance(const Eigen::VectorXd &XY,int one, int two);
double dldxi(const Eigen::VectorXd &XY,int one, int two);
double dldyi(const Eigen::VectorXd &XY,int one,int two);
Eigen::VectorXd gradL(const Eigen::VectorXd &XY,int one,int two);
double distance1(const double x1, const double y1, const double x2,const double y2);
double Energynetwork(const std::vector<spring> &springlist, const Eigen::VectorXd &XY,const std::vector<anchor> &Anchor);
double Ebend(const std::vector<std::vector<int>> &springpairs,const std::vector<spring> &springlist,const Eigen::VectorXd &XY);

Eigen::VectorXd Gradient(const std::vector<spring> &springlist,const Eigen::VectorXd &XY,
			 const std::vector<anchor> &Anchor);
Eigen::VectorXd gradEbend(const std::vector<triplet> &tripl,const Eigen::VectorXd &XY);

double dEda(const Eigen::VectorXd &XY,const std::vector<anchor> &Anchor,const Eigen::VectorXd &s0,
            const std::vector<spring> &springlist);//,const std::vector<triplet> &tripl,const int BendOn); 


#endif // MAKEMIKADONETWORK_H