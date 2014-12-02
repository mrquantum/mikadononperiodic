#ifndef EMERGYANDGRADIENTS_H
#define EMERGYANDGRADIENTS_H

#include<vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>

struct anchor{
  int label; 
  double xpos,ypos;
};

int SIGN(double a,double b);
int sgn(double x);
double ROSENBROCK(const Eigen::VectorXd &XY);
Eigen::VectorXd GRAD_rosen(const Eigen::VectorXd &XY);
double dROSENdA(const Eigen::VectorXd &XY,const Eigen::VectorXd &s);


double distance(const Eigen::VectorXd &XY,int one, int two);
double dldxi(const Eigen::VectorXd &XY,int one, int two);
double dldyi(const Eigen::VectorXd &XY,int one,int two);
Eigen::VectorXd gradL(const double x1,const double y1,const double x2, const double y2, 
               const int springnr, const std::vector<spring> &springlist,const int num);


double distance1(const double x1, const double y1, const double x2,const double y2);

double Energynetwork(const std::vector<spring> &springlist, 
                     const Eigen::VectorXd &XY);

double Ebend(const std::vector<std::vector<int>> &springpairs,
             const std::vector<spring> &springlist,
             const Eigen::VectorXd &XY,
             const double kappa);

Eigen::VectorXd Gradient(const std::vector<spring> &springlist,
                         const Eigen::VectorXd &XY);

Eigen::VectorXd gradEbend(const std::vector<std::vector<int>> &springpairs,
                          const std::vector<spring> &springlist,
                          const Eigen::VectorXd &XY,
                          double kappa);

double dEda(const Eigen::VectorXd &XY,
            const Eigen::VectorXd &s0,
            const std::vector<spring> &springlist,
            const std::vector<std::vector<int>> &springpairs,
            double kappa); //,const std::vector<triplet> &tripl,const int BendOn); 

double quad(double x);









#endif // MAKEMIKADONETWORK_H