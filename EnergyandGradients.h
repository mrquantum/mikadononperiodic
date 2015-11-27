#ifndef EMERGYANDGRADIENTS_H
#define EMERGYANDGRADIENTS_H

#include<vector>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
//#include "makemikadonetwork.h"
#include "structs.h"


struct anchor{
  int label; 
  double xpos,ypos;
};

double Dist(double x1,double y1,double x2,double y2,
            double g11,double g12,double g22);

double Energynetwork(const std::vector<spring> &springlist, 
                     const Eigen::VectorXd &XY,
                     const double g11,
                     const double g12,
                     const double g22);


double EnergyNetworkn(double *XY,networkinfo parameters);
void HarmonicGradientn(double *p,double *xi,networkinfo params);

double Ebend(const std::vector<std::vector<int>> &springpairs,
             const std::vector<spring> &springlist,
             const Eigen::VectorXd &XY,
             const double g11,
             const double g12,
             const double g22,
             const double kappa);

Eigen::VectorXd HarmonicGradient(const std::vector<spring> &springlist,
                        const Eigen::VectorXd &XY,
                        const double g11,
                        const double g12,
                        const double g22);

Eigen::VectorXd gradEbend(const std::vector<std::vector<int>> &springpairs, 
                    const std::vector<spring> &springlist, 
                    const Eigen::VectorXd &XY,
                    double g11,double g12, double g22, double kappa);

double dEda(const Eigen::VectorXd &XY,
            const Eigen::VectorXd &s0,
            const std::vector<spring> &springlist,
            const std::vector<std::vector<int>> &springpairs,
            double kappa,
            const double g11,
            const double g12,
            const double g22); //,const std::vector<triplet> &tripl,const int BendOn); 










#endif // MAKEMIKADONETWORK_H