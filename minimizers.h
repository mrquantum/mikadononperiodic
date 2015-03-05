#ifndef MINIMIZERS_H
#define MINIMIZERS_H

#include<vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include "BendingGrad.h"

int doBracketfind(double &a1,
                  double &a2,
                   const Eigen::VectorXd &XY,
                   const Eigen::VectorXd &s0,
                   const std::vector<spring> &springlist,
                   const std::vector<std::vector<int>> &springpairs, 
                   double kappa,
                   const double g11,
                   const double g12,
                   const double g22);


void doFalsePosition(double &xl,double &xh,double &root,
                 const Eigen::VectorXd &XY,
                 const Eigen::VectorXd &s0,
                 const std::vector<spring> &springlist,
                 const std::vector<std::vector<int>> &springpairs, 
                 double kappa,
                 const double g11,
                 const double g12,
                 const double g22);


void doSecant(   double &root,
                 const Eigen::VectorXd &XY,
                 const Eigen::VectorXd &s0,
                 const std::vector<spring> &springlist,
                 const std::vector<std::vector<int>> &springpairs, 
                 double kappa,
                 const double g11,
                 const double g12,
                 const double g22);


void doConjStep(Eigen::VectorXd &XY,
                Eigen::VectorXd &s0,
                Eigen::VectorXd &gradE,
                const std::vector<spring> &springlist,
                const std::vector<std::vector<int>> &springpairs,
                double kappa,
                int conjstep,
                double g11,
                double g12,
                double g22);



#endif // MINIMIZERS_H