#ifndef MINIMIZERS_H
#define MINIMIZERS_H

#include<vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>





int doBracketfind(double &a1,double &a2,
                   const Eigen::VectorXd &XY,
                   const Eigen::VectorXd &s0,
                   const std::vector<spring> &springlist,
                   const std::vector<std::vector<int>> &springpairs, 
                   double kappa,
                   const double g11,
                   const double g12,
                   const double g22);

void doBisection(double a1,double a2,double &root,
                 const Eigen::VectorXd &XY,
                 const Eigen::VectorXd &s0,
                 const std::vector<spring> &springlist,
                 const std::vector<std::vector<int>> &springpairs, 
                 double kappa,
                 const double g11,
                 const double g12,
                 const double g22
                );


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
                std::vector<spring> &springlist,
                std::vector<std::vector<int>> &springpairs,
                double &root,
                double kappa,
                int conjstep,
                const double g11,
                const double g12,
                const double g22);

void doSteepestDescent(Eigen::VectorXd &XY,
                Eigen::VectorXd &s0,
                Eigen::VectorXd &gradE,
                std::vector<spring> &springlist,
                std::vector<std::vector<int>> &springpairs,
                double &root,
                double kappa,
                const double g11,
                const double g12,
                const double g22);
                //Eigen::VectorXd &b);

Eigen::VectorXd Hessianapprox(const Eigen::VectorXd &XY,
                              const Eigen::VectorXd &XYm1,
                              const Eigen::VectorXd &g0,
                              const Eigen::VectorXd &g0m1);







#endif // MINIMIZERS_H