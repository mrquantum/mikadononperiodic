#ifndef SIMPLEBENDINGGRAD_H
#define SIMPLEBENDINGGRAD_H

#include "structs.h"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"


double StretchEnergy(const std::vector<spring> &springlist,
                    const Eigen::VectorXd &XY,
                    double sheardeformation);

double BendEnergy(const std::vector<std::vector<int>> &springpairs,
                const std::vector<spring> &springlist,
                const Eigen::VectorXd &XYphys,
                double kappa,
                double sheardeformation);


Eigen::VectorXd EffKappa(const std::vector<std::vector<int>> &springpairs,
                  const std::vector<spring> &springlist,
                  const Eigen::VectorXd &XY,
                  double kappa);

Eigen::VectorXd SimpleBendingGrad(const triplet &Springpairs,
                           const std::vector<spring> &springlist,
                           const Eigen::VectorXd &XY,
                           double sheardeformation); //XY in


void HarmonicGradPhys(const std::vector<spring> &springlist,
                         double* xy,
                         double* xi,
                         int size,
                         double sheardeformation);




#endif