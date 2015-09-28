#ifndef CWOAGONY_H
#define CWOAGONY_H


#include <iostream>
#include <ctime>
#include "random.h"
#include "makemikadonetwork.h"
#include "EnergyandGradients.h"
#include "minimizers.h"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <math.h>
#include <functional>
#include "importparam.h"
#include "BendingEnergy.h"
#include "BendingGrad.h"
#include "clusters.h"
#include <stdio.h>
#include <stdlib.h>
#include "writefunctions.h"
#include "shearleader.h"



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