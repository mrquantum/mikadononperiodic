#ifndef MAKESPRINGNETWORKS_H
#define MAKESPRINGNETWORKS_H

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
#include "structs.h"

Eigen::VectorXd makeSquareNetwork(int Number, vector<spring> &springlist);
Eigen::VectorXd triangularnetwork(int Number, vector<spring> &springlist);
spring makespring(int one, int two,Eigen::VectorXd &X, Eigen::VectorXd &Y,int wlr,int wud);




#endif
