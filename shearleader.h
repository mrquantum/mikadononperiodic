//In this header file we put all the shear protocol
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
//#include "exportfiles.h"
#include "writefunctions.h"
#include "cgwoagony.h"

using namespace std;

void shearsteps(double deltaboxdx,
               int NumberStepsRight,
               int NumberStepsLeft,
	       vector<spring> &springlist,
	       vector<vector<int>> &springpairs,
	       Eigen::VectorXd &XY,
	       int bendingOn,
	       double kappa,
	       int Nit,
	       double tolGradE,
	       ofstream &shearcoordinates,
	       ofstream &shearenergy);