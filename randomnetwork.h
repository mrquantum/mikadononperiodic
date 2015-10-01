#ifndef RANDOMNETWORK_H
#define RANDOMNETWORK_H

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
#include "makeSpringNetworks.h"

using namespace std;

Eigen::VectorXd initiateRandomNetwork(vector<spring> &springlist,
                                vector<vector<int>> &springpairs,
                                vector<stick> &mikado,
                                vector<stick> &mikorig,
                                vector<elonstick> &ELONSTICK,
                                vector<connected> &Connection,
                                vector<node> &nodes,
                                vector<node> &singleNodes,
                                vector<vector<int>> &conmatr,
                                vector<vector<int>> &Clusterv,
                                vector<vector<int>> &numberdistribution,
                                int NumberMikado,
                                vector<spring> &background,
                                Eigen::VectorXd &XYb,
                                int SEED,
                                double LStick,
                                double rlenshort,
                                double rlenlong,
                                double k1,
                                double k2,
                                double stretchf,
                                ofstream &springfile,
                                ofstream &anglefile,
                                ofstream &mikadofile,
                                ofstream &clusterdistribution,
                                ofstream &cluster,
                                ofstream &clusterdata,
                                ofstream &nodefile                                
                               );


#endif