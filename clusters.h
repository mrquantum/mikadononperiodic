#ifndef CLUSTERS_H
#define CLUSTERS_H




#include <iostream>
#include <ctime>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <math.h>
#include "makemikadonetwork.h"
using namespace std;

struct site{int number; 
    int beenthere;};


vector<vector<int>> connectivitymatrix(vector<spring> springs);
vector<vector<int>> connectivitymatrix(vector<vector<int>> connectedSticks,int numbermikado);

void addelement(int element, vector<int> &v);

vector< vector< int> > clusters(vector<vector<int> > M);

vector<vector<int>> Numberdistribution(vector<vector<int>> Clusters,int NSticks);

    
    


#endif

