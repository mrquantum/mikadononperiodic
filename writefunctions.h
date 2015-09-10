#ifndef WRITEFUNCTIONS_H
#define WRITEFUNCTIONS_H


#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"
#include <fstream>
#include "clusters.h"
#include "makemikadonetwork.h"
using namespace std;

void Write_Mikado_2txt(ofstream &mikadofile,vector<stick> &mikado);
void Write_Springs_2txt(ofstream &springfile,vector<spring> &springlist);
void Write_Clusterdistribution_2txt(ofstream &clusterdistribution,vector<vector<int>> &numberdistribution);
void Write_Clusterdata_2txt(ofstream &clusterdata,int SEED,int NumberMikado,double LStick,vector<vector<int>> &C);
void Write_Nodes_2txt(ofstream &nodefile,vector<node> &singleNodes);
void Write_ShearCoordinates_2txt(ofstream &shearcoordinates,Eigen::VectorXd &XY);
void Write_ShearEnergy_2txt(ofstream &shearenergy,double boxdx,double ETOT,double ESTRETCH, double EBEND, double lenGrad, int conjsteps);
void Write_Angles_2txt(ofstream &anglefile,int node1,int node2,int node3);
void Write_IndividualClusters_2txt(ofstream &cluster,vector<vector<int>> Clusterv);

#endif