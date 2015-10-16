#ifndef MAKEMIKADONETWORK_H
#define MAKEMIKADONETWORK_H

#include <iostream>
#include <ctime>
#include "random.h"
//#include "EnergyandGradients.h"
// #include "minimizers.h"
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
#include "randomnetwork.h"
#include "makeSpringNetworks.h"
#include "structs.h"


using namespace std;

// struct stick{
//   double x,y,th;
//   int nr, wlr, wud;
// };
// 
// struct connected{
//   int first,second,nrCon,recur;
//   double s1,s2;
//   int type; //0 for mik mik connections 1 for mik-backgroundconnections
//   int backgroundspring[2];
//   
// };
// 
// 
// struct elonstick{
//   int sticki;
//   //vector<numberpos> NumberAndPos;
//   std::vector<int> nr;
//   std::vector<double> S;
// };
// struct spring{
//   int one, two, wlr, wud, sticki;
//   double rlen,k;
// };
// struct node{
//   //int number, sticki;
//   int number;
//   //int stick1;
//   //int stick2;
//   double x;
//   double y;
// };




std::vector<stick> make_sticks(int N);
std::vector<stick> make_ghost_lr( const std::vector<stick> &m, double LStick, int NumberMikado);
std::vector<stick> make_ghost_ud(const std::vector<stick> &m, double LStick, int NumberMikado);

void make_connections(std::vector<connected> &Connection,const std::vector<stick> &m, double LStick,
    const vector<spring> &background,
    const Eigen::VectorXd &XYb);

void sortELEMENTSperMIKADO(std::vector<elonstick> &ELONSTICK,std::vector<connected> &Connection);
void orderElonstick(std::vector<int> &order,std::vector<elonstick> &ELONSTICK);

void makeSpringsAndNodes(const std::vector<elonstick> &ELONSTICK,
		     const std::vector<stick> &mikorig, std::vector<spring> &springlist, std::vector<node> &nodes,
		     double rlenshort,double rlenlong,double k1, double k2,double stretchf,vector<spring> &background,
                     Eigen::VectorXd &XYb);
bool operator<(const node& first,const node& second);
bool operator<(const elonstick &first, const elonstick &second);
bool operator<(const stick &first, const stick &second);

void makeSticks(std::vector<stick> &mikado,std::vector<stick> &mikorig,const int NumberMikado,const double LStick);

void makeConnections(std::vector<connected> &Connection,
                     const std::vector<stick> &mikado,
                     const double LStick,
                     const std::vector<spring> &background,
                     const Eigen::VectorXd &XYb);



double inbox(double x,double boxsize);
void makeSpringpairs(std::vector<std::vector<int>> &springpairs,const std::vector<spring> &springlist);
spring makespring(int node1,int node2,double x1,double x2, double y1, double y2,int stick,double k,double stretchf);
void makeanddeletebondsonbackground(vector<spring> &springlist,const vector<elonstick> &ELONSTICK,
                        stick &CURRENTSTICK,vector<double> &posonsticki,Eigen::VectorXd &Xb,Eigen::VectorXd &Yb,
                        int background_size,int sticknr,int elnr);

#endif // MAKEMIKADONETWORK_H