#ifndef STRUCTS_H
#define STRUCTS_H
#include <vector>
#include <array>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"



struct stick{
  double x,y,th,length;
  int nr, wlr, wud;
};

struct connected{
  int first,second,nrCon,recur;
  double s1,s2;
  int type; //0 for mik mik connections 1 for mik-backgroundconnections
  int backgroundspring[2];
  
};


struct elonstick{
  int sticki;
  //vector<numberpos> NumberAndPos;
  std::vector<int> nr;
  std::vector<double> S;
  std::vector<int> type;
  std::vector<std::array<int,2>> backgroundspringn;
 
  //include what spring it is. node 1 node j1
};
struct spring{
  int one, two, wlr, wud, sticki;
  double rlen,k;
};
struct node{
  //int number, sticki;
  int number;
  //int stick1;
  //int stick2;
  double x;
  double y;
};

class networkinfo {
public:
    networkinfo( std::vector<spring> &springlist_in, std::vector<std::vector<int>> &springpairs_in): springlist(springlist_in), springpairs(springpairs_in){};
    std::vector<spring>& springlist;
    std::vector<std::vector<int>>& springpairs;
    double g11,g12,g22;
    double sheardeformation;
    double kappa;
    int bendingon;
    int size;
//     Eigen::VectorXd effkappa;
};

struct triplet{
    std::vector<std::vector<int>> springpairs;
    Eigen::VectorXd EffKappa; 
};

struct gradstruct{//this struct denotes the partiel derivative of any function wrth the ith variable. 
    double ddxi;
    double ddyi;
    int i;
};

struct stresstensor{
    double sxx;
    double syy;
    double sxy;
};


#endif