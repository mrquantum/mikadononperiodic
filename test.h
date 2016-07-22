#ifndef TEST_H
#define TEST_H

#include "structs.h"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"
#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>
using namespace std;


Eigen::VectorXd XYtest()
{
    Eigen::VectorXd XY(10);
    XY<<0.5,0.5,0.5,0.2,0.8,0.2,0.5,0.8,0.5,0.5;
    return XY;
}

vector<spring> testsprings()
{
    vector<spring> springlist(0);
    spring newspring;
    newspring.one=0;
    newspring.two=1;
    newspring.wlr=0;
    newspring.wud=0;
    newspring.sticki=0;
    newspring.k=1;
    newspring.rlen=0.3;
    springlist.push_back(newspring);
    
    newspring.one=1;
    newspring.two=2;
    newspring.wlr=0;
    newspring.wud=0;
    newspring.sticki=0;
    newspring.k=1;
    newspring.rlen=0.3;
    springlist.push_back(newspring);
    
    newspring.one=2;
    newspring.two=0;
    newspring.wlr=0;
    newspring.wud=1;
    newspring.sticki=0;
    newspring.k=1;
    newspring.rlen=0.4;
    springlist.push_back(newspring);
    
    newspring.one=3;
    newspring.two=1;
    newspring.wlr=0;
    newspring.wud=0;
    newspring.sticki=1;
    newspring.k=1;
    newspring.rlen=0.3;
    springlist.push_back(newspring);
    
    newspring.one=1;
    newspring.two=4;
    newspring.wlr=0;
    newspring.wud=0;
    newspring.sticki=1;
    newspring.k=1;
    newspring.rlen=0.3;
    springlist.push_back(newspring);
    
    newspring.one=4;
    newspring.two=3;
    newspring.wlr=1;
    newspring.wud=0;
    newspring.sticki=1;
    newspring.k=1;
    newspring.rlen=0.4;
    springlist.push_back(newspring);
    
    newspring.one=3;
    newspring.two=2;
    newspring.wlr=0;
    newspring.wud=0;
    newspring.sticki=2;
    newspring.k=1;
    newspring.rlen=sqrt(2*pow(0.3,2));
    springlist.push_back(newspring);
    
    return springlist;
}

vector<vector<int       >       > testpairs()
{
    vector<vector<int>> pairs(0);
    vector<int> pair1(2);
    pair1[0]=3;
    pair1[1]=1;
    pairs.push_back(pair1);
    pair1[0]=1;
    pair1[1]=4;
    pairs.push_back(pair1);
    
    return pairs;
}












#endif
