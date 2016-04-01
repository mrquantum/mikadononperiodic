#include <iostream>
#include "structs.h"

#include <algorithm>
#include <vector>
#include <math.h>
#include <functional>
#include <stdio.h>
#include <stdlib.h>
//#include "exportfiles.h"
#include "writefunctions.h"
//#include <nlopt.h>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"

using namespace std;
using namespace Eigen;

stresstensor StressTensor(vector<spring> &springlist, VectorXd &XY,double strain){
    double sxx=0;
    double sxy=0;
    double syy=0;
    
    double k,rlen;
    int one,two,wud,wlr;
    int num=XY.size()/2;
    double xb1,yb1,xb2,yb2; //the box coordinates
    double x1,y1,x2,y2; //the real coordinataes
    double dx,dy,len,force, costh, sinth;
    for(int i=0;i<springlist.size();i++){
        rlen=springlist[i].rlen;
        k=springlist[i].k;
        one=springlist[i].one;
        two=springlist[i].two;
        wlr=springlist[i].wlr;
        wud=springlist[i].wud;
        
        //calculate r_one_two and f_one_two
        xb1=XY(one);
        xb2=XY(two)+wlr;
        yb1=XY(one+num);
        yb2=XY(two+num)+wud;
        
        x1=xb1+strain*yb1;
        y1=yb1;
        x2=xb2+strain*yb2;
        y2=yb2;
        
        dx=x2-x1;
        dy=y2-y1;
        
        len=sqrt(dx*dx+dy*dy);
        costh=dx/len;//cos of angle between the vector describing the spring and the x-axis;
        //sinth=sqrt(1-costh*costh);
        sinth=dy/len;
        force=-k*(len-rlen);
        //cout<<force<<"  "<<springlist[i].sticki<<endl;
        //some outputs to check and debug.
        //cout<<"c        "<<costh<<"     s       "<<sinth<<"     s^2+c^2      "<<costh*costh+sinth*sinth<<"   sticki  "<<springlist[i].sticki<<endl;
        sxx+=force*len*costh*costh;
        sxy+=force*len*sinth*costh;
        syy+=force*len*sinth*sinth;
        
        
    }

    stresstensor S;
    S.sxx=sxx;
    S.sxy=sxy;
    S.syy=syy;
    return S;
}