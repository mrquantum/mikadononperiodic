#include <iostream>
#include <vector>
#include <math.h>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"

using namespace std;
using namespace Eigen;

void bendinggradnew(double *XY,double* grad,int size,vector<spring> &springlist,vector<vector<int> > &springpairs,double b_rigid,double gamma)
//This is the BendingGradient in Physical coordinates, but the input XY is in BOX coordinates
{
    
    for(int i=0;i<size;i++){
        grad[i]=0.0;
    }
    
    double kappa=b_rigid;
    
    int num=size/2;
    int one, two, three;
    int wlr1,wlr2,wud1,wud2;
    int springone, springtwo;
    double xb1,xb2,xb3,yb1,yb2,yb3; //the box coordinates 
    double xp1,xp2,xp3,yp1,yp2,yp3; //the phys coordinates
    double L12, L23;
    double pref;
    for(int i=0;i<springpairs.size();i++){
        springone=springpairs[i][0];
        springtwo=springpairs[i][1];
        one=springlist[springone].one;
        two=springlist[springone].two;
        three=springlist[springtwo].two;
        wlr1=springlist[springone].wlr;
        wlr2=springlist[springtwo].wlr;
        wud1=springlist[springone].wud;
        wud2=springlist[springtwo].wud;
        
        xb1=XY[one];
        xb2=XY[two]+wlr1;
        xb3=XY[three]+wlr1+wlr2;
        yb1=XY[one+num];
        yb2=XY[two+num]+wud1;
        yb3=XY[three+num]+wud1+wud2;
        
        xp1=xb1+gamma*yb1;
        xp2=xb2+gamma*yb2;
        xp3=xb3+gamma*yb3;
        yp1=yb1;
        yp2=yb2;
        yp3=yb3;
        
        // and then here the componets of the gradient.
        L12=sqrt(pow(xp2-xp1,2)+pow(yp2-yp1,2));
        L23=sqrt(pow(xp3-xp2,2)+pow(yp3-yp2,2));
        pref=pow(L12,2)*pow(L23,2)*pow((L12+L23),2);
        
        
        //Tedious but can't get it nicer... 
        
        grad[one]+=
        kappa*(2*L23*(L12+L23)*((pow(L12,2)-pow(xp1-xp2,2))*(L23*(xp1-xp2)+L12*(-xp2 + xp3)) + 
        (xp1-xp2)*(yp1-yp2)*(L23*(-yp1+yp2)+L12*(yp2-yp3)))+L12*(-xp1+xp2)*(pow(L23*(xp1-xp2)+L12*(-xp2+xp3),2)+pow(L23*(yp1-yp2)
        +L12*(-yp2+yp3),2)))/(pow(L12,2)*pref);
        
        grad[two]+=
        kappa*(2*(L12+L23)*((-(pow(L12,2)*pow(L23,3))+pow(L23,3)*pow(xp1 - xp2,2) +pow(L12,3)*(-pow(L23,2)+pow(xp2 - xp3,2)))*
        (L23*(xp1-xp2)+L12*(-xp2 + xp3))+(pow(L23,3)*(xp1-xp2)*(yp1-yp2)+pow(L12,3)*(xp2-xp3)*(yp2 - yp3))*
        (L23*(yp1-yp2)+L12*(-yp2 + yp3)))-L12*L23*(L23*(-xp1 + xp2) + L12*(xp2 - xp3))*
        (pow(L23*(xp1-xp2)+L12*(-xp2 + xp3),2)+pow(L23*(yp1-yp2)+L12*(-yp2 + yp3),2)))/(pow(L12,2)*pow(L23,2)*pref);
        
        grad[three]+=
        kappa*(2*L12*(L12 + L23)*((L23*(-xp1 + xp2) + L12*(xp2 - xp3))*(-pow(L23,2) + pow(xp2 - xp3,2)) + 
        (xp2 - xp3)*(L23*(-yp1 + yp2) + L12*(yp2 - yp3))*(yp2 - yp3)) - 
        L23*(-xp2 + xp3)*(pow(L23*(xp1 - xp2) + L12*(-xp2 + xp3),2) + pow(L23*(yp1 - yp2) + L12*(-yp2 + yp3),2)))/(pow(L23,2)*pref);
        
        grad[one+num]+=
        kappa*(2*L23*(L12 + L23)*((xp1 - xp2)*(L23*(-xp1 + xp2) + L12*(xp2 - xp3))*(yp1 - yp2) + 
        (pow(L12,2) - pow(yp1 - yp2,2))*(L23*(yp1 - yp2) + L12*(-yp2 + yp3))) + 
        L12*(-yp1 + yp2)*(pow(L23*(xp1 - xp2) + L12*(-xp2 + xp3),2) + pow(L23*(yp1 - yp2) + L12*(-yp2 + yp3),2)))/(pow(L12,2)*pref);
        
        grad[two+num]+=
        kappa*(2*(L12 + L23)*((L23*(xp1 - xp2) + L12*(-xp2 + xp3))*(pow(L23,3)*(xp1 - xp2)*(yp1 - yp2) + pow(L12,3)*(xp2 - xp3)*(yp2 - yp3)) + 
        (-(pow(L12,2)*pow(L23,3)) + pow(L23,3)*pow(yp1 - yp2,2) + pow(L12,3)*(-pow(L23,2) + pow(yp2 - yp3,2)))*
        (L23*(yp1 - yp2) + L12*(-yp2 + yp3))) - L12*L23*(L23*(-yp1 + yp2) + L12*(yp2 - yp3))*
        (pow(L23*(xp1 - xp2) + L12*(-xp2 + xp3),2) + pow(L23*(yp1 - yp2) + L12*(-yp2 + yp3),2)))/(pow(L12,2)*pow(L23,2)*pref);
        
        grad[three+num]+=
        kappa*(2*L12*(L12 + L23)*((L23*(-yp1 + yp2) + L12*(yp2 - yp3))*(-pow(L23,2) + pow(yp2 - yp3,2)) + 
        (L23*(-xp1 + xp2) + L12*(xp2 - xp3))*(xp2 - xp3)*(yp2 - yp3)) - 
        L23*(-yp2 + yp3)*(pow(L23*(xp1 - xp2) + L12*(-xp2 + xp3),2) + pow(L23*(yp1 - yp2) + L12*(-yp2 + yp3),2)))/(pow(L23,2)*pref);
        
    }
    
    //transforms back to box coordinates
    for(int i=0;i<num;i++){
        grad[i+num]+=grad[i]*gamma;
    }

}