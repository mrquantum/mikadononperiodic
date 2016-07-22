#include <iostream>
#include "structs.h"

#include <algorithm>
#include <vector>
#include <math.h>
#include <functional>
#include <stdio.h>
#include <stdlib.h>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"

using namespace std;
using namespace Eigen;





double StretchEnergy(const vector<spring> &springlist,
                    const VectorXd &XY,
                    double sheardeformation)
{
    int one,two,wlr,wud;
    double x1,x2,y1,y2;
    double X1,X2,Y1,Y2;
    int num=XY.size()/2;
    double gamma=sheardeformation;
    double k,rlen;
    double L;
    double E=0.0;
    
    for(int i=0;i<springlist.size();i++){
        one=springlist[i].one;
        two=springlist[i].two;
        wlr=springlist[i].wlr;
        wud=springlist[i].wud;
        k=springlist[i].k;
        rlen=springlist[i].rlen;
        
        x1=XY(one);
        x2=XY(two)+wlr;
        y1=XY(one+num);
        y2=XY(two+num)+wud;
        
        X1=x1+gamma*y1;
        X2=x2+gamma*y2;
        Y1=y1;
        Y2=y2;

        L=sqrt(pow((X1-X2),2)+pow((Y1-Y2),2));
        E+= k*pow((L-rlen),2);
    }
    return E;
}

double BendEnergy(const vector<vector<int>> &springpairs,
                const vector<spring> &springlist,
                const VectorXd &XY,
                double kappa,
                double sheardeformation)
{
    //This function calculates the Bending energy w. an input of physical coordinates.
    int one,two,three,springone,springtwo,wud1,wud2,wlr1,wlr2;
    double x1,x2,x3,y1,y2,y3;
    double X1,X2,X3,Y1,Y2,Y3;
    double L12,L23;
    int num=XY.size()/2;
    double gamma=sheardeformation;
    double dE;
    double E=0;
    for(int i=0;i<springpairs.size();i++){
        springone=springpairs[i][0];
        springtwo=springpairs[i][1];
        
        one=springlist[springone].one;
        two=springlist[springone].two;
        three=springlist[springtwo].two;
        wlr1=springlist[springone].wlr;
        wud1=springlist[springone].wud;
        wlr2=springlist[springtwo].wlr;
        wud2=springlist[springtwo].wud;
        
        x1=XY(one);
        x2=XY(two)+wlr1;
        x3=XY(three)+wlr1+wlr2;
        
        y1=XY(one+num);
        y2=XY(two+num)+wud1;
        y3=XY(three+num)+wud1+wud2;
        
        X1=x1+y1*sheardeformation;
        X2=x2+y2*sheardeformation;
        X3=x3+y3*sheardeformation;
        Y1=y1;
        Y2=y2;
        Y3=y3;
        
        
        
        L12=sqrt(pow((X1-X2),2)+pow((Y1-Y2),2));
        L23=sqrt(pow((x2-x3),2)+pow((Y2-Y3),2));

        dE=(kappa/(L12+L23))*(pow((-X1 + X2)/L12 - (-X2 + X3)/L23,2)+pow((-Y1 + Y2)/L12 - (-Y2 + Y3)/L23,2));
        E+=dE;        
    }
    return E;
}



VectorXd EffKappa(const vector<vector<int>> &springpairs,
                  const vector<spring> &springlist,
                  const VectorXd &XY,
                  double kappa)

//This function calculates the effective bendingrigidity := kappa/(L12+L23) for all 
//spring - pairs. This only works in physical coordinates! So better do it before any
//shear deformation.

{
    int size=springpairs.size();
    int one,two,three,wud1,wud2,wlr1,wlr2,springone,springtwo;
    double x1,x2,x3,y1,y2,y3,L12,L23;
    
    VectorXd effkappa(springpairs.size());
    
    
    int num=XY.size()/2;
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
        
        x1=XY(one);
        x2=XY(two)+wlr1;
        x3=XY(three)+wlr1+wlr2;
        y1=XY(one+num);
        y2=XY(two+num)+wud1;
        y3=XY(three+num)+wud1+wud2;
        
        L12=sqrt(pow((x1-x2),2)+pow((y1-y2),2));
        L23=sqrt(pow((x2-x3),2)+pow((y2-y3),2));
        
        effkappa(i)=kappa/(L12+L23);
        
    }
    return effkappa;
}





VectorXd SimpleBendingGrad(const triplet &Springpairs,
                           const vector<spring> &springlist,
                           const VectorXd &XY,
                           double sheardeformation) //XY in
{
    int one,two,three;
    int wlr1,wlr2,wud1,wud2;
    double x1,x2,x3,y1,y2,y3; //Boxcoordinates
    double X1,X2,X3,Y1,Y2,Y3; //REAL coordinates
    //double X1phys,X2phys,X3phys,Y1phys,Y2phys,Y3phys;
    double gradx1,gradx2,gradx3,grady1,grady2,grady3;
    double L12,L23,kappa;
    int num=XY.size()/2;
    int springone, springtwo;
    double gamma=sheardeformation;
    vector<vector<int>> triplets=Springpairs.springpairs;
    VectorXd effkappa=Springpairs.EffKappa;
    
    VectorXd BGrad(XY.size());
    for(int i=0;i<XY.size();i++){
        BGrad(i)=0;
    }
    

    for(int i=0;i<triplets.size();i++){
        kappa=effkappa[i];
        springone=triplets[i][0];
        springtwo=triplets[i][1];
        
        one=springlist[springone].one;
        two=springlist[springone].two;
        three=springlist[springtwo].two;
        
        wlr1=springlist[springone].wlr;
        wlr2=springlist[springtwo].wlr;
        wud1=springlist[springone].wud;
        wud2=springlist[springtwo].wud;
        
        //These are the box-coordinates
        x1=XY(one);
        x2=XY(two)+wlr1;
        x3=XY(three)+wlr1+wlr2;
        y1=XY(one+num);
        y2=XY(two+num)+wud1;
        y3=XY(three+num)+wud1+wud2;
        
        //Here the physical coordinates
        X1=x1+gamma*y1;
        X2=x2+gamma*y2;
        X3=x3+gamma*y3;
        Y1=y1;
        Y2=y2;
        Y3=y3;

        L12=sqrt(pow((X1-X2),2)+ pow((Y1-Y2),2));
        L23=sqrt(pow((X2-X3),2)+ pow((Y2-Y3),2));
        
        gradx1=((2*(Y1 - Y2)*(X3*(Y1 - Y2) + X1*(Y2 - Y3) + X2*(-Y1 + Y3)))/(pow(L12,3)*L23));
        gradx2=((2*(X3*(-Y1 + Y2) + X2*(Y1 - Y3) + X1*(-Y2 + Y3))*(pow(X3,2)*(Y1 - Y2) + 
                            pow(X2,2)*(Y1 - Y3) + (pow(X1,2) + (Y1 - Y2)*(Y1 - Y3))*(Y2 - Y3) + 2*X2*(X3*(-Y1 + Y2) + X1*(-Y2 + Y3))))/
                            (pow(L12,3)*pow(L23,3)));
        gradx3=((2*(Y2 - Y3)*(X3*(Y1 - Y2) + X1*(Y2 - Y3) + X2*(-Y1 + Y3)))/(L12*pow(L23,3)));
        grady1=((-2*(X1 - X2)*(X3*(Y1 - Y2) + X1*(Y2 - Y3) + X2*(-Y1 + Y3)))/(pow(L12,3)*L23));
        grady2=((2*(X3*(Y1 - Y2) + X1*(Y2 - Y3) + X2*(-Y1 + Y3))*(pow(X1,2)*(X2 - X3) +
                        pow(X2,2)*X3 - X3*pow(Y1 - Y2,2) + X1*(-pow(X2,2) + pow(X3,2) + pow(Y2 - Y3,2)) - X2*(pow(X3,2) -
                        (Y1 - Y3)*(Y1 - 2*Y2 + Y3))))/(pow(L12,3)*pow(L23,3)));
        grady3=((2*(X2 - X3)*(X3*(-Y1 + Y2) + X2*(Y1 - Y3) + X1*(-Y2 + Y3)))/(L12*pow(L23,3)));
        
        
        //Transform the gradient in box-coordinates
        BGrad(one)+=effkappa(i)*gradx1;
        BGrad(two)+=effkappa(i)*gradx2;
        BGrad(three)+=effkappa(i)*gradx3;
        BGrad(one+num)+=effkappa(i)*(grady1+gamma*gradx1);
        BGrad(two+num)+=effkappa(i)*(grady2+gamma*gradx2);
        BGrad(three+num)+=effkappa(i)*(grady3+gamma*gradx3);
    }
    return BGrad;
}



void HarmonicGradPhys(const vector<spring> &springlist,
                         double *XY,
                         double *xi,
                         int size,
                         double sheardeformation)
//xi is the Harmonic gradient in physical coordinates, XY is in boxcoordinates.
{
    int one,two,wlr,wud;
    double x1,x2,y1,y2;
    double X1,X2,Y1,Y2; //Physcoordinates
    double gradx1,gradx2,grady1,grady2;
    double L,rlen,k;
    double gradx,grady;
    double gamma=sheardeformation;
    int num=size/2;
    
    
    for(int i =0;i<springlist.size();i++){
        one=springlist[i].one;
        two=springlist[i].two;
        wlr=springlist[i].wlr;
        wud=springlist[i].wud;
        k=springlist[i].k;
        rlen=springlist[i].rlen;
        
        x1=XY[one];
        y1=XY[one+num];
        x2=XY[two]+wlr;
        y2=XY[two+num]+wud;
        
        X1=x1+gamma*y1;
        X2=x2+gamma*y2;
        Y1=y1;
        Y2=y2;
        
        L=sqrt(pow((X1-X2),2)+pow((Y1-Y2),2));
        
        gradx=(2*k*(L - rlen)*(X1 - X2))/L;
        grady=(2*k*(L - rlen)*(Y1 - Y2))/L;
        
        gradx1=gradx;
        gradx2=-gradx;
        grady1=grady+gamma*gradx1;
        grady2=-grady+gamma*gradx2;
        
        xi[one]+=gradx1;
        xi[two]+=gradx2;
        xi[one+num]+=grady1;
        xi[two+num]+=grady2;
    }

}








