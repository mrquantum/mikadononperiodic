#include "random.h"
#include <algorithm>
#include "makemikadonetwork.h"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include <vector>
#include "EnergyandGradients.h"
#include "BendingGrad.h"
#include<iostream>
using namespace Eigen;
using namespace std;
const double pi=4.0*atan(1.0);


double Dist(double x1,double y1,double x2,double y2,double g11,double g12,double g22)
{
 double d;
 d=g11*(x2-x1)*(x2-x1)+2*g12*(x2-x1)*(y2-y1)+g22*(y2-y1)*(y2-y1);
 d=sqrt(d);
 return d;   
}

double Energynetwork(const vector<spring> &springlist, const VectorXd &XY,
                     const double g11, const double g12, const double g22)
{
  double Energy=0;
  int num=XY.size()/2;
  double k,L;
  double x1,x2,y1,y2;
  int one,two;
  double dE;

    for(int i=0;i<springlist.size();i++){
        //cout<<"nrs"<<springlist[i].one<<" "<<springlist[i].one+num<<"     "<<springlist[i].two<<" "<<springlist[i].two+num<<endl;
        k=springlist[i].k;
        L=springlist[i].rlen;
        x1=XY(springlist[i].one);
        x2=XY(springlist[i].two)+springlist[i].wlr;
        y1=XY(springlist[i].one+num);
        y2=XY(springlist[i].two+num)+springlist[i].wud;
        dE= 0.5*k*pow(sqrt(
            g11*(x1-x2)*(x1-x2)+
            g22*(y1-y2)*(y1-y2)+
            2*g12*(x1-x2)*(y1-y2))-L,2);
        Energy=Energy+dE;
    }
  return Energy;  
}

double Ebend(const vector<vector<int>> &springpairs,
             const vector<spring> &springlist,
             const VectorXd &XY,
             const double g11,
             const double g12,
             const double g22,
             const double kappa)
{
 double Energy=0;
 int num=XY.size()/2;
 double l1,l2,costh;
 double x1,y1,x21,y21,x23,y23,x3,y3;
 for(int i=0;i<springpairs.size();i++){
 int springone=springpairs[i][0];
    int springtwo=springpairs[i][1];
 
    int coordNRone=springlist[springone].one;
    int coordNRtwo=springlist[springone].two;
    int coordNRthree=springlist[springtwo].two;
 
    x1=XY(coordNRone);
    y1=XY(coordNRone+num);
 
    x21=XY(coordNRtwo)+springlist[springone].wlr; //version of (x2,y2) that lies in on spring 1, so possibly outside of the box
    y21=XY(coordNRtwo+num)+springlist[springone].wud;
 
    x23=XY(coordNRtwo);                 //version of (x2,y2) that is on spring 2, so MUST be inside the box
    y23=XY(coordNRtwo+num);
 
    x3=XY(coordNRthree)+springlist[springtwo].wlr;
    y3=XY(coordNRthree+num)+springlist[springtwo].wud;

    l1=Dist(x1,y1,x21,y21,g11,g12,g22);
    l2=Dist(x23,y23,x3,y3,g11,g12,g22);
    costh=(g11*(x21-x1)*(x3-x23)+g12*(x21-x1)*(y3-y23)+g12*(y21-y1)*(x3-x23)+g22*(y21-y1)*(y3-y23))/(l1*l2);
    if(costh>1.0) costh=1.0;
    if(costh<-1.0) costh=-1.0;

    Energy=Energy+kappa*pow(pi-acos(costh),2)/(l1+l2);
}

return Energy; 
}


// work on a version of the energy AND gradient that the new 
// version of the CG-algorithm accepts 
double EnergyNetworkn(double *XY,networkinfo parameters){
    const vector<spring> springlist=parameters.springlist;
    double g11=parameters.g11;
    double g12=parameters.g12;
    double g22=parameters.g22;
    int size=parameters.size;
    
    double Energy=0;
    int num=size/2;
    double k,L;
    double x1,x2,y1,y2;
    int one,two;
    double dE;

    for(int i=0;i<springlist.size();i++){
        //cout<<"nrs"<<springlist[i].one<<" "<<springlist[i].one+num<<"     "<<springlist[i].two<<" "<<springlist[i].two+num<<endl;
        k=springlist[i].k;
        L=springlist[i].rlen;
        x1=XY[springlist[i].one];
        x2=XY[springlist[i].two]+springlist[i].wlr;
        y1=XY[springlist[i].one+num];
        y2=XY[springlist[i].two+num]+springlist[i].wud;
        dE= 0.5*k*pow(sqrt(
            g11*(x1-x2)*(x1-x2)+
            g22*(y1-y2)*(y1-y2)+
            2*g12*(x1-x2)*(y1-y2))-L,2);
        Energy=Energy+dE;
    }
  return Energy; 

}


void HarmonicGradientn(double *p,double *xi,networkinfo params)
//The harmonic-gradient as the new cgmethod wants it
{
    double g11=params.g11;
    double g12=params.g12;
    double g22=params.g22;
    vector<spring> springlist=params.springlist;
    int size=params.size;
    
    int one,two,num=size/2;
    double dx,dy,k,L,dist;
    double gradx,grady;
        
    double *x=p;
    double *y=p+num;
    
    //make sure that the gradient is zero to begin with;
    for(int i=0;i<size;i++){
        xi[i]=0;
    }
    
  
    for(int i=0;i<springlist.size();i++){
        one=springlist[i].one;
        two=springlist[i].two;
        dx=x[one]-(x[two]+springlist[i].wlr);
        dy=y[one]-(y[two]+springlist[i].wud);
        k=springlist[i].k;
        L=springlist[i].rlen;
        dist=sqrt( g11*dx*dx+ 2*g12*dx*dy+ g22*dy*dy );
   
        gradx= k*(dist-L)*(g11*dx+g12*dy)/dist;
        grady= k*(dist-L)*(g22*dy+g12*dx)/dist;

        xi[one] += gradx;
        xi[two] -= gradx;
        xi[one+num] += grady;
        xi[two+num] -= grady;
    }
}









VectorXd HarmonicGradient(const vector<spring> &springlist,
                  const VectorXd &XY,
                  const double g11,
                  const double g12,
                  const double g22)
{
  VectorXd gradE(XY.size());
  for(int i=0;i<gradE.size();i++){
    gradE(i)=0;
  }
  VectorXd X(XY.size()/2);
  VectorXd Y(XY.size()/2);
  X=XY.head(XY.size()/2);
  Y=XY.tail(XY.size()/2);

  for(int i=0;i<springlist.size();i++){
    int one=springlist[i].one;
    int two=springlist[i].two;
    int num=XY.size()/2;
    double dx=X(one)-(X(two)+springlist[i].wlr);
    double dy=Y(one)-(Y(two)+springlist[i].wud);
    double k=springlist[i].k;
    double L=springlist[i].rlen;
    double dist=sqrt( g11*dx*dx+ 2*g12*dx*dy+ g22*dy*dy );
   
    double gradx= k*(dist-L)*(g11*dx+g12*dy)/dist;
    double grady= k*(dist-L)*(g22*dy+g12*dx)/dist;

    gradE(one) += gradx;
    gradE(two) -= gradx;
    gradE(one+num) += grady;
    gradE(two+num) -= grady;
  }
return gradE;
}

VectorXd gradEbend(const vector<vector<int>> &springpairs, 
                    const vector<spring> &springlist, 
                    const VectorXd &XY,
                    double g11,
                    double g12, 
                    double g22, 
                    double kappa)
{
    VectorXd grad(XY.size()); //Total gradient
    VectorXd gradL1(XY.size());
    VectorXd gradL2(XY.size());
    VectorXd gradC(XY.size()), gradNum(XY.size()), gradDenom(XY.size());
    VectorXd firstpart(XY.size()),secondpart(XY.size());
    
    
    int num=XY.size()/2;
    int springone,springtwo;
    int coordNRone, coordNRtwo, coordNRthree;
    double x1,y1,x21,y21,x23,y23,x3,y3;
    double l1,l2;
    double costh, sinth;
    double numerator, denominator;
    
    for(int i=0;i<grad.size();i++){
        grad(i)=0;
        firstpart(i)=0;
        secondpart(i)=0;
        
    }
    
    //Loop over all springpairs (=angles) to find the bending energy.
    for(int i=0;i<springpairs.size();i++){
        
        //Make sure gradL are zeros
        for(int j=0;j<gradL1.size();j++){
            gradL1(j)=0;
            gradL2(j)=0;
            gradNum(j)=0;
            gradDenom(j)=0;
        }
        
        springone=springpairs[i][0];
        springtwo=springpairs[i][1];
        coordNRone=springlist[springone].one;
        coordNRtwo=springlist[springone].two;
        coordNRthree=springlist[springtwo].two;
        
        x1=XY(coordNRone);
        y1=XY(coordNRone+num);
        x21=XY(coordNRtwo)+springlist[springone].wlr;
        y21=XY(coordNRtwo+num)+springlist[springone].wud;
        x23=XY(coordNRtwo);
        y23=XY(coordNRtwo+num);
        x3=XY(coordNRthree)+springlist[springtwo].wlr;
        y3=XY(coordNRthree)+springlist[springtwo].wud;
        
        l1=Dist(x1,y1,x21,y21,g11,g12,g22);
        l2=Dist(x23,y23,x3,y3,g11,g12,g22);
        
        gradL1(coordNRone)=g11*(x1-x21)+g12*(y1-y21);
        gradL1(coordNRone+num)=g12*(x1-x21)+g22*(y1-y21);
        gradL1(coordNRtwo)=g11*(x21-x1)+g12*(y21-y1);
        gradL1(coordNRtwo+num)=g12*(x21-x1)+g22*(y21-y1);
        gradL1=gradL1/l1;
        
        gradL2(coordNRtwo)=g11*(x23-x3)+g12*(y23-y3);
        gradL2(coordNRtwo+num)=g12*(x23-x3)+g22*(y23-y3);
        gradL2(coordNRthree)=g11*(x3-x23)+g12*(y3-y23);
        gradL2(coordNRthree+num)=g12*(x3-x23)+g22*(y3-y23);
        gradL2=gradL2/l2;
        
        //the cos of the angle with the metric tensor
        costh=(g11*(x21-x1)*(x3-x23)+g12*(x21-x1)*(y3-y23)+g12*(y21-y1)*(x3-x23)+g22*(y21-y1)*(y3-y23))/(l1*l2);
        
        if(costh>1) costh=1;
        if(costh<-1) costh=-1;
        
        sinth=sqrt(1-costh*costh);
        if(sinth<0.0001) sinth=0.0001;
        sinth=1/sinth;

        firstpart=firstpart-pow(pi-acos(costh),2)*((gradL1+gradL2)/((l1+l2)*(l1+l2)));
        
        numerator=g11*(x21-x1)*(x3-x21)+
            g12*(x21-x1)*(y3-y23)+
            g12*(y21-y1)*(x3-x23)+
            g22*(y21-y1)*(y3-y23);
        denominator=l1*l2;

        //Calculate the gradient of costh=num/denom;
        gradNum(coordNRone)=g11*(x23-x3)+g12*(y23-y3);
        gradNum(coordNRone+num)=g12*(x23-x3)+g22*(y23-y3);
        gradNum(coordNRtwo)=g11*(x3-x23)+g11*(x1-x21)+g12*(y3-y23)+g12*(y1-y21);
        gradNum(coordNRtwo+num)=g12*(x1-x21)+g12*(x3-x23)+g22*(y3-y23)+g22*(y1-y21);
        gradNum(coordNRthree)=g11*(x21-x1)+g12*(y21-y1);
        gradNum(coordNRthree+num)=g12*(x21-x1)+g22*(y21-y1);
        gradDenom=l1*gradL2+l2*gradL1;
        gradC=(denominator*gradNum-numerator*gradDenom)/(denominator*denominator);
        
        secondpart=secondpart+2*(pi-acos(costh))*sinth*gradC/(l1+l2);
    }
    grad=kappa*(secondpart+firstpart);
    return grad;
    
}

double dEda(const VectorXd &XY,
            const VectorXd &s0,
            const vector<spring> &springlist,
            const vector<vector<int>> &springpairs,
            double kappa,
            const double g11,
            const double g12,
            const double g22)
{  
    double out;
    out=s0.dot(HarmonicGradient(springlist,XY,g11,g12,g22)+BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22));
    return out;
}

  






