#include "random.h"
#include <algorithm>
#include "makemikadonetwork.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <vector>
#include "EnergyandGradients.h"
#include<iostream>
using namespace Eigen;
using namespace std;
const double pi=4.0*atan(1.0);
   
int SIGN(double a,double b)
{
    double c=a*b;
    if(c>=0) {
        return 1;
    } else{
        return -1;    
    }
}

int sgn(double x)
{
    if(x>=0) {
        return 1;
    } else{
        return -1;}  
}

double ROSENBROCK(const Eigen::VectorXd &XY)
{
 double f=pow((1-XY(0)),2)+pow((XY(1)-XY(0)*XY(0)),2);
 return f;
}

VectorXd GRAD_rosen(const Eigen::VectorXd &XY)
{
 VectorXd GRAD(2);
 GRAD<<-2*(1-XY(0))-4*XY(0)*(XY(1)-XY(0)*XY(0)),2*(XY(1)-XY(0)*XY(0));
 return GRAD;
}

double dROSENdA(const Eigen::VectorXd &XY,const Eigen::VectorXd &s)
{
 double dda=s.dot(GRAD_rosen(XY));
 return dda;
    
}

double Energynetwork(const vector<spring> &springlist, const VectorXd &XY,
                     const double g11, const double g12, const double g22)
{

  double Energy=0;
  int num=XY.size()/2;
  double k,L;
  double x1,x2,y1,y2;
  int one,two;
    for(int i=0;i<springlist.size();i++){
        k=springlist[i].k;
        L=springlist[i].rlen;
        x1=XY(springlist[i].one);
        x2=XY(springlist[i].two)+springlist[i].wlr;
        y1=XY(springlist[i].one+num);
        y2=XY(springlist[i].two+num)+springlist[i].wud;
        
      Energy=Energy+
      0.5*k*pow(sqrt(
          g11*(x1-x2)*(x1-x2)+
          g22*(y1-y2)*(y1-y2)+
          2*g12*(x1-x2)*(y1-y2))-L,2);
    }

  return Energy;  
}

double distance1(const double x1, const double y1, const double x2,const double y2)
{
 double dist=sqrt(pow((x2-x1),2)+pow((y2-y1),2));
 return dist;
}

double Dist(double x1,double y1,double x2,double y2,double g11,double g12,double g22)
{
 double d;
 d=g11*(x2-x1)*(x2-x1)+2*g12*(x2-x1)*(y2-y1)+g22*(y2-y1)*(y2-y1);
 d=sqrt(d);
 return d;   
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


VectorXd Gradient(const vector<spring> &springlist,
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
  
  double x1,x2,y1,y2;
  double k; 
  double L;
  int one; int two; 
  int num=XY.size()/2;
for(int i=0;i<springlist.size();i++){
    
    one=springlist[i].one;
    two=springlist[i].two;
    x1=X(one);
    x2=X(two)+springlist[i].wlr;
    y1=Y(one);
    y2=Y(two)+springlist[i].wud;
    k=springlist[i].k;
    L=springlist[i].rlen;
    
 
    gradE(one)=gradE(one)+
    k*(sqrt(
        g11*(x1-x2)*(x1-x2)+
        g12*(x1-x2)*(y1-y2)+
        g22*(y1-y2)*(y1-y2))-L)/
        sqrt(
        g11*(x1-x2)*(x1-x2)+
        g12*(x1-x2)*(y1-y2)+
        g22*(y1-y2)*(y1-y2))*
        (g11*(x1-x2)+g12*(y1-y2));
        
    gradE(two)=gradE(two)+
      k*(sqrt(
        g11*(x1-x2)*(x1-x2)+
        g12*(x1-x2)*(y1-y2)+
        g22*(y1-y2)*(y1-y2))-L)/
        sqrt(
        g11*(x1-x2)*(x1-x2)+
        g12*(x1-x2)*(y1-y2)+
        g22*(y1-y2)*(y1-y2))*
        (g11*(x2-x1)+g12*(y2-y2));
        
        
    gradE(one+num)=gradE(one+num)+
        k*(sqrt(
        g11*(x1-x2)*(x1-x2)+
        g12*(x1-x2)*(y1-y2)+
        g22*(y1-y2)*(y1-y2))-L)/
        sqrt(
        g11*(x1-x2)*(x1-x2)+
        g12*(x1-x2)*(y1-y2)+
        g22*(y1-y2)*(y1-y2))*
        (g22*(y1-y2)+g12*(x1-x2));
        
    gradE(two+num)=gradE(two+num)+
        k*(sqrt(
        g11*(x1-x2)*(x1-x2)+
        g12*(x1-x2)*(y1-y2)+
        g22*(y1-y2)*(y1-y2))-L)/
        sqrt(
        g11*(x1-x2)*(x1-x2)+
        g12*(x1-x2)*(y1-y2)+
        g22*(y1-y2)*(y1-y2))*
        (g22*(y2-y1)+g12*(x2-x1));
}
return gradE;  
}

VectorXd gradEbendn(const vector<vector<int>> &springpairs, 
                    const vector<spring> &springlist, 
                    const VectorXd &XY,
                    double g11,double g12, double g22, 
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


double dEda(const VectorXd &XY,const VectorXd &s0,const vector<spring> &springlist,
    const vector<vector<int>> &springpairs,double kappa,
    const double g11,
    const double g12,
    const double g22)
{  
    double out;
    out=s0.dot(Gradient(springlist,XY,g11,g12,g22)+gradEbendn(springpairs,springlist,XY,g11,g12,g22,kappa));
    return out;  
}

  






