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
   
double distance(const Eigen::VectorXd &XY,int one, int two)
{
  double l;
  int num=XY.size()/2;
  l=sqrt(pow(XY(two)-XY(one),2)+pow(XY(two+num)-XY(one+num),2));
  return l;
}
  
double dldxi(const Eigen::VectorXd &XY,int one, int two)
{
   double deri;
   deri=-1*(1/distance(XY,one,two))*(XY(two)-XY(one));
return deri;
}  
  
double dldyi(const Eigen::VectorXd &XY,int one,int two)
{
  double deri;
  int num=XY.size()/2;
  deri=-1*(1/distance(XY,one,two))*(XY(two+num)-XY(one+num));
return deri;
}

VectorXd gradL(const Eigen::VectorXd &XY,int one,int two)
{
 VectorXd grad(XY.size());
 for(int i=0;i<XY.size();i++){
     grad(i)=0;
 }
 
 int num=XY.size()/2;
 double x1=XY(one);
 double x2=XY(two);
 double y1=XY(one+num);
 double y2=XY(two+num);
 
    grad(one)=grad(one)+(1/distance(XY,one,two))*(-(x2-x1));
    grad(two)=grad(two)+(1/distance(XY,one,two))*(x2-x1);
    grad(one+num)=grad(one+num)+(1/distance(XY,one,two))*(-(y2-y1));
    grad(two+num)=grad(two+num)+(1/distance(XY,one,two))*(y2-y1);
return grad;
}  
    
double Energynetwork(const vector<spring> &springlist, const VectorXd &XY,const vector<anchor> &Anchor)
{
  VectorXd X(XY.size()/2);
  VectorXd Y(XY.size()/2);
  X=XY.head(XY.size()/2);
  Y=XY.tail(XY.size()/2);
  double Energy=0;
    for(int i=0;i<springlist.size();i++){
      Energy=Energy+
      0.5*springlist[i].k*pow(sqrt(pow(X(springlist[i].one)-(X(springlist[i].two)+springlist[i].wlr),2)
      +pow(Y(springlist[i].one)-(Y(springlist[i].two)+springlist[i].wud),2))-springlist[i].rlen,2);
    }
//  for(int j=0;j<Anchor.size();j++){
//   Energy=Energy+0.5*k*pow(sqrt(pow(X(Anchor[j].label)-Anchor[j].xpos,2)+pow(Y(Anchor[j].label)-Anchor[j].ypos,2))-L,2);
//  }
  return Energy;  
}

double distance1(const VectorXd &XY,const double x1, const double y1, const double x2,const double y2)
{
 double dist=sqrt(pow((x2-x1),2)+pow((y2-y1),2));
 return dist;
}




double Ebend(const vector<vector<int>> &springpairs,
             const vector<spring> &springlist,
             const VectorXd &XY)  
{
 double Energy=0;
 double kappa=1;
 int num=XY.size()/2;
 
 for(int i=0;i<springpairs.size();i++){
 int springone=springpairs[i][0];
 int springtwo=springpairs[i][1];
 int coordNRone=springlist[springone].one;
 int coordNRtwo=springlist[springone].two;
 int coordNRthree=springlist[springtwo].two;
 
 double x1=XY(coordNRone);
 double y1=XY(coordNRone+num);
 double x2=XY(coordNRtwo);
 double y2=XY(coordNRtwo+num);
 double x3=XY(coordNRthree);
 double y3=XY(coordNRthree+num);
 
 double c=x2*x3-x2*x2-x1*x3+x1*x2+y2*y3-y2*y2-y1*y3+y1*y2; //cos theta
 if(c<-1) c=-1;
 if(c>1) c=1;
 
 Energy=Energy+
    (0.5*kappa/(distance1(XY,x1,y1,x2,y2)+distance1(XY,x2,y2,x3,y3)))*pow(acos(c),2);
 }
 
 return Energy; 
}
  


VectorXd Gradient(const vector<spring> &springlist,const VectorXd &XY,const vector<anchor> &Anchor)
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
 gradE(one)=gradE(one)+
 springlist[i].k*(sqrt(pow(X(one)-(X(two)+springlist[i].wlr),2)
 +pow(Y(one)-(Y(two)+springlist[i].wud),2))-springlist[i].rlen)*
 (X(one)-(X(two)+springlist[i].wlr))/
 sqrt(pow(X(one)-(X(two)+springlist[i].wlr),2)
 +pow(Y(one)-(Y(two)+springlist[i].wud),2));
 
 gradE(two)=gradE(two)-
 springlist[i].k*(sqrt(pow(X(one)-(X(two)+springlist[i].wlr),2)
 +pow(Y(one)-(Y(two)+springlist[i].wud),2))-springlist[i].rlen)*
 (X(one)-(X(two)+springlist[i].wlr))/
 sqrt(pow(X(one)-(X(two)+springlist[i].wlr),2)
 +pow(Y(one)-(Y(two)+springlist[i].wud),2));

 gradE(one+num)=gradE(one+num)+
  springlist[i].k*(sqrt(pow(X(one)-(X(two)+springlist[i].wlr),2)
 +pow(Y(one)-(Y(two)+springlist[i].wud),2))-springlist[i].rlen)*
 (Y(one)-(Y(two)+springlist[i].wud))/
 sqrt(pow(X(one)-(X(two)+springlist[i].wlr),2)
 +pow(Y(one)-(Y(two)+springlist[i].wud),2));
 
  gradE(two+num)=gradE(two+num)-
  springlist[i].k*(sqrt(pow(X(one)-(X(two)+springlist[i].wlr),2)
 +pow(Y(one)-(Y(two)+springlist[i].wud),2))-springlist[i].rlen)*
 (Y(one)-(Y(two)+springlist[i].wud))/
 sqrt(pow(X(one)-(X(two)+springlist[i].wlr),2)
 +pow(Y(one)-(Y(two)+springlist[i].wud),2));
}
return gradE;  
}  


VectorXd gradEbend(const vector<triplet> &tripl,const Eigen::VectorXd &XY)
{
  VectorXd grad(XY.size());
  VectorXd firstpart(grad.size());
  VectorXd secondpart(grad.size());

  //Initiate the gradient with zero
  for(int i=0;i<grad.size();i++){
      grad(i)=0;
      firstpart(i)=0;
      secondpart(i)=0;
  }  

  int num=grad.size()/2;
  double kappa=1;
  double numerator,denumerator;
  double c; //c=cos(theta)
  //add the first term of the gradient arccos(...)^2*grad(1/l12+l23)
  for(int i=0;i<tripl.size();i++)
  {
    int one=tripl[i].one;
    int two=tripl[i].two;
    int three=tripl[i].three;
    numerator=XY(two)*XY(three)-pow(XY(two),2)-XY(one)*XY(three)+XY(one)*XY(two)+
		XY(two+num)*XY(three+num)-pow(XY(two+num),2)-XY(one+num)*XY(three+num)+
		XY(one+num)*XY(two+num);
    denumerator=distance(XY,one,two)*distance(XY,two,three); 

   c=numerator/denumerator;
   if(c>1) c=1;
   if(c<-1) c=-1;
   firstpart=firstpart+pow(acos(c),2)*(-1/pow((distance(XY,one,two)+distance(XY,two,three)),2))*
        (gradL(XY,one,two)+gradL(XY,two,three));
   } 
  
  //now the second part = 1/(l1+l2) * grad theta^2;
  for(int i=0;i<tripl.size();i++){
   int one=tripl[i].one;
   int two=tripl[i].two;
   int three=tripl[i].three;
   numerator=XY(two)*XY(three)-pow(XY(two),2)-XY(one)*XY(three)+XY(one)*XY(two)+
                XY(two+num)*XY(three+num)-pow(XY(two+num),2)-XY(one+num)*XY(three+num)+
                XY(one+num)*XY(two+num);
   denumerator=distance(XY,one,two)*distance(XY,two,three); 
   
   //This is the argument of the arccos. Arccos(c)=arccos(cos(theta))=theta --> c=cos(theta)
   double c=numerator/denumerator;
   if(c>1) c=1;
   if(c<-1) c=-1;
   //This is the sin(theta);
   double s=sqrt(1-c*c);
   if(s<0.0001) s=0.0001;
   
   VectorXd GRADc(secondpart.size());
   VectorXd GRADnumerator(secondpart.size());
   VectorXd GRADdenumerator(secondpart.size());
   
   for(int j=0;j<GRADnumerator.size();j++){
       GRADnumerator(j)=0;
       GRADdenumerator(j)=0;
    }
   GRADnumerator(one)=(XY(two)-XY(three));
   GRADnumerator(two)=(XY(three)-2*XY(two)+XY(one));
   GRADnumerator(three)=(XY(two)-XY(one));
   GRADnumerator(one+num)=(XY(two+num)-XY(three+num));
   GRADnumerator(two+num)=(XY(three+num)-2*XY(two+num)+XY(one+num));
   GRADnumerator(three+num)=(XY(two+num)-XY(one+num));

   GRADdenumerator=distance(XY,two,three)*gradL(XY,one,two)+distance(XY,one,two)*gradL(XY,one,two);
   GRADc=(denumerator*GRADnumerator-numerator*GRADdenumerator)/(denumerator*denumerator);
  
   secondpart=secondpart+
   (1/(distance(XY,one,two)+distance(XY,two,three)))*
   (2*acos(c)*(-1/s))*GRADc;
  
  }
  grad=.5*kappa*(firstpart+secondpart);
   return grad; 
}



double dEda(const VectorXd &XY,const vector<anchor> &Anchor,const VectorXd &s0,const vector<spring> &springlist)
{  
    double out;
    out=s0.dot(Gradient(springlist,XY,Anchor));
    return out;  
}

