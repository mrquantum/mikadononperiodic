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

VectorXd gradL(const double x1,const double y1,const double x2, const double y2, 
               const int springnr, const vector<spring> &springlist,const int num)
{
 VectorXd grad(2*num);
 for(int i=0;i<grad.size();i++){
     grad(i)=0;
 }
 
 int coordNRone=springlist[springnr].one;
 int coordNRtwo=springlist[springnr].two;
 
    grad(coordNRone)=grad(coordNRone)+(1/distance1(x1,y1,x2,y2))*(-(x2-x1));
    grad(coordNRtwo)=grad(coordNRtwo)+(1/distance1(x1,y1,x2,y2))*(x2-x1);
    grad(coordNRone+num)=grad(coordNRone+num)+(1/distance1(x1,y1,x2,y2))*(-(y2-y1));
    grad(coordNRtwo+num)=grad(coordNRtwo+num)+(1/distance1(x1,y1,x2,y2))*(y2-y1);
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

double distance1(const double x1, const double y1, const double x2,const double y2)
{
 double dist=sqrt(pow((x2-x1),2)+pow((y2-y1),2));
 return dist;
}


double Ebend(const vector<vector<int>> &springpairs,
             const vector<spring> &springlist,
             const VectorXd &XY,
            const double kappa)  
{
 double Energy=0;
 int num=XY.size()/2;
 
 for(int i=0;i<springpairs.size();i++){
 int springone=springpairs[i][0];
    int springtwo=springpairs[i][1];
 
    int coordNRone=springlist[springone].one;
    int coordNRtwo=springlist[springone].two;
    int coordNRthree=springlist[springtwo].two;
 
    double x1=XY(coordNRone);
    double y1=XY(coordNRone+num);
 
    double x21=XY(coordNRtwo)+springlist[springone].wlr; //version of (x2,y2) that lies in on spring 1, so possibly outside of the box
    double y21=XY(coordNRtwo+num)+springlist[springone].wud;
 
    double x23=XY(coordNRtwo);                 //version of (x2,y2) that is on spring 2, so MUST be inside the box
    double y23=XY(coordNRtwo+num);
 
    double x3=XY(coordNRthree)+springlist[springtwo].wlr;
    double y3=XY(coordNRthree+num)+springlist[springtwo].wud;

    Vector2d v1,v2;
    v1<<(x21-x1),(y21-y1);
    v2<<(x3-x23),(y3-y23);
    double dotv1v2=v1.dot(v2);
    double lenv1v2=sqrt(v1.dot(v1))*sqrt(v2.dot(v2));
    double c=dotv1v2/lenv1v2;
    
    if(c<-1) c=-1;
    if(c>1) c=1;
    
    double dE=(kappa/(distance1(x1,y1,x21,y21)+distance1(x23,y23,x3,y3)))*pow(acos(c),2);
    Energy=Energy+dE;    
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


VectorXd gradEbend(const vector<vector<int>> &springpairs,const vector<spring> &springlist,const VectorXd &XY,double kappa)
{
  VectorXd grad(XY.size()); //Total gradient
  VectorXd firstpart(grad.size());  //(pi-acos(th))^2 *grad(1/(L1+L2))
  VectorXd secondpart(grad.size()); //1/(L1+L2) grad(pi-acos(th))^2
  VectorXd gradL1L2m1(grad.size()); //grad( 1/(L1L2))
  VectorXd GRADc(grad.size()); //grad cos th in terms of v1.v2/v1v2
  VectorXd GRADnumerator(grad.size()); //grad v1.v2
  VectorXd GRADdenumerator(grad.size());//grad v1v2 , just the lengths
  
  
  
  //Initiate the gradient with zero's
  for(int i=0;i<grad.size();i++){
      grad(i)=0;
      firstpart(i)=0;
      secondpart(i)=0;
  }  

  int num=grad.size()/2;
  double numerator,denumerator;
  double c,s; //c=cos(theta)
 
 
 //add the first term of the gradient arccos(...)^2*grad(1/l12+l23)
  for(int i=0;i<springpairs.size();i++)
  {
    int springone=springpairs[i][0];
    int springtwo=springpairs[i][1];
    int coordNRone=springlist[springone].one;
    int coordNRtwo=springlist[springone].two;
    int coordNRthree=springlist[springtwo].two;
   
   double x1=XY(coordNRone);
   double y1=XY(coordNRone+num);
   double x21=XY(coordNRtwo)+springlist[springone].wlr;
   double y21=XY(coordNRtwo+num)+springlist[springone].wud;

   double x23=XY(coordNRtwo);
   double y23=XY(coordNRtwo+num);

   
   double x3=XY(coordNRthree)+springlist[springtwo].wlr;
   double y3=XY(coordNRthree+num)+springlist[springtwo].wud;
   
for(int j=0;j<gradL1L2m1.size();j++){
       gradL1L2m1(j)=0;
   }
   
   Vector2d v1,v2;
   v1<<(x1-x21),(y1-y21); 
   v2<<(x3-x23),(y3-y23);
   
   //When shearing the box we are going to change this to customized functions:
   numerator=v1.dot(v2);
   denumerator=sqrt(v1.dot(v1))*sqrt(v2.dot(v2));
  
   c=numerator/denumerator;
   if(c>1) c=1;
   if(c<-1) c=-1;
      
   double d12=distance1(x1,y1,x21,y21);
   double d23=distance1(x3,y3,x23,y23);
   
   gradL1L2m1(coordNRone)= (x1-x21)/d12;
   gradL1L2m1(coordNRone+num)=(y1-y21)/d12;
   gradL1L2m1(coordNRtwo)=(x21-x1)/d12+(x23-x3)/d23;
   gradL1L2m1(coordNRtwo+num)=(y21-y1)/d12+(y23-y3)/d23;
   gradL1L2m1(coordNRthree)=(x3-x23)/d23;
   gradL1L2m1(coordNRthree+num)=(y3-y23)/d23;
   
   gradL1L2m1=gradL1L2m1*(-1/(d12+d23)*(d12+d23));
   
   
    firstpart=firstpart+pow((pi-acos(c)),2)*gradL1L2m1;
    
  }
  
  //now the second part = 1/(l1+l2) * grad (pi-theta)^2;
  
  for(int i=0;i<springpairs.size();i++){
   
   
   for(int j=0;j<GRADnumerator.size();j++){ //Set the new gradients to zero
       GRADnumerator(j)=0;
       GRADdenumerator(j)=0;
       GRADc(j)=0;
    }
      
  int springone=springpairs[i][0];
  int springtwo=springpairs[i][1];
            
  int coordNRone=springlist[springone].one;
  int coordNRtwo=springlist[springone].two;
  int coordNRthree=springlist[springtwo].two;
  
  double x1=XY(coordNRone);
  double y1=XY(coordNRone+num);
  
  double x21=XY(coordNRtwo)+springlist[springone].wlr;
  double y21=XY(coordNRtwo+num)+springlist[springone].wud;
  
  double x23=XY(coordNRtwo);
  double y23=XY(coordNRtwo+num);
  
  double x3=XY(coordNRthree)+springlist[springtwo].wlr;
  double y3=XY(coordNRthree+num)+springlist[springtwo].wud;

  //cout<<x1<<"  "<<y1<<"  "<<x21<<"  "<<y21<<"  "<<x23<<"  "<<y23<<"  "<<x3<<"  "<<y3<<endl;
  
  Vector2d v1,v2; 
  v1<<(x1-x21),(y1-y21);
  v2<<(x3-x23),(y3-y23);
  
   double d12=distance1(x1,y1,x21,y21);
   double d23=distance1(x3,y3,x23,y23); 
  
   //numerator=v1.dot(v2);
   numerator=(x1-x21)*(x3-x23)+(y1-y21)*(y3-y23);
   denumerator=d12*d23;
   double c=numerator/denumerator;
   if(c>1) c=1;
   if(c<-1) c=-1;
   //This is the sin(theta);
   s=sqrt(1-c*c);
   if(s<0.0001) s=0.0001;
   s /= s;

    
    
    GRADnumerator(coordNRone)=GRADnumerator(coordNRone)-x23+x3;
    GRADnumerator(coordNRone+num)=GRADnumerator(coordNRone+num)-y23+y3;
    
    GRADnumerator(coordNRtwo)=GRADnumerator(coordNRtwo)-x3+x21+x23-x1;
    GRADnumerator(coordNRtwo+num)=GRADnumerator(coordNRtwo+num)-y3+y21+y23-y1;
    
    GRADnumerator(coordNRthree)=GRADnumerator(coordNRthree)-x21+x1;
    GRADnumerator(coordNRthree+num)=GRADnumerator(coordNRthree+num)-y21+y1;

    
    GRADdenumerator(coordNRone)=GRADdenumerator(coordNRone)+(d23/d12)*(x1-x21);
    GRADdenumerator(coordNRone+num)=GRADdenumerator(coordNRone+num)+(d23/d12)*(y1-y21);
    
    GRADdenumerator(coordNRtwo)=GRADdenumerator(coordNRtwo)+(d23/d12)*(x21-x1);
    GRADdenumerator(coordNRtwo+num)=GRADdenumerator(coordNRtwo+num)+(d23/d12)*(y21-y1);
    
    GRADdenumerator(coordNRtwo)=GRADdenumerator(coordNRtwo)+(d12/d23)*(x23-x3);
    GRADdenumerator(coordNRtwo+num)=GRADdenumerator(coordNRtwo+num)+(d12/d23)*(y23-y3);
    
    GRADdenumerator(coordNRthree)=GRADdenumerator(coordNRthree)+(d12/d23)*(x3-x23);
    GRADdenumerator(coordNRthree+num)=GRADdenumerator(coordNRthree+num)+(d12/d23)*(y3-y23);

    GRADc=(denumerator*GRADnumerator-numerator*GRADdenumerator)/(pow(denumerator,2));
    //cout<<GRADc.dot(GRADc)<<endl;

    double dacos=-1-c*c/2-3*c*c*c*c/8;
    if(GRADc.dot(GRADc)<.1) GRADc=1/s*GRADc;
    secondpart=secondpart+(1.0/(d12+d23))*2*(pi-acos(c))*s*GRADc;
  
  }

  grad=kappa*(firstpart+secondpart);
  return grad; 
}


double dEda(const VectorXd &XY,const vector<anchor> &Anchor,const VectorXd &s0,const vector<spring> &springlist,
    const vector<vector<int>> &springpairs,double kappa)
{  
    double out;
    out=s0.dot((Gradient(springlist,XY,Anchor)+gradEbend(springpairs,springlist,XY,kappa)));
    return out;  
}

