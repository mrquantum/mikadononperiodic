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
  
for(int i=0;i<springlist.size();i++){
    
    int one=springlist[i].one;
    int two=springlist[i].two;
    int num=XY.size()/2;
    double x1,x2,y1,y2;
    x1=X(one);
    x2=X(two)+springlist[i].wlr;
    y1=Y(one);
    y2=Y(two)+springlist[i].wud;
    double k=springlist[i].k;
    double L=springlist[i].rlen;
    
//     gradE(one)=gradE(one)+
//         springlist[i].k*
//         (sqrt(
//             pow((x1-x2),2)+pow(y1-y2,2)-springlist[i].rlen)*
//             (x1-x2)/sqrt(pow(x1-x2,2)+pow(y1-y2,2)) );
 
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

double Dist(double x1,double y1,double x2,double y2,double g11,double g12,double g22)
{
 double d;
 d=g11*(x2-x1)*(x2-x1)+2*g12*(x2-x1)*(y2-y1)+g22*(y2-y1)*(y2-y1);
 d=sqrt(d);
 return d;   
    
}

VectorXd gradEbendn(const vector<vector<int>> &springpairs, const vector<spring> &springlist, const VectorXd &XY,
                    double g11,double g12, double g22, double kappa)
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

    //double dacos=-1-c*c/2-3*c*c*c*c/8;
    //if(GRADc.dot(GRADc)<.1) GRADc=1/s*GRADc;
    secondpart=secondpart+(1.0/(d12+d23))*2*(pi-acos(c))*s*GRADc;
  
  }

  grad=kappa*(firstpart+secondpart);
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

double quad(double x)
{
 double out=x*x-2*x-1;   
 return out;
}    
    






