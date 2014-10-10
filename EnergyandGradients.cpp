#include "random.h"
#include <algorithm>
#include "makemikadonetwork.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <vector>
using namespace Eigen;
using namespace std;
const double pi=4.0*atan(1.0);
   
  
double Energynetwork(const vector<spring> &springlist, const VectorXd &XY,double k,double L)
{
  VectorXd X(XY.size()/2);
  VectorXd Y(XY.size()/2);
  X=XY.head(XY.size()/2);
  Y=XY.tail(XY.size()/2);
  double Energy=0;
    for(int i=0;i<springlist.size();i++){
      Energy=Energy+
      0.5*k*pow(sqrt(pow(X(springlist[i].one)-(X(springlist[i].two)+springlist[i].wlr),2)
      +pow(Y(springlist[i].one)-(Y(springlist[i].two)+springlist[i].wud),2))-L,2);
    }
  return Energy;  
}
  
  
VectorXd Gradient(const vector<spring> &springlist,const VectorXd &XY,double k, double L)
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
 gradE(springlist[i].one)=gradE(springlist[i].one)+
 k*(sqrt(pow(X(springlist[i].one)-(X(springlist[i].two)+springlist[i].wlr),2)
 +pow(Y(springlist[i].one)-(Y(springlist[i].two)+springlist[i].wud),2))-L)*
 (X(springlist[i].one)-(X(springlist[i].two)+springlist[i].wlr))/
 sqrt(pow(X(springlist[i].one)-(X(springlist[i].two)+springlist[i].wlr),2)
 +pow(Y(springlist[i].one)-(Y(springlist[i].two)+springlist[i].wud),2));
 
 
 gradE(springlist[i].two)=gradE(springlist[i].two)-
 k*(sqrt(pow(X(springlist[i].one)-(X(springlist[i].two)+springlist[i].wlr),2)
 +pow(Y(springlist[i].one)-(Y(springlist[i].two)+springlist[i].wud),2))-L)*
 (X(springlist[i].one)-(X(springlist[i].two)+springlist[i].wlr))/
 sqrt(pow(X(springlist[i].one)-(X(springlist[i].two)+springlist[i].wlr),2)
 +pow(Y(springlist[i].one)-(Y(springlist[i].two)+springlist[i].wud),2));

 
 gradE(springlist[i].one+X.size())=gradE(springlist[i].one+X.size())+
  k*(sqrt(pow(X(springlist[i].one)-(X(springlist[i].two)+springlist[i].wlr),2)
 +pow(Y(springlist[i].one)-(Y(springlist[i].two)+springlist[i].wud),2))-L)*
 (Y(springlist[i].one)-(Y(springlist[i].two)+springlist[i].wud))/
 sqrt(pow(X(springlist[i].one)-(X(springlist[i].two)+springlist[i].wlr),2)
 +pow(Y(springlist[i].one)-(Y(springlist[i].two)+springlist[i].wud),2));
 
 
  gradE(springlist[i].two+X.size())=gradE(springlist[i].two+X.size())-
  k*(sqrt(pow(X(springlist[i].one)-(X(springlist[i].two)+springlist[i].wlr),2)
 +pow(Y(springlist[i].one)-(Y(springlist[i].two)+springlist[i].wud),2))-L)*
 (Y(springlist[i].one)-(Y(springlist[i].two)+springlist[i].wud))/
 sqrt(pow(X(springlist[i].one)-(X(springlist[i].two)+springlist[i].wlr),2)
 +pow(Y(springlist[i].one)-(Y(springlist[i].two)+springlist[i].wud),2));
}

return gradE;  
}  


double dEda(const VectorXd &XY,const VectorXd &s0,const vector<spring> &springlist,double k,double L) 
{  
  double dEda=s0.dot(Gradient(springlist,XY,k,L));
  return dEda;
}

