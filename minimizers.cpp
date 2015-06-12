#include "random.h"
#include <algorithm>
#include "makemikadonetwork.h"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include <vector>
#include "EnergyandGradients.h"
#include "BendingGrad.h"
#include<iostream>
#include "minimizers.h"

using namespace Eigen;
using namespace std;


// void doFalsePosition(double &a1,double &a2,double &root,
//                     const VectorXd &XY,
//                     const VectorXd &s0, 
//                     const vector<spring> &springlist,
//                     const vector<vector<int>> &springpairs, 
//                     double kappa,
//                     const double g11,
//                     const double g12,
//                     const double g22)
// {
//  double fl,fh,xl,xh,swap,dx,del,f;
//  double xacc=.00001;
//  
//  int Maxit=100;
//  fl=dEda(XY+a1*s0,s0,springlist,springpairs,kappa,g11,g12,g22);
//  fh=dEda(XY+a2*s0,s0,springlist,springpairs,kappa,g11,g12,g22);
// 
//  
//  if(fl<0.0){  //xl =xlow and xh=xhigh --> f(xl)<f(xh);
//      xl=a1; 
//      xh=a2;
//  }
//  else{
//      xl=a2;
//      xh=a1;
//      swap=fl;
//      fl=fh;
//      fh=swap;
// }
// dx=xh-xl;
// //int ii=0;
// for(int i=0;i<Maxit;i++){
//     root=xl+dx*fl/(fl-fh); //This is a secant step
//     f=dEda(XY+root*s0,s0,springlist,springpairs,kappa,g11,g12,g22);
//     if(f<0.0){
//      del=xl-root;
//      xl=root;
//      fl=f;
//     }
//     else{
//         del=xh-root;
//         xh=root;
//         fh=f;
//     }
//     dx=xh-xl;
// //    ii++;
//     if(fabs(del)<xacc || f==0.0) break;
// }
// 
// //  cout<<"i=    "<<ii<<endl;
// }

void doSecant(double &root,
              const VectorXd &XY,
              const VectorXd &s0, 
              const vector<spring> &springlist,
              const vector<vector<int>> &springpairs, 
              double kappa,
              const double g11,
              const double g12,
              const double g22)
{
 double an2=0;
 double an1=1e-5;
 
 double an;
 double tol=1e-11;
 int q=0; 
 double dEda2,dEda1;
 
dEda2=dEda(XY+an2*s0,s0,springlist,springpairs,kappa,g11,g12,g22);

 do{ 
    dEda1=dEda(XY+an1*s0,s0,springlist,springpairs,kappa,g11,g12,g22);
    an=an1-dEda1*(an1-an2)/(dEda1-dEda2);
    an2=an1;
    an1=an;
    dEda2=dEda1;
    q++;
 }while(q<1000 && fabs(an2-an1)>tol);
 root=an;
 //cout<<"Q IS  "<<"\t"<<q<<endl;
}

void doConjStep(VectorXd &XY,
                VectorXd &s0,
                VectorXd &gradE,
                const vector<spring> &springlist,
                const vector<vector<int>> &springpairs,
                double kappa,
                int conjsteps,
                double g11,
                double g12,
                double g22)
{
    double a1=0.0;
    double a2=1.0;
    double betan;
    VectorXd gradEn(gradE.size());
    VectorXd sn(s0.size());
    functor network(XY,s0,springlist,springpairs,kappa,g11,g12,g22);
   doBracketfind(a1,a2,network);
    if(doBracketfind(a1,a2,network))
    {
        // double an=doFalsePosition(network,a1,a2);
        //double an=Ridder(network,a1,a2);
        double an=Brent(network,a1,a2,1e-12);
        //Update the positions.
        XY=XY+an*s0;
        gradEn=HarmonicGradient(springlist,XY,g11,g12,g22)+BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22);
        betan=(gradEn-gradE).dot(gradEn)/(gradE.dot(gradE));
        //Did not find bracket, reset CG-method    
    } else{
        betan=0.0;
        gradE=HarmonicGradient(springlist,XY,g11,g12,g22)+BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22);
        s0=-gradE;
        //cout<<"Bracket failed, Reset CG"<<endl;
        return;
    }
    if(conjsteps%5 ==0.0) betan=0.0;
    if(betan<0.0) betan=0; //max(betan,0)
    if(abs(gradEn.dot(gradE))>.5*gradE.dot(gradE)) betan=0.0; 
    if(-2*gradE.dot(gradE)>gradE.dot(s0) && gradE.dot(s0) >-.2*gradE.dot(gradE)) betan=0.0;
    //cout<<"\r"<<"** Beta "<<betan<<flush;
    sn=-gradEn+betan*s0;    
    gradE=gradEn;
    s0=sn;
}









