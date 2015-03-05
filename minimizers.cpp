#include "random.h"
#include <algorithm>
#include "makemikadonetwork.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <vector>
#include "EnergyandGradients.h"
#include "BendingGrad.h"
#include<iostream>
using namespace Eigen;
using namespace std;

int doBracketfind(double &a1,double &a2,
                   const VectorXd &XY,
                   const VectorXd &s0, 
                   const vector<spring> &springlist,
                   const vector<vector<int>> &springpairs, 
                   double kappa,
                   const double g11,
                   const double g12,
                   const double g22)
                 
//This function finds the inteval on which a mathematical function passes through zero.
//that is [x1,x2] where f(x1)*f(x2)<0.0;
{
 int maxit=50;   
 double f1,f2,FACTOR;
 f1=dEda(XY+a1*s0,s0,springlist,springpairs,kappa,g11,g12,g22);
 
 if(f1>0) return 0;
 
 f2=dEda(XY+a2*s0,s0,springlist,springpairs,kappa,g11,g12,g22);

 FACTOR=1.6;
  
 if(a1==a2){ //We need two different points
    cout<<"Bad initial range in bracketfinder"<<endl;
  }
 int jj=0;
  for(int j=0;j<maxit;j++){ //Make a bracket.
      if(f1*f2<0.0) break; 
    jj++;
//     if(abs(f1)<abs(f2)){
//         a1=a1+FACTOR*(a1-a2);
//         f1=dEda(XY+a1*s0,s0,springlist,springpairs,kappa);
//     }
//     else{
        a2=a2+FACTOR*(a2-a1);
        f2=dEda(XY+a2*s0,s0,springlist,springpairs,kappa,g11,g12,g22);
    if(j==49){
        cout<<"not found"<<endl; 
        return 0;
//     }
        
    }
  }  
   
//   cout<<jj<<endl;
   if(f1>f2) cout<<"EXTERMINATE  "<<dEda(XY,s0,springlist,springpairs,kappa,g11,g12,g22)<<endl;
   return 1;
}


void doFalsePosition(double &a1,double &a2,double &root,
                    const VectorXd &XY,
                    const VectorXd &s0, 
                    const vector<spring> &springlist,
                    const vector<vector<int>> &springpairs, 
                    double kappa,
                    const double g11,
                    const double g12,
                    const double g22)
{
 double fl,fh,xl,xh,swap,dx,del,f;
 double xacc=.00001;
 
 int Maxit=100;
 fl=dEda(XY+a1*s0,s0,springlist,springpairs,kappa,g11,g12,g22);
 fh=dEda(XY+a2*s0,s0,springlist,springpairs,kappa,g11,g12,g22);

 
 if(fl<0.0){  //xl =xlow and xh=xhigh --> f(xl)<f(xh);
     xl=a1; 
     xh=a2;
 }
 else{
     xl=a2;
     xh=a1;
     swap=fl;
     fl=fh;
     fh=swap;
}
dx=xh-xl;
//int ii=0;
for(int i=0;i<Maxit;i++){
    root=xl+dx*fl/(fl-fh); //This is a secant step
    f=dEda(XY+root*s0,s0,springlist,springpairs,kappa,g11,g12,g22);
    if(f<0.0){
     del=xl-root;
     xl=root;
     fl=f;
    }
    else{
        del=xh-root;
        xh=root;
        fh=f;
    }
    dx=xh-xl;
//    ii++;
    if(fabs(del)<xacc || f==0.0) break;
}

//  cout<<"i=    "<<ii<<endl;
}

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
 double an1=0.00001;
 
 double an;
 double tol=0.000000000001;
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
 cout<<"Q IS  "<<"\t"<<q<<endl;
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
    double a2=.0001;
    double betan;
    double root=0.0;
    VectorXd gradEn(gradE.size());
    VectorXd sn(s0.size());

    //Did find bracket
    doBracketfind(a1,a2,XY,s0,springlist,springpairs,kappa,g11,g12,g22);
    if(doBracketfind(a1,a2,XY,s0,springlist,springpairs,kappa,g11,g12,g22)){
    doFalsePosition(a1,a2,root,XY,s0,springlist,springpairs,kappa,g11,g12,g22);
    //cout<<root<<endl;
    //    doSecant(root,XY,s0,springlist,springpairs,kappa,g11,g12,g22); //Do Linesearch;
        double an=root;
        XY=XY+an*s0; //Update positions
        gradEn=HarmonicGradient(springlist,XY,g11,g12,g22)+BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22);
        betan=(gradEn-gradE).dot(gradEn)/(gradE.dot(gradE));
    
    //Did not find bracket, reset CG-method    
    } else{
        betan=0;
        gradE=HarmonicGradient(springlist,XY,g11,g12,g22)+BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22);
        s0=-gradE;
        cout<<"Bracket failed, Reset CG"<<endl;
        return;
    }
    if(conjsteps%100 ==0) betan=0;
    if(betan<0) betan=0; //max(betan,0)
    if(abs(gradEn.dot(gradE))>.5*gradE.dot(gradE)) betan=0; 
    if(-2*gradE.dot(gradE)>gradE.dot(s0) && gradE.dot(s0) >-.2*gradE.dot(gradE)) betan=0;
    //cout<<"\r"<<"** Beta "<<betan<<flush;
    sn=-gradEn+betan*s0;    
    gradE=gradEn;
    s0=sn;
}









