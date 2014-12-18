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

 FACTOR=1.1;
  
 if(a1==a2){ //We need two different points
    cout<<"Bad initial range in bracketfinder"<<endl;
  }
 
  for(int j=0;j<maxit;j++){ //Make a bracket.
      if(f1*f2<0.0) break; 
    
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

   if(f1>f2) cout<<"EXTERMINATE  "<<dEda(XY,s0,springlist,springpairs,kappa,g11,g12,g22)<<endl;
   return 1;
}
    
    
    
    
void doBisection(double a1,double a2,double &root,
    const VectorXd &XY,
    const VectorXd &s0, 
    const vector<spring> &springlist,
    const vector<vector<int>> &springpairs, 
    double kappa,
    const double g11,
    const double g12,
    const double g22
    )
{
 double f1,f2,fc,c;
 int q=0;
 do{
 c=0.5*(a1+a2);
 f1=dEda(XY+a1*s0,s0,springlist,springpairs,kappa,g11,g12,g22);
 f2=dEda(XY+a2*s0,s0,springlist,springpairs,kappa,g11,g12,g22);
 fc=dEda(XY+c*s0,s0,springlist,springpairs,kappa,g11,g12,g22);
 
 if(f1*fc>0.0) a1=c;
 else a2=c;
 q++;
}while(abs(fc)>.000001 && q<50);
 root=c;
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
 double xacc=.0000000000001;
 
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
    if(abs(del)<xacc || f==0.0) break;
}
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
 double an2=-0.00000001;
 double an1=0.0000001;
 double an;
 double tol=0.0000001;
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
 }while(q<50 && abs(an2-an1)>tol);   
 root=an;   
}

void doConjStep(VectorXd &XY,
                VectorXd &s0,
                VectorXd &gradE,
                vector<spring> &springlist,
                vector<vector<int>> &springpairs,
                double &root,
                double kappa,
                int conjsteps,
                const double g11,
                const double g12,
                const double g22)
{
  
//double a1=-.000001;
    double a1=0.0;
    double a2=.000000001;
    double betan;
    VectorXd gradEn(gradE.size());
    VectorXd sn(s0.size());
    
    

    //Did find breacket
    if(doBracketfind(a1,a2,XY,s0,springlist,springpairs,kappa,g11,g12,g22)){
    cout<<a1<<"\t"<<a2<<"\t"<<dEda(XY+a1*s0,s0,springlist,springpairs,kappa,g11,g12,g22)<<"\t"<<
    dEda(XY+a2*s0,s0,springlist,springpairs,kappa,g11,g12,g22)<< endl;   
    
    doFalsePosition(a1,a2,root,XY,s0,springlist,springpairs,kappa,g11,g12,g22);
    //doSecant(root,XY,s0,springlist,springpairs,kappa,g11,g12,g22); //Do Linesearch;
    double an=root;
    XY=XY+an*s0; //Update positions
    
    gradEn=Gradient(springlist,XY,g11,g12,g22);//+gradEbend(springpairs,springlist,XY,kappa);
    //double betan=gradEn.dot(gradEn)/(gradE.dot(gradE));
    betan=(gradEn-gradE).dot(gradEn)/(gradE.dot(gradE));
    
    //Did not find bracket, reset CG-method    
    } else{
        betan=0;
        cout<<"Bracket failed, Reset CG"<<endl;
    }
    if(betan<0.0000) betan=0;   //beta is max(beta,0)
    if(conjsteps%100 ==0) betan=0;
    if(abs(gradEn.dot(gradE))>.3*gradE.dot(gradE)) betan=0; 
    //if(-2*gradE.dot(gradE)>gradE.dot(s0) && gradE.dot(s0) >-.2*gradE.dot(gradE)) betan=0;
    cout<<"**  "<<betan<<endl;
    
    sn=-gradEn+betan*s0;    
    gradE=gradEn;
    s0=sn;
}

void doSteepestDescent(VectorXd &XY,
                VectorXd &s0,
                VectorXd &gradE,
                vector<spring> &springlist,
                vector<vector<int>> &springpairs,
                double &root,
                double kappa,
                const double g11,
                const double g12,
                const double g22
                //VectorXd &b
                      )
{
    doSecant(root,XY,s0,springlist,springpairs,kappa,g11,g12,g22);
    double an=root;
    XY=XY+an*s0;
//    s0=-gradE;
    s0=-Gradient(springlist,XY,g11,g12,g22);
    //-gradEbend(springpairs,springlist,XY,kappa);
    gradE=-s0;   
       
}

VectorXd Hessianapprox(const VectorXd &XY,const VectorXd &XYm1,const VectorXd &g0,const VectorXd &g0m1)
{
    VectorXd H(XY.size());
    VectorXd dXY(XY.size());
    VectorXd dg(g0.size());
        for(int i=0;i<XY.size();i++){
            dXY(i)=XY(i)-XYm1(i);
            dg(i)=g0(i)-g0m1(i);
            double hi=dg(i)/dXY(i);
            H(i)=hi > 0 ? hi : 1;
         }
    
return H; 
   
}