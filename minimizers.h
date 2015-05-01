#ifndef MINIMIZERS_H
#define MINIMIZERS_H

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define UNUSED (-1.11e30);
#define ITMAX 400
#define EPS 2.0e-8

#include<vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include "BendingGrad.h"
#include "EnergyandGradients.h"
//#include "nrutil.h"

struct functor{
            const Eigen::VectorXd &XY;
            const Eigen::VectorXd &s0;
            const std::vector<spring> &springlist;
            const std::vector<std::vector<int>> &springpairs;
            const double &kappa;
            const double &g11;
            const double &g12;
            const double &g22;
            
            functor(const Eigen::VectorXd &XY,
                const Eigen::VectorXd &s0,
                const std::vector<spring> &springlist,
                const std::vector<std::vector<int>> &springpairs,
                double kappa,
                const double g11,
                const double g12,
                const double g22) : XY(XY), s0(s0), springlist(springlist), springpairs(springpairs), 
                                    kappa(kappa), g11(g11), g12(g12), g22(g22){}
            double eval( double alpha )
            {
                return dEda(XY+alpha*s0,s0,springlist,springpairs,kappa,g11,g12,g22);
            }
};


template <typename struct_type>
inline int doBracketfind(double &a1,double &a2,struct_type functor)
                 
//This function finds the inteval on which a mathematical function passes through zero.
//that is [x1,x2] where f(x1)*f(x2)<0.0;
{
 int maxit=50;   
 double f1,f2,FACTOR;
 f1=functor.eval(a1);
 if(f1>0) return 0;
 f2=functor.eval(a2);

 FACTOR=1.618034;

 if(a1==a2){ //We need two different points
    cout<<"Bad initial range in bracketfinder"<<endl;
  }
 int jj=0;
  for(int j=0;j<maxit;j++){ //Make a bracket.
      if(f1*f2<0.0) break; 
    jj++;
        a2=a2+FACTOR*(a2-a1);
        f2=functor.eval(a2);
    if(j==49){
        cout<<"Bracket not found"<<endl; 
        return 0;
    }
  }  
   if(f1>f2) cout<<"EXTERMINATE  "<<functor.eval(0)<<endl;
   return 1;
}

template <typename struct_type>
inline double doFalsePosition(struct_type functor,double a1,double a2)
{
 double fl,fh,xl,xh,swap,dx,del,f, root;
 double xacc=1e-6;
 
 int Maxit=100;
 fl=functor.eval(a1);
 fh=functor.eval(a2);

 
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
    f=functor.eval(root);
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
return root;
//  cout<<"i=    "<<ii<<endl;
}

template <typename struct_type>
inline double Ridder(struct_type functor, double x1, double x2)
{
    int j;
    double ans, fh, fl,fm, fnew, s, xh, xl, xm,xnew;
    fl=functor.eval(x1);
    fh=functor.eval(x2);
    int maxit=60;
    double xacc=1e-12;

    if((fl>0.0 && fh<0.0) || (fl<0.0 && fh>0.0)){
        xl=x1;
        xh=x2;
        ans=UNUSED;

    for(j=1;j<=maxit;j++){
        std::cout<<j<<endl;
       xm=0.5*(xl+xh); 
       fm=functor.eval(xm);
       s=sqrt(fm*fm-fl*fh);
       
       if (s==0) return ans;
       xnew=xm+(xm-xl)*((fl>=fh? 1.0 :-1.0)*fm/s);
       if(fabs(xnew-ans)<=xacc) return ans;
       ans=xnew;
       fnew=functor.eval(ans);
       
       if(fnew==0.0) return ans;
       if( SIGN(fm,fnew) != fm){
        xl=xm;
        fl=fm;
        xh=ans;
        fh=fnew;
       } else if(SIGN(fh,fnew) !=fh){
           xl=ans;
           fl=fnew;
       } else cout<<"never get there"<<endl;
    if(fabs(xh-xl)<=xacc) return ans;
    }
   //cout<<"exeded count"<<endl;
 }
 else{
     if(fl==0.0) return x1;
     if(fh==0.0) return x2;
    }
    
    return 0.0;
    
}


template <typename struct_type>
inline double Brent(struct_type functor,double x1, double x2, double tol)
{
    int iter;
    double a=x1; 
    double b=x2; 
    double c=x2; 
    double d, e, min1, min2;
    double fa=functor.eval(a);
    double fb=functor.eval(b); 
    double fc, p, q, r, s, tol1, xm;
    
    if((fa>0.0 && fb >0.0) || (fa <0.0 && fb <0.0))
        std::cout<<"root must be bracketed"<<endl;
    
    fc=fb;
    for(iter =1;iter<ITMAX;iter++) {
        if((fb>0.0 && fc>0.0) || (fb <0.0 && fc<0.0)){
            c=a;
           fc=fa;
           e=d=b-a;
        }
        
        if(fabs(fc) <fabs(fb)){
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        
        tol1=2.0*EPS*fabs(b)+0.5*tol;
        xm=0.5*(c-b);
        if(fabs(xm) <=tol1 ||fb==0.0) return b;
        if(fabs(e)>=tol1 && fabs(fa)>fabs(fb)) {
            s=fb/fa;
            if(a==c){
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if(p>0.0) q=-q;
            p=fabs(p);
            min1=3.0*xm*q-fabs(tol1*q);
            min2=fabs(e*q);
            if(2.0*p<(min1<min2 ? min1 : min2)){
                e=d;
                d=p/q;
            } else{
                d=xm;
                e=d;
            }
        }else {
            d=xm;
            e=d;
        }
        a=b;
        fa=fb;
        if(fabs(d)>tol1)
            b+=d;
        else
            b+=SIGN(tol1,xm);
        fb=functor.eval(b);
    }

   // nrerror("Maximum number of iterations exceeded in zbrent");
    return 0.0;
};


int doBracketfind(double &a1,
                  double &a2,
                   const Eigen::VectorXd &XY,
                   const Eigen::VectorXd &s0,
                   const std::vector<spring> &springlist,
                   const std::vector<std::vector<int>> &springpairs, 
                   double kappa,
                   const double g11,
                   const double g12,
                   const double g22);


void doFalsePosition(double &xl,double &xh,double &root,
                 const Eigen::VectorXd &XY,
                 const Eigen::VectorXd &s0,
                 const std::vector<spring> &springlist,
                 const std::vector<std::vector<int>> &springpairs, 
                 double kappa,
                 const double g11,
                 const double g12,
                 const double g22);


void doSecant(   double &root,
                 const Eigen::VectorXd &XY,
                 const Eigen::VectorXd &s0,
                 const std::vector<spring> &springlist,
                 const std::vector<std::vector<int>> &springpairs, 
                 double kappa,
                 const double g11,
                 const double g12,
                 const double g22);


void doConjStep(Eigen::VectorXd &XY,
                Eigen::VectorXd &s0,
                Eigen::VectorXd &gradE,
                const std::vector<spring> &springlist,
                const std::vector<std::vector<int>> &springpairs,
                double kappa,
                int conjstep,
                double g11,
                double g12,
                double g22);



#endif // MINIMIZERS_H