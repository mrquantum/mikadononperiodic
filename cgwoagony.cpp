#include <iostream>
#include <ctime>
#include "random.h"
#include "makemikadonetwork.h"
#include "EnergyandGradients.h"
#include "minimizers.h"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <math.h>
#include <functional>
#include "importparam.h"
#include "BendingEnergy.h"
#include "BendingGrad.h"
#include "clusters.h"
#include <stdio.h>
#include <stdlib.h>
#include "writefunctions.h"
#include "shearleader.h"
#include "simpleBendingGrad.h"

using namespace std;
using namespace Eigen;

void CGAGONY(VectorXd &XY,
             vector<spring> &springlist, 
             vector<vector<int>> &springpairs, 
             int bendingon, 
             double kappa, 
             double g11, double g12, double g22,double sheardeformation){
    int i=0;
    int k=0;
    int j;
    int brk; //break boolean
    int jmax=1000;
    int imax=1e6;
    VectorXd r(XY.size());
    VectorXd s,d;
    double deltanew,delta0,deltad,deltaold,deltamid;
    double eps=1e-10;
    double sigma0=.1;
    double alpha;
    double eta,etaprev,deta,beta;
    
    //Calculate the gradient
    if(bendingon==0){
        //r=HarmonicGradient(springlist,XY,g11,g12,g22);
        r=HarmonicGradPhys(springlist,XY,sheardeformation);
    } else {
        r=HarmonicGradient(springlist,XY,g11,g12,g22)+BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22);
    }
    
    s=r;
    d=s;
    
    deltanew=r.dot(d);
    deltaold=deltanew;
    alpha=-sigma0;

    //Linesearching
    while(i<imax && deltanew>eps*eps*delta0){
        j=0;
        sigma0=-alpha;
        deltad=d.dot(d);
        if(bendingon==0){
            //etaprev=(HarmonicGradient(springlist,XY+sigma0*d,g11,g12,g22)).dot(d);
            etaprev=(HarmonicGradPhys(springlist,XY+sigma0*d,sheardeformation)).dot(d);
        } else{
            etaprev=(HarmonicGradient(springlist,XY+sigma0*d,g11,g12,g22)+
            BendingGrad(springpairs,springlist,XY+sigma0*d,kappa,g11,g12,g22)).dot(d);
        }
        
        do{
            if(bendingon==0){
                //eta=(HarmonicGradient(springlist,XY,g11,g12,g22)).dot(d);
                eta=(HarmonicGradPhys(springlist,XY,sheardeformation)).dot(d);
             } else{
                eta=(HarmonicGradient(springlist,XY,g11,g12,g22)+
                BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22)).dot(d);
            }
            double dsq=d.dot(d);
            deta=etaprev-eta;
            if(fabs(deta)>1e-22){
                alpha=alpha*eta/(etaprev-eta);
                //cout<<alpha<<"      alpha"<<endl;
                XY=XY+alpha*d;
                //cout<<"L2= "<<XY.dot(XY)<<endl;
                etaprev=eta;
                j++;
                brk=0;
            } else{
                brk=1;
                break;
                
            }
        }while(brk==0 && (j<jmax && alpha*alpha*deltad>eps*eps));
        if(bendingon==0){
            //r=-HarmonicGradient(springlist,XY+sigma0*d,g11,g12,g22);
            r=-HarmonicGradPhys(springlist,XY+sigma0*d,sheardeformation);
        } else{
            r=-(HarmonicGradient(springlist,XY+sigma0*d,g11,g12,g22)+
            BendingGrad(springpairs,springlist,XY+sigma0*d,kappa,g11,g12,g22));
        }
        deltaold=deltanew;
        deltamid=r.dot(s);
        
        s=r;
        deltanew=r.dot(s);
        beta=(deltanew-deltamid)/deltaold;
        k++;
        
        if(brk==1){ //len grad doesnt change anymore
            break;
        }
        
        if(k==10 || beta<0){
            d=s;
            k=0;
        } else {
            d=s+beta*d;
        }
        i++;
    }
    cout<<i<<"  steps"<<endl;

}


void CGAGONY2(VectorXd &XY,
             vector<spring> &springlist, 
             vector<vector<int>> &springpairs, 
             int bendingon, 
             double kappa, 
             double g11, double g12, double g22){

    int i=0;
    int k=0;
    int j;
    int imax=1e6;
    double eps=1e-6;
    VectorXd r,d;
    double deltanew,delta0,deltad;
    double eta,etaprev;
    double alpha,beta;
    double sigma0=1e-8;
    
    if(bendingon==0){
        r=-HarmonicGradient(springlist,XY,g11,g12,g22);
    } else{
        r=-HarmonicGradient(springlist,XY,g11,g12,g22)-BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22);
    }
    
    d=r;
    deltanew=r.dot(r);
    delta0=deltanew;
    
    while(i<imax && deltanew>eps*eps*delta0){
        j=0;
        deltad=d.dot(d);
        alpha=-sigma0;
        
        if(bendingon==0){
            etaprev=(HarmonicGradient(springlist,XY+sigma0*d,g11,g12,g22)).dot(d);
        } else{
            etaprev=(HarmonicGradient(springlist,XY+sigma0*d,g11,g12,g22)+
            BendingGrad(springpairs,springlist,XY+sigma0*d,kappa,g11,g12,g22)).dot(d);
        }
    
        do{ //This is the secant method, maybe replace it w. some other line search/
            if(bendingon==0){
                eta=(HarmonicGradient(springlist,XY,g11,g12,g22)).dot(d);
            } else{
                eta=(HarmonicGradient(springlist,XY,g11,g12,g22)+
                BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22)).dot(d);
            }
            alpha=alpha*eta/(etaprev-eta);
            XY=XY+alpha*d;
            etaprev=eta;
            j++;
        }while(j<100 && alpha*alpha*deltad>eps*eps);
        
        //calculate new gradient
        if(bendingon==0){
            r=-HarmonicGradient(springlist,XY,g11,g12,g22);
        } else{
            r=-HarmonicGradient(springlist,XY,g11,g12,g22)-BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22);
        }
        
        delta0=deltanew;
        deltanew=r.dot(r);
        beta=delta0/deltanew;
        d=r+beta*d;
        k++;
        if(k==5 || r.dot(d)<0){
            d=r;
            k=0;
        }
        i++;
    }
    
}