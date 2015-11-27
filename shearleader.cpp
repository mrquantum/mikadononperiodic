//In this header file we put all the shear protocol
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
//#include "exportfiles.h"
#include "writefunctions.h"
#include "cgwoagony.h"
#include <nlopt.h>
#include "cgmethod.h"

using namespace std;
using namespace Eigen;

void shearsteps(double deltaboxdx,
               int NumberStepsRight,
               int NumberStepsLeft,
               vector<spring> &springlist,
               vector<vector<int>> &springpairs,
               VectorXd &XY,
                int bendingOn,
               double kappa,
               int Nit,
               double tolGradE,
               ofstream &shearcoordinates,
               ofstream &shearenergy)
              {
  
    double g11,g12,g22; //The components of the metric tensor
  
    double boxdx=0.0;
    double EBEND,ESTRETCH,ETOT;
    VectorXd gradE(XY.size());
    VectorXd s0(XY.size());
    int conjsteps;
    double lenGrad;
    int iter;
    double fret;
    networkinfo info;
    info.g11=1.0;
    info.g12=0.0;
    info.g22=1.0;
    info.springlist=springlist;
    info.springpairs=springpairs;
    info.size=XY.size();
    //first do a calibtration
    frprmn(XY.data(),XY.size(),1.0e-10,&iter,&fret,EnergyNetworkn,HarmonicGradientn,info);
    
    //CGAGONY(XY,springlist,springpairs,bendingOn,kappa,1.0,0.0,1.0);
    //then shake
    //XY=XY+shake(XY.size(),0.000001);
    //minimize again;
    //CGAGONY(XY,springlist,springpairs,bendingOn,kappa,1.0,0.0,1.0);

    for(int k=0;k<(NumberStepsLeft+NumberStepsRight);k++){
        g11=1.0;
        g12=boxdx;
        g22=1.0+boxdx*boxdx;
        info.g11=1.0;
        info.g12=boxdx;
        info.g22=1.0+boxdx*boxdx;
        //The right CG algorithm, works only w. harmonic springs now
        cout<<"length before Harm\t"<<sqrt(HarmonicGradient(springlist,XY,g11,g12,g22).dot(HarmonicGradient(springlist,XY,g11,g12,g22)))<<endl;
        cout<<"length before Bend\t"<<sqrt(BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22)
            .dot(BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22)))<<endl;
        
        frprmn(XY.data(),XY.size(),1.0e-6,&iter,&fret,EnergyNetworkn,HarmonicGradientn,info);
        cout<<"Es= "<<fret<<endl;

           cout<<"length after Harm\t"<<HarmonicGradient(springlist,XY,g11,g12,g22).dot(HarmonicGradient(springlist,XY,g11,g12,g22))<<endl;
           cout<<"length after Bend\t"<<BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22).dot(BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22))<<endl;
        
        //cout<<"lengrad="<<lenGrad<<endl;
        Write_ShearCoordinates_2txt(shearcoordinates,XY);
        Write_ShearEnergy_2txt(shearenergy,boxdx,fret,ESTRETCH,EBEND,lenGrad,conjsteps);

        if(k<NumberStepsRight){
            boxdx=boxdx+deltaboxdx;
        } else {
            boxdx=boxdx-deltaboxdx;
        }
        cout<<"This was a shearstep     "<<k<<endl;
    }
}


VectorXd shake(int size, double Temperature){
    my_random::get_gre(2);
    VectorXd dXY(size);
    for(int i=0 ;i<size;        i++){
        dXY(i)=Temperature*(randf()-0.5);
    }
    
    return dXY;
    
}

double EnergyNLopt(unsigned n, const double *XY, double *gradE, void *my_func_data)
{       
   VectorXd XYEigen(n);
   for(int i=0  ;i<n;i++){
       XYEigen(i)=XY[i];
   }
   
   functor *F = static_cast<functor*>(my_func_data);
   
   if( gradE ){
        VectorXd Gr=HarmonicGradient(F->springlist,XYEigen,F->g11,F->g12,F->g22);
   
        for(int i=0;i<n;i++){
            gradE[i]=Gr(i);
        }
   }
   double Energy=Energynetwork(F->springlist,XYEigen,F->g11,F->g12,F->g22);
   return Energy;
   
}





