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
//#include <nlopt.h>
#include "cgmethod.h"
#include "stresstensor.h"
#include "simpleBendingGrad.h"

using namespace std;
using namespace Eigen;

double Efunc2use(double *p,networkinfo parameters){
  double E;
  int bendingon=parameters.bendingon;
  if(bendingon==0){
      E=EnergyNetworkn(p,parameters);
  } else if(bendingon ==1){
      E=EnergyNetworkn(p,parameters)+EbendingCn(p,parameters);
  }
  return E;
}

void grad2use(double *p,double *xi, networkinfo parameters)
{
  int bendingOn=parameters.bendingon;
  if(bendingOn==0){
      HarmonicGradientn(p,xi,parameters);
  }else if(bendingOn==1){
      HarmonicGradientn(p,xi,parameters);
      int size=parameters.size;
      double xb[size];
      BendinggradN(p,xb,parameters);
      for(int j=0;j<size;j++){
	xi[j]=xi[j]+xb[j];
      }
  }
}

void shearsteps(double deltaboxdx,
               int NumberStepsRight,
               int NumberStepsLeft,
               vector<spring> &springlist,
               vector<vector<int>> &springpairs,
               VectorXd &XY,
                int bendingOn,
               double kappa,
               VectorXd effkappa, 
               int Nit,
               double tolGradE,
               ofstream &shearcoordinates,
               ofstream &shearenergy,
               ofstream &Strestens)
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
    info.bendingon=bendingOn;
    info.sheardeformation=boxdx;
    info.effkappa=effkappa;
    frprmn(XY.data(),XY.size(),1e-20,&iter,&fret,Efunc2use,grad2use,info); //Do one calibration before shearing
    //CGAGONY(XY,springlist,springpairs,bendingOn,kappa,1.0,0.0,1.0);
    stresstensor S=StressTensor(springlist,XY,0.0);  
    //cout<<"sigma_xx="<<S.sxx<<"\t sigma_xy="<<S.sxy<<"\t sigma_yy="<<S.syy<<endl;
    Strestens<<boxdx<<"\t"<<S.sxx<<"\t"<<S.sxy<<"\t"<<S.syy<<"\t"<<0.5*(S.sxx+S.syy)<<endl;

    for(int k=0;k<(NumberStepsLeft+NumberStepsRight);k++){
        g11=1.0;
        g12=boxdx;
        g22=1.0+boxdx*boxdx;
        info.g11=1.0;
        info.g12=boxdx;
        info.g22=1.0+boxdx*boxdx;
        info.sheardeformation=boxdx;
	
        //CGAGONY(XY,springlist,springpairs,bendingOn,kappa,g11,g12,g22,boxdx); //fast
        frprmn(XY.data(),XY.size(),1e-30,&iter,&fret,Efunc2use,grad2use,info); //slow but better THIS ONE
        S=StressTensor(springlist,XY,boxdx);

        Write_ShearCoordinates_2txt(shearcoordinates,XY);
	//ESTRETCH=EnergyNetworkn(XY.data(),info);
	ESTRETCH=StretchEnergy(springlist,XY,boxdx);
        EBEND=BendEnergy(springpairs,springlist,XY,kappa,boxdx);

        Write_ShearEnergy_2txt(shearenergy,boxdx,ESTRETCH+EBEND,ESTRETCH,EBEND,lenGrad,iter);

        if(k<NumberStepsRight){
            boxdx=boxdx+deltaboxdx;
        } else {
            boxdx=boxdx-deltaboxdx;
        }
        cout<<"This was a shearstep     "<<k<<" with    "<<iter<<endl;
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





