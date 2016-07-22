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
#include "newbendinggrad.h"

using namespace std;
using namespace Eigen;

double Efunc2use(double *p,networkinfo parameters){
  double E;
//   cout<<"test ef"<<endl;
  int bendingon=parameters.bendingon;
  if(bendingon==0){
      E=EnergyNetworkn(p,parameters);
  } else if(bendingon ==1){
//      cout<<"test ef"<<endl;

      E=EnergyNetworkn(p,parameters)+EbendingCn(p,parameters);

//         cout<<"test ef"<<endl;

}
  return E;
}

void grad2use(double *p,double *xi, networkinfo parameters)
{
//   cout<<"test g222"<<endl;
  int bendingOn=parameters.bendingon;
  int size=parameters.size;
//         cout << parameters.springlist.size() << endl;

//   vector<spring> springlist=parameters.springlist;
//   vector<vector<int > > springpairs=parameters.springpairs;
  double sheardeformation=parameters.sheardeformation;
  double kappa=parameters.kappa;
  //Setting the gradient to 0
  for(int i=0;i<size;i++){
      xi[i]=0.0;
  }
  
  if(bendingOn==0){
    HarmonicGradPhys(parameters.springlist,p,xi,size,sheardeformation);
  }else if(bendingOn==1){
//       cout << parameters.springlist.size() << endl;
    HarmonicGradPhys(parameters.springlist,p,xi,size,sheardeformation);
//         cout<<"test g"<<endl;

      double xb[size]; //define bendinggrad and set to 0.0
      for(int j=0;j<size;j++){
          xb[j]=0.0;
      }

      bendinggradnew(p,xb,size,parameters.springlist,parameters.springpairs,kappa,sheardeformation);
//         cout<<"test g"<<endl;

      //Contruct the total gradient.
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
               int Nit,
               double tolGradE,
               ofstream &shearcoordinates,
               ofstream &shearenergy,
               ofstream &Strestens,
               ofstream &EN)
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
    
    networkinfo info(springlist,springpairs);
    info.g11=1.0;
    info.g12=0.0;
    info.g22=1.0;
//     info.springlist=springlist;
//     info.springpairs=springpairs;
    info.size=XY.size();
    info.bendingon=bendingOn;
    info.sheardeformation=boxdx;
    info.kappa=kappa;

    //XY is in BOX coordinates 
   frprmn(XY.data(),XY.size(),1e-20,&iter,&fret,Efunc2use,grad2use,info,EN); //Do one calibration before shearing
//     frprmn(XY.data(),XY.size(),1e-12,&iter,&fret,NULL,NULL,info); //Do one calibration before shearing
//     CGAGONY(XY,springlist,springpairs,bendingOn,kappa,1.0,0.0,1.0,0.0);
    stresstensor S=StressTensor(springlist,XY,0.0);  
    Strestens<<boxdx<<"\t"<<S.sxx<<"\t"<<S.sxy<<"\t"<<S.syy<<"\t"<<0.5*(S.sxx+S.syy)<<endl;

    for(int k=0;k<(NumberStepsLeft+NumberStepsRight);k++){
        g11=1.0;
        g12=boxdx;
        g22=1.0+boxdx*boxdx;
        info.g11=1.0;
        info.g12=boxdx;
        info.g22=1.0+boxdx*boxdx;
        info.sheardeformation=boxdx;
	
//         CGAGONY(XY,springlist,springpairs,bendingOn,kappa,g11,g12,g22,boxdx); //fast
        frprmn(XY.data(),XY.size(),1e-15,&iter,&fret,Efunc2use,grad2use,info,EN); //slow but better THIS ONE

        S=StressTensor(springlist,XY,boxdx);
        Write_ShearCoordinates_2txt(shearcoordinates,XY);
	//ESTRETCH=EnergyNetworkn(XY.data(),info);
	ESTRETCH=StretchEnergy(springlist,XY,boxdx);
//         cout << "Now EBEND: " << XY(600) << "\t" << kappa << "\t" << boxdx << "\t" << springpairs[80][0] << "\t" << springlist[80].one <<  endl;
        EBEND=BendEnergy(springpairs,springlist,XY,kappa,boxdx);
//         cout << "RESULT: " << EBEND << endl;

        Write_ShearEnergy_2txt(shearenergy,boxdx,ESTRETCH+EBEND,ESTRETCH,EBEND,lenGrad,iter);

        if(k<NumberStepsRight){
            boxdx=boxdx+deltaboxdx;
        } else {
            boxdx=boxdx-deltaboxdx;
        }
        cout<<"This was a shearstep     "<<k<<" with    "<<iter<<endl;
    }
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





