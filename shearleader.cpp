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

using namespace std;
using namespace Eigen;

void shearsteps(double deltaboxdx,
               int NumberStepsRight,
               int NumberStepsLeft,
               vector<spring> &springlist,
               vector<vector<int>> &springpairs,
               VectorXd &XY,
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
    
    for(int k=0;k<(NumberStepsLeft+NumberStepsRight);k++){
        g11=1.0;
        g12=boxdx;
        g22=1.0+boxdx+boxdx;
        //Check wheather there is bending rigidity or not. If not, we shall not
        //have to compute the bending gradient and multiply it by zero each step
        //which will save massive amounts of time (hopefully)
        if(kappa==0.0){
            ESTRETCH=Energynetwork(springlist,XY,g11,g12,g22);
            ETOT=ESTRETCH;
            gradE=HarmonicGradient(springlist,XY,g11,g12,g22);
            cout<<"kappa was zero and we are in the first case"<<endl;
        } else{
            ESTRETCH=Energynetwork(springlist,XY,g11,g12,g22);
            EBEND=EbendingC(springpairs,springlist,XY,kappa,g11,g12,g22);
            ETOT=ESTRETCH+EBEND;
            gradE=HarmonicGradient(springlist,XY,g11,g12,g22)+BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22);
            cout<<"kappa was not zero and we are in the second case"<<endl;
        }
        s0=-gradE;
  
        //Now we initialize the nonlin. conjgrad-method
        //Maybe we are going to change this in the future 
        do{
            conjsteps++;
            doConjStep(XY,s0,gradE,springlist,springpairs,kappa,conjsteps,g11,g12,g22);
            if(kappa==0.0){
                ESTRETCH=Energynetwork(springlist,XY,g11,g12,g22);
                ETOT=ESTRETCH;
            } else{
                ESTRETCH=Energynetwork(springlist,XY,g11,g12,g22);
                EBEND=EbendingC(springpairs,springlist,XY,kappa,g11,g12,g22);  
                ETOT=ESTRETCH+EBEND;
            }
            lenGrad=sqrt(gradE.dot(gradE));
        }while(conjsteps<Nit && lenGrad>tolGradE);
 
        Write_ShearCoordinates_2txt(shearcoordinates,XY);
        Write_ShearEnergy_2txt(shearenergy,boxdx,ETOT,ESTRETCH,EBEND,lenGrad,conjsteps);

        if(k<NumberStepsRight){
            boxdx=boxdx+deltaboxdx;
        } else {
            boxdx=boxdx-deltaboxdx;
        }
        cout<<"This was a shearstep     "<<k<<flush;
    }
}