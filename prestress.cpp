#include "random.h"
#include <algorithm>
#include "makemikadonetwork.h"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include <vector>
#include "EnergyandGradients.h"
#include "minimizers.h"
#include<iostream>

using namespace std;
using namespace Eigen;

double prestress(const vector<spring> &springlist,
                 const vector<vector<int>> &springpairs,
                 VectorXd &XY,
                 double kappa,int Nit)
{
        double eps=0.000001;
        double g11=1+2*eps;
        double g22=1+2*eps;
        double g12=0;
        double Est0=Energynetwork(springlist,XY,1,0,1);
        double Eben0=Ebend(springpairs,springlist,XY,1,0,1,kappa);
        double Etot0=Est0+Eben0;
        cout<<"E0="<<Etot0<<endl;
        double EstN;
        double EbenN;
        double EtotN;
        
        VectorXd gradE(XY.size());
        VectorXd s0(XY.size());
        gradE=HarmonicGradient(springlist,XY,g11,g12,g22)+gradEbend(springpairs,springlist,XY,g11,g12,g22,kappa);
        s0=-gradE;
        double lengrad0=sqrt(s0.dot(s0));
        double lengradN;
        int count=0;
        do{
            doConjStep(XY,s0,gradE,springlist,springpairs,kappa,count,g11,g12,g22);
            count++;
            lengradN=sqrt(gradE.dot(gradE));
            cout<<count<<"  "<<lengradN/lengrad0<<endl;

        }while(count<Nit && lengrad0<200*lengradN);
        
        
       EstN=Energynetwork(springlist,XY,g11,g12,g22);
       EbenN=Ebend(springpairs,springlist,XY,g11,g12,g22,kappa);
       EtotN=EstN+EbenN;
       cout<<"En="<<EtotN<<endl;
       double dE=EtotN-Etot0;
       double dV=2*eps;
       double P=-dE/dV;
       
       cout<<"THE ENERGY DIFFERENCE "<<dE<<endl;
       return P;       
}

