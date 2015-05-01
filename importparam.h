#ifndef IMPORTPARAM_H
#define IMPORTPARAM_H

#include <iostream>
#include <ctime>
#include "random.h"
#include "makemikadonetwork.h"
#include "EnergyandGradients.h"
#include "minimizers.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <math.h>
#include <functional>

using namespace std;

class param
{
public: double rlenshort,rlenlong,kappa,k1,k2,LStick,tolGradE,stretchf,StepSize;
        int NumberMikado,Nit,NumberStepsRight,NumberStepsLeft; //Number of mikado, max number of conjugate iterations
    param()
    {
     rlenshort=.0;
     rlenlong=.0;
     kappa=.0;
     stretchf=1.0;
     k1=.0;
     k2=.0;
     LStick=.0;
     NumberMikado=0;
     Nit=0;
     tolGradE=.0;        
    }  
    
    ~param(){         
    }
    

    int validate() //Checks the parameters 
    {
            if(rlenshort<=0.0) {
                cout<<"restlength should be >0"<<endl;
                return 0;
            }
            if(rlenlong<0.0) {
                cout<<"restlength should be >0"<<endl;
                return 0;
            }
             if(kappa<0.0) {
                cout<<"bending rigidity should be >0"<<endl;
                return 0;
            }
                 if(k1<0.0) {
                cout<<"spring constant >0"<<endl;
                return 0;
            }
                if(k2<0.0) {
                cout<<"spring constant >0"<<endl;
                return 0;
            }
                if(LStick<0.0) {
                cout<<"Length stick is <0"<<endl;
                return 0;
            }
               if(NumberMikado<2) {
                cout<<"Make a better network"<<endl;
                return 0;
            }
                   if(Nit<=0) {
                cout<<"No conjugate steps are made"<<endl;
                return 0;
            } else{
                return 1; 
            }
        
    }
    
    int parse_param_assignment(char *str)
    {
        const int buffsize=101;
        char lvalue[buffsize], rvalue[buffsize];
        int rv=sscanf(str,"%100[^=]=%s",lvalue,rvalue);
        if(rv!=2){
            cout<<"error in param assignment"<<endl;
        }
        //rv=sscanf(lvalue,"%s",lvalue);
        cout<<lvalue<<endl;
        
        if(strcmp(lvalue,"rlenshort")==0){
            rlenshort=atof(rvalue);
        }
        else if(strcmp(lvalue,"rlenlong")==0){
            rlenlong=atof(rvalue);
        }
        else if(strcmp(lvalue,"kappa")==0){
            kappa=atof(rvalue);
        }
        else if(strcmp(lvalue,"k1")==0){
            k1=atof(rvalue);
        }
        else if(strcmp(lvalue,"k2")==0){
            k2=atof(rvalue);
        }
        else if(strcmp(lvalue,"Lstick")==0){
            LStick=atof(rvalue);
        }
        else if(strcmp(lvalue,"NumberMikado")==0){
            NumberMikado=atoi(rvalue);
        }
        else if(strcmp(lvalue,"Nit")==0){
            Nit=atoi(rvalue);
        }
        else if(strcmp(lvalue,"tolGradE")==0){
            tolGradE=atof(rvalue);
        }
        else if(strcmp(lvalue,"stretchf")==0){
            stretchf=atof(rvalue);
        }
        else if(strcmp(lvalue,"StepSize")==0){
            StepSize=atof(rvalue);
        }
        else if(strcmp(lvalue,"NumberStepsLeft")==0){
            NumberStepsLeft=atoi(rvalue);
        }
         else if(strcmp(lvalue,"NumberStepsRight")==0){
            NumberStepsRight=atoi(rvalue);
         }
        else{
            cout<<"Did not recognise param"<<endl;
            return -1;            
            
        }
        return 1;   
        
    }
    
    
};













#endif // IMPORTPARAM_H