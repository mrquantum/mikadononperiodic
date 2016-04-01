#include<iostream>
#include<vector>
#include<algorithm>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "makemikadonetwork.h"
#include "simpleBendingGrad.h"

const double pi=4.0*atan(1.0);

using namespace std;
using namespace Eigen;
double Ebending(const vector<vector<int>> &springpairs,
                const vector<spring> &springlist,
                const VectorXd &XY,
                const double kappa,
                double g11,
                double g12,
                double g22)
{
    double E=0.0;
    for(int i=0;i<springpairs.size();i++){
        int springone=springpairs[i][0];
        int springtwo=springpairs[i][1];
        int num=XY.size()/2;
        
        int one=springlist[springone].one;
        int two=springlist[springone].two;
        int three=springlist[springtwo].two;
        
        double x1=XY(one);
        double y1=XY(one+num);
        double x2=XY(two)+springlist[springone].wlr;
        double y2=XY(two+num)+springlist[springone].wud;
        double x3=XY(three)+springlist[springone].wlr+springlist[springtwo].wlr;
        double y3=XY(three+num)+springlist[springone].wud+springlist[springtwo].wud;
        double l1=sqrt(g11*(x2-x1)*(x2-x1)+2*g12*(x2-x1)*(y2-y1)+g22*(y2-y1)*(y2-y1));
        
        double l2=sqrt(g11*(x3-x2)*(x3-x2)+2*g12*(x3-x2)*(y3-y2)+g22*(y3-y2)*(y3-y2));
        double v1dotv2=g11*(x2-x1)*(x3-x2)+
                       g12*(x2-x1)*(y3-y2)+
                       g12*(y2-y1)*(x3-x2)+
                       g22*(y2-y1)*(y3-y2);
        double costh=v1dotv2/(l1*l2);
        if(costh>1.0) costh=1.0;
        if(costh<-1.0) costh=-1.0;
        E=E+kappa*pow(acos(costh),2)/(l1+l2);
        //cout<<l1<<"\t"<<l2<<"\t"<<costh<<"\t"<<E<<endl;

    }
    return E;
}




double EbendingC(const vector< vector< int > > &springpairs, 
                 const vector< spring > &springlist, 
                 const VectorXd &XY, 
                 const double kappa, 
                 const double g11, const double g12, const double g22)
{

    double E=0.0;
    for(int i=0;i<springpairs.size();i++){
        int springone=springpairs[i][0];
        int springtwo=springpairs[i][1];
        int num=XY.size()/2;
        
        int one=springlist[springone].one;
        int two=springlist[springone].two;
        int three=springlist[springtwo].two;
        
        double x1=XY(one);
        double y1=XY(one+num);
        double x2=XY(two)+springlist[springone].wlr;
        double y2=XY(two+num)+springlist[springone].wud;
        double x3=XY(three)+springlist[springone].wlr+springlist[springtwo].wlr;
        double y3=XY(three+num)+springlist[springone].wud+springlist[springtwo].wud;
        
        double l1=sqrt(g11*(x2-x1)*(x2-x1)+2*g12*(x2-x1)*(y2-y1)+g22*(y2-y1)*(y2-y1));
        double l2=sqrt(g11*(x3-x2)*(x3-x2)+2*g12*(x3-x2)*(y3-y2)+g22*(y3-y2)*(y3-y2));
        
        double thetasq=pow(((g11*(x2-x1)+g12*(y2-y1))/l1)-((g11*(x3-x2)+g12*(y3-y2))/l2),2)+
                       pow(((g22*(y2-y1)+g12*(x2-x1))/l1)-((g22*(y3-y2)+g12*(x3-x2))/l2),2);
    
        E=E+kappa*thetasq/(l1+l2);
    }
    
    return E;
    
    
    
}

double EbendingCn(double *XY, networkinfo parameters) //This is the wrapper E funct that the CG frmp accepts
{
    const vector<spring> springlist=parameters.springlist;
    const vector<vector<int>> springpairs=parameters.springpairs;
    double g11=parameters.g11;
    double g12=parameters.g12;
    double g22=parameters.g22;
    double kappa=parameters.kappa;
    int size=parameters.size;
    int num=size/2;
    double sheardeformation=parameters.sheardeformation;
    double E=0.0;
    
    Eigen::Map<Eigen::VectorXd> XY_temp( XY, size );
    E=BendEnergy(springpairs,springlist,XY_temp,kappa,sheardeformation);
    
    return E;
}