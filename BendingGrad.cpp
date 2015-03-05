#include<iostream>
#include<vector>
#include<algorithm>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include "makemikadonetwork.h"

const double pi=4.0*atan(1.0);

using namespace std;
using namespace Eigen;

VectorXd BendingGrad(const vector<vector<int>> &springpairs,
                     const vector<spring> &springlist,
                     const VectorXd &XY,
                     double kappa,
                     double g11,
                     double g12,
                     double g22)
{
    int num=XY.size()/2;
    VectorXd bendinggrad(2*num);
    VectorXd gradp1(2*num);
    VectorXd gradp2(2*num);
    VectorXd nom(2*num);
    VectorXd den(2*num);
    for(int i=0;i<2*num;i++){
        bendinggrad(i)=0.0;
        gradp1(i)=0.0;
        gradp2(i)=0.0;
        nom(i)=0.0;
        den(i)=0.0;
    }
    
    for(int i=0; i<springpairs.size(); i++){
        int springone=springpairs[i][0];
        int springtwo=springpairs[i][1];
        int one=springlist[springone].one;
        int two=springlist[springone].two;
        int three=springlist[springtwo].two;
        
        double x1=XY(one);
        double y1=XY(one+num);
        double x2=XY(two)+springlist[springone].wlr;
        double y2=XY(two+num)+springlist[springone].wud;
        double x3=XY(three)+springlist[springone].wlr+springlist[springtwo].wlr;
        double y3=XY(three+num)+springlist[springone].wud+springlist[springtwo].wud;
        
        double l1=sqrt(g11*(x2-x1)*(x2-x1)+2*g12*(x2-x1)*(y2-y1)+g22*(y2-y1)*(y2-y1)); //len of v1
        double l2=sqrt(g11*(x3-x2)*(x3-x2)+2*g12*(x3-x2)*(y3-y2)+g22*(y3-y2)*(y3-y2)); //len of v2
        double v1dotv2=g11*(x2-x1)*(x3-x2)+
                       g12*(x2-x1)*(y3-y2)+
                       g12*(y2-y1)*(x3-x2)+
                       g22*(y2-y1)*(y3-y2);
        double costh=v1dotv2/(l1*l2);
        if(costh<-1.0) costh=-1.0;
        if(costh>1.0) costh=1.0;
        double sinth=sqrt(1-costh*costh);
        if(sinth<0.0001) sinth=0.0001;
        double invsin=1.0/sinth;        
        //cout<<invsin<<endl;
        double th=acos(costh);
        //cout<<th<<"\t"<<invsin<<endl;
        //firt the part w. grad(1/l1+l2))
        double dgradL1x=(g11*(x2-x1)+g12*(y2-y1))/l1;
        double dgradL1y=(g22*(y2-y1)+g12*(x2-x1))/l1;
        double dgradL2x=(g11*(x3-x2)+g12*(y3-y2))/l2;
        double dgradL2y=(g22*(y3-y2)+g12*(x3-x2))/l2;
        //Update gradient components
        gradp1(one)-=-pow((th/(l1+l2)),2)*dgradL1x;
        gradp1(two)+=-pow((th/(l1+l2)),2)*dgradL1x;
        gradp1(one+num)-=-pow((th/(l1+l2)),2)*dgradL1y;
        gradp1(two+num)+=-pow((th/(l1+l2)),2)*dgradL1y;
        
        gradp1(two)-=-pow((th/(l1+l2)),2)*dgradL2x;
        gradp1(three)+=-pow((th/(l1+l2)),2)*dgradL2x;
        gradp1(two+num)-=-pow((th/(l1+l2)),2)*dgradL2y;
        gradp1(three+num)+=-pow((th/(l1+l2)),2)*dgradL2y;
    
        //Now we must add 1/l1+l2 * grad th^2
        double nominator=v1dotv2;
        double denuminator=l1*l2;
        
        nom(one)=-g11*(x3-x2)-g12*(y3-y2);
        nom(two)=g11*(x3-x2)+g12*(y3-y2)-g11*(x2-x1)-g12*(y2-y1);
        nom(three)=g11*(x2-x1)+g12*(y2-y1);
        nom(one+num)=-g22*(y3-y2)-g12*(x3-x2);
        nom(two+num)=g22*(y3-y2)+g12*(x3-x2)-g22*(y2-y1)-g12*(x2-x1);
        nom(three+num)=g22*(y2-y1)+g12*(x2-x1);
        
        den(one)=-(l2)*dgradL1x;
        den(two)=(l2)*dgradL1x;
        den(one+num)=-(l2)*dgradL1y;
        den(two+num)=(l2)*dgradL1y;
        den(two)=-(l1)*dgradL2x;
        den(three)=(l1)*dgradL2x;
        den(two+num)=-(l1)*dgradL2y;
        den(three+num)=(l1)*dgradL2y;

        gradp2+=-2*(1.0/(l1+l2))*th*invsin*(nominator*den-denuminator*nom)/(denuminator*denuminator);
        
        //Reset the den and num components;
        nom(one)=0.0;
        nom(two)=0.0;
        nom(three)=0.0;
        nom(one+num)=0.0;
        nom(two+num)=0.0;
        nom(three+num)=0.0;
        den(one)=0.0;
        den(two)=0.0;
        den(three)=0.0;
        den(one+num)=0.0;
        den(two+num)=0.0;
        den(three+num)=0.0;

    }
    bendinggrad=kappa*(gradp1+gradp2);
    //bendinggrad=kappa*gradp2;
    return bendinggrad;
}