#include<iostream>
#include<vector>
#include<algorithm>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
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
    for(int i=0;i<2*num;i++){
        bendinggrad(i)=0.0;
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
        
        double dx12=x2-x1;
        double dy12=y2-y1;
        double dy23=y3-y2;
        double dx23=x3-x2;
        double dy13=y3-y1;
        
        
        double d12squared=g11*pow(dx12,2) + 2*g12*dx12*dy12 + g22*pow(dy12,2);
        double d23squared=g11*pow(dx23,2) + 2*g12*dx23*dy23 + g22*pow(dy23,2);
        double d12=sqrt(d12squared);
        double d23=sqrt(d23squared);
    
        bendinggrad(one)=bendinggrad(one)+
        (
            (
                2*(
                    -((-2*g11*dx12 - 2*g12*dy12)*
            (g11*dx12 + g12*dy12))/(2.*pow(d12squared,1.5))   -g11/d12
                )
            *   (
                    (g11*dx12 + g12*dy12)/d12 -(g11*dx23 + g12*dy23)/d23
                ) +
                2*(-((-2*g11*dx12 - 2*g12*dy12)*(g12*dx12 + g22*dy12))/(2.*pow(d12squared,1.5))   
                    -g12/d12)*((g12*dx12 + g22*dy12)/d12 -(g12*dx23 + g22*dy23)/d23)  
            )
            /(d12 + sqrt(g11*pow(dx23,2) + 2*g12*dx23*dy13 + g22*pow(dy13,2)))  
            -   (  (-2*g11*dx12 - 2*g12*dy12)*(  pow((g11*dx12 + g12*dy12)/d12 - (g11*dx23 + g12*dy23)
            /d23,2) +   pow((g12*dx12 + g22*dy12)/d12 -(g12*dx23 + g22*dy23)/d23,2))  )  
            /(2.*d12*pow(d12 + sqrt(g11*pow(dx23,2) + 2*g12*dx23*(dy13) + g22*pow(dy13,2)),2))
        );
        
   
        bendinggrad(two)=  bendinggrad(two)+( (2*  (  -(  (g11*dx12 + g12*dy12)*
        (2*g11*dx12 + 2*g12*dy12) )/(2.*pow(d12squared,1.5))   +g11/d12   +   
        (  (-2*g11*dx23 - 2*g12*dy23)*(g11*dx23 + g12*dy23) )/(2.*pow(d23squared,1.5))   +g11/d23  
        )  *  (  (g11*dx12 + g12*dy12)/d12 - (g11*dx23 + g12*dy23)/d23  )    +   2*(  
        -( (2*g11*dx12 + 2*g12*dy12)*(g12*dx12 + g22*dy12) )/(2.*pow(d12squared,1.5))   +   
        g12/d12 +   (  (-2*g11*dx23 - 2*g12*dy23)*(g12*dx23 + g22*dy23) )/(2.*pow(d23squared,1.5))   
        +   g12/d23  )  *   ( (g12*dx12 + g22*dy12)/d12 - (g12*dx23 + g22*dy23)/d23  )  )/  
        (d12 + sqrt(g11*pow(dx23,2) + 2*g12*dx23*dy13 + g22*pow(dy13,2)))  -    (   (  
        (2*g11*dx12 + 2*g12*dy12)/(2.*d12) +   
        (-2*g11*dx23 - 2*g12*dy13)/(2.*sqrt(g11*pow(dx23,2) + 2*g12*dx23*(dy13) + g22*pow(dy13,2)))   
        )*   (   pow((g11*dx12 + g12*dy12)/d12 -(g11*dx23 + g12*dy23)/d23,2)   
        + pow((g12*dx12 + g22*dy12)/d12 -(g12*dx23 + g22*dy23)/d23,2)   )   )   
        /pow(d12 +sqrt(g11*pow(dx23,2) + 2*g12*dx23*dy13 + g22*pow(dy13,2)),2));
        
        
        
        bendinggrad(three)=bendinggrad(three)+(  (2*(  (  (g11*dx23 + g12*dy23)*
        (2*g11*dx23 + 2*g12*dy23) )/(2.*pow(d23squared,1.5))   -   
        g11/d23)*( (g11*dx12 + g12*dy12)/d12 - (g11*dx23 + g12*dy23)/d23  ) +   
        2*( ( (2*g11*dx23 + 2*g12*dy23)*(g12*dx23 + g22*dy23) )  
        /(2.*pow(d23squared,1.5)) - g12/d23  )*   
        ( (g12*dx12 + g22*dy12)/d12 - (g12*dx23 + g22*dy23)/d23  )  )  /  
        (d12 + sqrt(g11*pow(dx23,2) + 2*g12*dx23*(dy13) + g22*pow(dy13,2))  )   -   
        (  (2*g11*dx23 + 2*g12*(dy13))*  (pow((g11*dx12 + g12*dy12)/d12 - (g11*dx23 + g12*dy23)
        /d23,2) +   pow((g12*dx12 + g22*dy12)/d12 - (g12*dx23 + g22*dy23)/d23,2)  )  )  /  
        (2.*sqrt(g11*pow(dx23,2) + 2*g12*dx23*dy13 + g22*pow(dy13,2))*  
        pow(d12 + sqrt(g11*pow(dx23,2) + 2*g12*dx23*dy13 + g22*pow(dy13,2)),2)  ));

        
        bendinggrad(one+num)=bendinggrad(one+num)+(  (2*  (-( (g11*dx12 + g12*dy12)*
        (-2*g12*dx12 - 2*g22*dy12) )/(2.*pow(d12squared,1.5))   
        - g12/d12)*( (g11*dx12 + g12*dy12)/d12 - (g11*dx23 + g12*dy23)/d23  )   +   
        2*(-(  (-2*g12*dx12 - 2*g22*dy12)*
        (g12*dx12 + g22*dy12) )/(2.*pow(d12squared,1.5)) - g22/d12)  *   
        ( (g12*dx12 + g22*dy12)/d12 - (g12*dx23 + g22*dy23)/d23) )/  
        (d12 + sqrt(g11*pow(dx23,2) + 2*g12*dx23*dy13 + g22*pow(dy13,2)))   -   (   
        ( (-2*g12*dx12 - 2*g22*dy12)/(2.*d12) +   (-2*g12*dx23 - 2*g22*
        (dy13))/(2.*sqrt(g11*pow(dx23,2) + 2*g12*dx23*(dy13) + g22*pow(dy13,2)))  )  *  (   
        pow((g11*dx12 + g12*dy12)/d12 - (g11*dx23 + g12*dy23)/d23,2)    
        + pow((g12*dx12 + g22*dy12)/d12 - (g12*dx23 + g22*dy23)/d23,2)  )    )/  
        pow(d12 +sqrt(g11*pow(dx23,2) + 2*g12*dx23*(dy13) + g22*pow(dy13,2)),2));

        
        bendinggrad(two+num)= bendinggrad(two+num)+( (2*  (  -((g11*dx12 + g12*dy12)*
        (2*g12*dx12 + 2*g22*dy12))/(2.*pow(d12squared,1.5)) +   g12/d12 +   ((g11*dx23 + g12*dy23)*
        (-2*g12*dx23 - 2*g22*dy23))/   (2.*pow(d23squared,1.5)) +   g12/d23)*         
        ((g11*dx12 + g12*dy12)/d12 -   (g11*dx23 + g12*dy23)/d23) +    2*(-((g12*dx12 + g22*dy12)*
        (2*g12*dx12 + 2*g22*dy12))/   (2.*pow(d12squared,1.5)) +   g22/d12 +   
        ((-2*g12*dx23 - 2*g22*dy23)*(g12*dx23 + g22*dy23))/   (2.*pow(d23squared,1.5)) +   
        g22/d23)*   ((g12*dx12 + g22*dy12)/d12 -   (g12*dx23 + g22*dy23)/d23))/  (d12 +    
        sqrt(g11*pow(dx23,2) + 2*g12*dx23*(dy13) + g22*pow(dy13,2))) -    
        ((2*g12*dx12 + 2*g22*dy12)*(pow((g11*dx12 + g12*dy12)/  d12 -    
        (g11*dx23 + g12*dy23)/d23,2) +   pow((g12*dx12 + g22*dy12)/d12 -    
        (g12*dx23 + g22*dy23)/d23,2)))/  (2.*d12*   pow(d12 +   sqrt(g11*pow(dx23,2) + 2*g12*dx23*
        (dy13) + g22*pow(dy13,2)),2)));
        
        bendinggrad(three+num)=bendinggrad(three+num)+( (2*(((g11*dx23 + g12*dy23)*(2*g12*dx23 + 2*g22*dy23))/   
        (2.*pow(d23squared,1.5)) -   g12/d23)*   ((g11*dx12 + g12*dy12)/d12 -   
        (g11*dx23 + g12*dy23)/d23) +    2*(((g12*dx23 + g22*dy23)*(2*g12*dx23 + 2*g22*dy23))/   
        (2.*pow(d23squared,1.5)) -   g22/d23)*   ((g12*dx12 + g22*dy12)/d12 -   
        (g12*dx23 + g22*dy23)/d23))/  (d12 +    sqrt(g11*pow(dx23,2) + 2*g12*dx23*
        dy13 + g22*pow(dy13,2))) -    ((2*g12*dx23 + 2*g22*(dy13))*(pow((g11*dx12 + g12*dy12)/  d12 -    
        (g11*dx23 + g12*dy23)/d23,2) +   pow((g12*dx12 + g22*dy12)/d12 -    
        (g12*dx23 + g22*dy23)/d23,2)))/  (2.*sqrt(g11*pow(dx23,2) + 2*g12*dx23*(dy13) + g22*pow(dy13,2))*   
        pow(d12 +   sqrt(g11*pow(dx23,2) + 2*g12*dx23*(dy13) + g22*pow(dy13,2)),2)));



    }

    return kappa*bendinggrad;
}

double l_ij(int springnr, VectorXd &XY,vector<spring> &springlist,double g11, double g12,double g22){
    double l=0.0;
    int NUM=XY.size()/2;
    int one=springlist[springnr].one;
    int two=springlist[springnr].two;
    int wlr=springlist[springnr].wlr;
    int wud=springlist[springnr].wud;
    
    double x1,x2,y1,y2;
    x1=XY(one);
    x2=XY(two)+1.0*wlr;
    y1=XY(one+NUM);
    y2=XY(two+NUM)+1.0*wud;
    
    l=sqrt(g11*pow((x2-x1),2)+2*g12*(x2-x1)*(y2-y1)+g22*pow((y2-y2),2));
    return l;
}




double F_for_bending(int springone, int springtwo,vector<spring> &springlist,VectorXd &XY,double g11,double g12,double g22){
    double F=0.0;
    int NUM=XY.size();
    int one,two,three;
    one=springlist[springone].one;
    two=springlist[springone].two;
    three=springlist[springtwo].two;
    
    double x1,x2,x3,y1,y2,y3;
    x1=XY(one);
    x2=XY(two)+springlist[springone].wlr*1.0;
    x3=XY(three)+1.0*springlist[springone].wlr+1.0*springlist[springtwo].wlr;
    
    y1=XY(one+NUM);
    y2=XY(two+NUM)+1.0*springlist[springone].wud;
    y3=XY(three+NUM)+1.0*springlist[springone].wud+1.0*springlist[springtwo].wud;
    
    double dx12=x2-x1;
    double dx23=x3-x2;
    double dy12=y2-y1;
    double dy23=y3-y2;
    double l12=l_ij(springone,XY,springlist,g11,g12,g22);
    double l23=l_ij(springtwo,XY,springlist,g11,g12,g22);
    
    F=pow(((g11*dx12+g12*dy12)/l12)-((g11*dx23+g12*dy23)/l23),2)+pow(((g22*dy12+g12*dx12)/l12)-((g22*dy23+g12*dx23)/l23),2);
    
    return F;
}









// VectorXd BendingGrad1(const vector<vector<int>> &springpairs,
//                      const vector<spring> &springlist,
//                      const VectorXd &XY,
//                      double kappa,
//                      double g11,
//                      double g12,
//                      double g22)
// {
//     
//     int num=XY.size()/2;
//     VectorXd bendinggrad(2*num);
//     
//     for(int i=0;i<2*num;i++){
//         bendinggrad(i)=0.0;
//     }
//     
//     
//     for(int i=0; i<springpairs.size(); i++){
//         int springone=springpairs[i][0];
//         int springtwo=springpairs[i][1];
//         int one=springlist[springone].one;
//         int two=springlist[springone].two;
//         int three=springlist[springtwo].two;
//         
//         double x1=XY(one);
//         double y1=XY(one+num);
//         double x2=XY(two)+springlist[springone].wlr;
//         double y2=XY(two+num)+springlist[springone].wud;
//         double x3=XY(three)+springlist[springone].wlr+springlist[springtwo].wlr;
//         double y3=XY(three+num)+springlist[springone].wud+springlist[springtwo].wud;
//         
//         double dx12=x2-x1;
//         double dy12=y2-y1;
//         double dy23=y3-y2;
//         double dx23=x3-x2;
//         double dy13=y3-y1;
//         
//         
//         double d12squared=g11*pow(dx12,2) + 2*g12*dx12*dy12 + g22*pow(dy12,2);
//         double d23squared=g11*pow(dx23,2) + 2*g12*dx23*dy23 + g22*pow(dy23,2);
//         double d12=sqrt(d12squared);
//         double d23=sqrt(d23squared);
    
        
        //The bending gradient is grad(1/(l12+l23)* F^2 + 1/(l12+l23)*grad(F^2)
        //Let us make a function F so we can call it later. 
        
        
        
        
        
    //}
    
    
    
    
    
    
    
    
    
    
//}






