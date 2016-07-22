#include<iostream>
#include<vector>
#include<algorithm>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "makemikadonetwork.h"
#include "simpleBendingGrad.h"
#include "newbendinggrad.h"

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


void BendinggradN(double *p,double *xi,networkinfo params)//wrapper function for the bending gradient
{
    double g11=params.g11;
    double g12=params.g12;
    double g22=params.g22;
    vector<spring> springlist=params.springlist;
    vector<vector<int>> springpairs=params.springpairs;
    int size=params.size;
    int one,two,num=size/2;
    double kappa=params.kappa;
    double sheardeformation=params.sheardeformation;
//     VectorXd gradB(size);
    triplet springpairsandkappa;
    //VectorXd effkappa=params.effkappa;
    //springpairsandkappa.EffKappa=effkappa;
    springpairsandkappa.springpairs=springpairs;
    
    
    Eigen::Map<Eigen::VectorXd> XY_temp( p, size );
    Eigen::VectorXd XY_ref = XY_temp;
    bendinggradnew(p,xi,size,springlist,springpairs,kappa,sheardeformation);

}













double L_ij(double x1,double y1,double x2,double y2){
    return sqrt(pow((x2-x1),2)+pow((y2-y1),2));
}


vector<gradstruct> gradL(int one, int two,double x1,double y1,double x2,double y2,double Lij)
//takes two interacting nodes connected by a spring in physical coordinates.
//also takes the positions and the networkinfo. 
{
    vector<gradstruct> gradb(0);
    
    gradstruct e1,e2;
    e1.ddxi=(x1-x2)/Lij;
    e1.ddyi=(y1-y2)/Lij;
    e1.i=one;
    
    e2.ddxi=-e1.ddxi;
    e2.ddyi=-e1.ddyi;
    e2.i=two;
    
    gradb.push_back(e1);
    gradb.push_back(e2);
    
    return gradb;
}

double thsq(double x1,double y1, double x2, double y2, double x3, double y3, double L12, double L23)
{
    double f=pow((  (x2-x1)/L12  - (x3-x2)/L23 ),2)+pow((  (y2-y1)/L12  - (y3-y2)/L23 ),2);
    return f;
    
}


//the bendinggradient has two components: 
//1) grad(1/l12+l23)*thsq, thsq is approximated by f
//2) grad(thsq)*1/(l12+l23)

vector<gradstruct> bendgradpart1(double f, double L12, double L23,vector<gradstruct> &gradL12,vector<gradstruct> &gradL23)
{
    vector<gradstruct> grad(0);
    gradstruct element;
    double prefactor=-f/pow((L12+L23),2);
    for(int j=0;j<gradL12.size();j++){
        element.ddxi=prefactor*gradL12[j].ddxi;
        element.ddyi=prefactor*gradL12[j].ddyi;
        element.i=gradL12[j].i;
        grad.push_back(element);
    }
    
    for(int j=0;j<gradL23.size();j++){
        element.ddxi=prefactor*gradL23[j].ddxi;
        element.ddyi=prefactor*gradL23[j].ddyi;
        element.i=gradL23[j].i;
        grad.push_back(element);
    }
    return grad;
}

vector<gradstruct> bendgradpart2(double L12, double L23, 
                             double x1,double y1, double x2, double y2, double x3, double y3,
                             int i1,int i2,int i3,
                             vector<gradstruct> &gradL12, vector<gradstruct> &gradL23)
{
    gradstruct el1,el2,el3,el4,element;
    vector<gradstruct> grad(0);
    double prefx=2*((x2-x1)/L12-(x3-x2)/L23);
    double prefy=2*((y2-y1)/L12-(y3-y2)/L23);
    double prefactor1=1.0/(L12+L23);
    double prefactor2=1.0/pow(L12,2);
    double prefactor3=1.0/pow(L23,2);
    
    
    el1.ddxi=prefactor1*prefx/L12;
    el1.ddyi=prefactor1*prefy/L12;
    el1.i=i2;
    
    el2.ddxi=-prefactor1*prefx/L12;
    el2.ddyi=-prefactor1*prefy/L12;
    el2.i=i1;
    
    el3.ddxi=-prefactor1*prefx/L23;
    el3.ddyi=-prefactor1*prefy/L23;
    el3.i=i3;
    
    el4.ddxi=prefactor1*prefx/L23;
    el4.ddyi=prefactor1*prefy/L23;
    el4.i=i2;
    
    grad.push_back(el1);
    grad.push_back(el2);
    grad.push_back(el3);
    grad.push_back(el4);

    double pg1=-prefactor1*prefx*(x2-x1)*prefactor2;
    for(int k=0;k<gradL12.size();k++){
        element.ddxi=gradL12[k].ddxi*pg1;
        element.ddyi=gradL12[k].ddyi*pg1;
        element.i=gradL12[k].i;
        grad.push_back(element);
    }
    double pg2=prefactor1*prefx*(x3-x2)*prefactor3;
    for(int k=0;k<gradL23.size();k++){
        element.ddxi=gradL23[k].ddxi*pg2;
        element.ddyi=gradL23[k].ddyi*pg2;
        element.i=gradL23[k].i;
        grad.push_back(element);
    }
    
    double pg3=-prefactor1*prefy*(y2-y1)*prefactor2;
    for(int k=0;k<gradL12.size();k++){
        element.ddxi=gradL12[k].ddxi*pg3;
        element.ddyi=gradL12[k].ddyi*pg3;
        element.i=gradL12[k].i;
        grad.push_back(element);
    }
    
    double pg4=-prefactor1*prefy*(y3-y2)*prefactor3;
    for(int k=0;k<gradL23.size();k++){
        element.ddxi=gradL23[k].ddxi*pg4;
        element.ddyi=gradL23[k].ddyi*pg4;
        element.i=gradL23[k].i;
        grad.push_back(element);
    }
    return grad;
    
}

void physbendinggradient(double *p, double *g,networkinfo &info)
//input: positions and gradient and the network info struct
{
    //loop over the springpairs
    vector<vector<int>> springpairs=info.springpairs;
    vector<spring> springlist=info.springlist;
    //loop over springpairs later
    int springone, springtwo;
    int one,two,three;
    double x1,y1,x2,y2,x3,y3,L12,L23;
    double thetasq;
    int num=(info.size)/2;
    double kappa=info.kappa;
    int wlr1,wlr2,wud1,wud2;
    vector<gradstruct> gradLen12(0);
    vector<gradstruct> gradLen23(0);
    vector<gradstruct> bendgrad1(0);
    vector<gradstruct> bendgrad2(0);
    
    //this is gonna be a loop. 
    for(int ii=0;ii<springpairs.size();ii++){
        springone=springpairs[ii][0];
        springtwo=springpairs[ii][1];
        
        one=springlist[springone].one;
        two=springlist[springone].two;
        three=springlist[springtwo].two;
        wlr1=springlist[springone].wlr;
        wud1=springlist[springone].wud;
        wlr2=springlist[springtwo].wlr;
        wud2=springlist[springtwo].wud;
        
        x1=p[one];
        y1=p[one+num];
        x2=p[two]+wlr1;
        y2=p[two+num]+wud1;
        x3=p[three]+wlr1+wlr2;
        y3=p[three+num]+wud1+wud2;
        
        L12=L_ij(x1,y1,x2,y2);
        L23=L_ij(x2,y2,x3,y3);
        
        thetasq=thsq(x1,y1,x2,y2,x3,y3,L12,L23);
        gradLen12=gradL(one,two,x1,y1,x2,y2,L12);
        gradLen23=gradL(two,three,x2,y2,x3,y3,L23);
        
        bendgrad1=bendgradpart1(thetasq,L12,L23,gradLen12,gradLen23);
        bendgrad2=bendgradpart2(L12,L23,x1,y1,x2,y2,x3,y3,one,two,three,gradLen12,gradLen23);
    
        //now fill the gradient vector
        int pos;
        for(int k=0;k<bendgrad1.size();k++){
            pos=bendgrad1[k].i;
            g[pos]+=kappa*bendgrad1[k].ddxi;
            g[pos+num]+=kappa*bendgrad1[k].ddyi;
        }
        for(int k=0;k<bendgrad2.size();k++){
            pos=bendgrad2[k].i;
            g[pos]+=kappa*bendgrad2[k].ddxi;
            g[pos+num]+=kappa*bendgrad2[k].ddyi;
        }
        
        cout<<"comps first part"<<endl;
         for(int m=0;m<bendgrad1.size();m++){
            cout<<bendgrad1[m].i<<"  "<<bendgrad1[m].ddxi<<"  "<<bendgrad1[m].ddyi<<endl;
        }
        
        cout<<"comps scnd part"<<endl;

        for(int m=0;m<bendgrad2.size();m++){
            cout<<bendgrad2[m].i<<"  "<<bendgrad2[m].ddxi<<"  "<<bendgrad2[m].ddyi<<endl;
        }
    }
    
    
}


