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
#include "writefunctions.h"
#include "shearleader.h"

using namespace std;
using namespace Eigen;

VectorXd makeSquareNetwork(int Number, vector<spring> &springlist){
      ofstream squarelatice("sqlattice.txt");
      double dx=1.0/(Number);
      double dy=dx;
      double x=0.5*dx;
      double y=x;
      VectorXd X(Number*Number);
      VectorXd Y(Number*Number);
      VectorXd XY(2*X.size());
      
      int nodenumber=0;
      //make the position vector
      for(int j=0;j<Number;j++){
	y=0.5*dy+j*dy;
	for(int i=0;i<Number;i++){
	   Y(nodenumber)=y;
	   X(nodenumber)=x+i*dx;
	   nodenumber++;
	}
    }
    for(int i =0; i<X.size();i++){
     //cout<<X(i)<<"\t"<<Y(i)<<endl;
     squarelatice<<X(i)<<"\t"<<Y(i)<<endl;
    }
   for(int i=0;i<X.size();i++){
       XY(i)=X(i);
       XY(i+X.size())=Y(i);
   }

   spring newspring;
    
    int sticki=0;
    for(int j=0;j<Number;j++){
      
      for(int i=0;i<Number-1;i++){
	//These are the horizontal springs
	newspring.one=(i)+(j*Number);
	newspring.two=(i+1)+(j*Number);
	newspring.sticki=-1;
	newspring.rlen=dx;
	newspring.k=1.0;
	newspring.wlr=0;
	newspring.wud=0;
	springlist.push_back(newspring);
      }
      newspring.one=(Number-1)+(j*Number);
      newspring.two=j*(Number);
      newspring.sticki=-1;
      newspring.rlen=dx;
      newspring.k=1.0;
      newspring.wlr=1;
      newspring.wud=0;
      springlist.push_back(newspring);
  
      sticki++;
    }

  for(int i=0;i<Number;i++){
    //these are the vertical springs
   for(int j=0;j<Number-1;j++){
    newspring.one=i+j*Number;
    newspring.two=i+(j+1)*Number;
    newspring.sticki=-1;
    newspring.rlen=dy;
    newspring.k=1.0;
    newspring.wlr=0;
    newspring.wud=0;
    springlist.push_back(newspring);
     
   }
    newspring.one=i+(Number-1)*Number;
    newspring.two=i;
    newspring.sticki=-1;
    newspring.rlen=dy;
    newspring.k=1.0;
    newspring.wlr=0;
    newspring.wud=1;
    springlist.push_back(newspring);
    sticki++; 
  }

    return XY;
}


VectorXd triangularnetwork(int Number, vector<spring> &springlist){
    double dx=1.0/Number;
    double dy=dx;
    double x1=1.0/(3*Number); //two types of rows
    double x2=5.0/(6.0*Number);
    double y1=x1;
    double x,y;
    VectorXd X(Number*Number);
    VectorXd Y(Number*Number); 
    VectorXd XY(2*X.size());
    ofstream TRINET("triangular.txt");
    int particle=0;
    int rownr=0;
    for(int i=0; i<Number;i++){
        rownr=i;
        y=y1+rownr*dy;
        for(int j=0; j<Number;j++){
            //check if we are in an even or an uneven row
            if(rownr%2==0){
                x=x1+j*dx;
            }else{
                x=x2+j*dx;
            }
            X(particle)=x;
            Y(particle)=y;
            particle++;
        }
    }

    for(int i=0;i<X.size();i++){
        XY(i)=X(i);
        XY(i+X.size())=Y(i);
    }

    for(int i=0;i<X.size();i++){
       TRINET<<X(i)<<"\t"<<Y(i)<<endl; 
    }

    //Now assign the springs
    int currentnode;
    spring newspring;
    double dist;
    for(int i=0;i<Number;i++){ //row;
       for(int j=0;j<Number;j++){
           currentnode=Number*i+j;
           if(i>0 && i<Number-1){//not the top and bottom row 
               if(j>0 && j<Number-1 && i%2!=0){ //and not the l/r row EVEN ROWS
                   springlist.push_back(makespring(currentnode,currentnode+1,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode-1,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode+Number,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode+Number+1,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode-Number,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode-Number+1,X,Y,0,0));
               }
               
                if(j>0 && j<Number-1 && i%2==0){ //and not the l/r row UNEVEN
                   springlist.push_back(makespring(currentnode,currentnode+1,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode-1,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode+Number,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode+Number-1,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode-Number,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode-Number-1,X,Y,0,0));
               }
               //the right boundary for the even rows
               if(j==Number-1 && i%2==0){
                   springlist.push_back(makespring(currentnode,currentnode-1,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode+Number-1,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode+Number,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode-Number-1,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode-Number,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode-Number+1,X,Y,1,0));
               }
               //Even Rows left boundary
               if(j==0 && i%2==0){
                   springlist.push_back(makespring(currentnode,currentnode+1,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode+Number-1,X,Y,-1,0));
                   springlist.push_back(makespring(currentnode,currentnode+Number,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode+2*Number-1,X,Y,-1,0));
                   springlist.push_back(makespring(currentnode,currentnode-1,X,Y,-1,0));
                   springlist.push_back(makespring(currentnode,currentnode-Number,X,Y,0,0));
               }
               //Uneven rows left boundary
               if(j==0 && i%2!=0){
                   springlist.push_back(makespring(currentnode,currentnode+1,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode+Number,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode+Number+1,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode+Number-1,X,Y,-1,0));
                   springlist.push_back(makespring(currentnode,currentnode-Number,X,Y,0,0));
                   springlist.push_back(makespring(currentnode,currentnode-Number+1,X,Y,0,0));
               }
           }
           //Top row
           if(i==Number-1){
                if(j!=0 && j!=Number-1){
                    springlist.push_back(makespring(currentnode,currentnode+1,X,Y,0,0));
                    springlist.push_back(makespring(currentnode,currentnode-1,X,Y,0,0));
                    springlist.push_back(makespring(currentnode,currentnode-Number,X,Y,0,0));
                    springlist.push_back(makespring(currentnode,currentnode-Number+1,X,Y,0,0));
                    springlist.push_back(makespring(currentnode,currentnode-(Number-1)*Number,X,Y,0,1));
                    springlist.push_back(makespring(currentnode,currentnode-(Number-1)*Number+1,X,Y,0,1));
                }
        }
           //bottom row
            if(i==0){
                if(j!=0 &&j!=Number-1){
                        springlist.push_back(makespring(currentnode,currentnode+1,X,Y,0,0));
                        springlist.push_back(makespring(currentnode,currentnode-1,X,Y,0,0));
                        springlist.push_back(makespring(currentnode,currentnode+Number-1,X,Y,0,0));
                        springlist.push_back(makespring(currentnode,currentnode+Number,X,Y,0,0));
                        springlist.push_back(makespring(currentnode,currentnode+(Number-1)*Number-1,X,Y,0,-1));
                        springlist.push_back(makespring(currentnode,currentnode+(Number-1)*Number,X,Y,0,-1));
                }
            }

        //finally the edges : i=0,j=0 ; i=0,j=N-1 ; i=N-1,j=0; i=N-1; j=N-1
            if(i==0 && j==0){
                        springlist.push_back(makespring(currentnode,currentnode+1,X,Y,0,0));
                        springlist.push_back(makespring(currentnode,currentnode+Number-1,X,Y,-1,0));
                        springlist.push_back(makespring(currentnode,currentnode+Number,X,Y,0,0));
                        springlist.push_back(makespring(currentnode,currentnode+2*Number-1,X,Y,-1,0));
                        springlist.push_back(makespring(currentnode,currentnode+(Number-1)*Number,X,Y,0,-1));
                        springlist.push_back(makespring(currentnode,currentnode+(Number)*Number-1,X,Y,-1,-1));
                
            }
            
            if(i==0 && j==Number-1){
                springlist.push_back(makespring(currentnode,currentnode-1,X,Y,0,0));
                springlist.push_back(makespring(currentnode,currentnode-Number+1,X,Y,-1,0));
                springlist.push_back(makespring(currentnode,currentnode+Number-1,X,Y,0,0));
                springlist.push_back(makespring(currentnode,currentnode+Number,X,Y,0,0));
                springlist.push_back(makespring(currentnode,currentnode+(Number-1)*Number,X,Y,0,-1));
                springlist.push_back(makespring(currentnode,currentnode+(Number-1)*Number-1,X,Y,0,-1));
                
            }
            
            if(i==Number-1 && j==0){
                springlist.push_back(makespring(currentnode,currentnode+1,X,Y,0,0));
                springlist.push_back(makespring(currentnode,currentnode+Number-1,X,Y,-1,0));
                springlist.push_back(makespring(currentnode,currentnode-Number,X,Y,0,0));
                springlist.push_back(makespring(currentnode,currentnode-Number+1,X,Y,0,0));
                springlist.push_back(makespring(currentnode,currentnode-(Number-1)*Number,X,Y,0,1));
                springlist.push_back(makespring(currentnode,currentnode-(Number-1)*Number+1,X,Y,0,1));
                
            }
            
            if(i==Number-1 && j==Number-1){
                springlist.push_back(makespring(currentnode,currentnode-1,X,Y,0,0));
                springlist.push_back(makespring(currentnode,currentnode-Number+1,X,Y,1,0));
                springlist.push_back(makespring(currentnode,currentnode-Number,X,Y,0,0));
                springlist.push_back(makespring(currentnode,currentnode-2*Number+1,X,Y,1,0));
                springlist.push_back(makespring(currentnode,currentnode-(Number-1)*Number,X,Y,0,1));
                springlist.push_back(makespring(currentnode,currentnode-(Number)*Number+1,X,Y,1,1));
                
            }
       }
    }
    
    //finally remove double springs
    int gotit;
    vector<spring> uniquesprings;
    uniquesprings.push_back(springlist[0]);
    for(int i=0;i<springlist.size();i++){
        for(int j=0;j<uniquesprings.size();j++){
            if((springlist[i].one==uniquesprings[j].one && springlist[i].two==uniquesprings[j].two)
                ||(springlist[i].one==uniquesprings[j].two && springlist[i].two==uniquesprings[j].one))
            {
             gotit=1;
             break;
            }
        }
        if(gotit==0){
            uniquesprings.push_back(springlist[i]);
        }
        gotit=0;
    }
    springlist=uniquesprings;
    
    
    
    
    
    return XY; 
}




spring makespring(int one, int two,VectorXd &X, VectorXd &Y,int wlr,int wud){
            double dist;
            spring newspring;
            newspring.one=one;
            newspring.two=two;
            newspring.k=1.0;
            dist=sqrt(pow((X(newspring.one)-(X(newspring.two)+wlr)),2)+pow((Y(newspring.one)-(Y(newspring.two)+wud)),2));
            newspring.rlen=dist;
            newspring.wlr=wlr;
            newspring.wud=wud;
            newspring.sticki=-1;
            
            return newspring;
}














