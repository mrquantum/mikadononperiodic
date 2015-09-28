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
	newspring.sticki=sticki;
	newspring.rlen=dx;
	newspring.k=1.0;
	newspring.wlr=0;
	newspring.wud=0;
	springlist.push_back(newspring);
      }
      newspring.one=(Number-1)+(j*Number);
      newspring.two=j*(Number);
      newspring.sticki=sticki;
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
    newspring.sticki=sticki;
    newspring.rlen=dy;
    newspring.k=1.0;
    newspring.wlr=0;
    newspring.wud=0;
    springlist.push_back(newspring);
     
   }
    newspring.one=i+(Number-1)*Number;
    newspring.two=i;
    newspring.sticki=sticki;
    newspring.rlen=dy;
    newspring.k=1.0;
    newspring.wlr=0;
    newspring.wud=1;
    springlist.push_back(newspring);
    sticki++; 
  }

    return XY;
}