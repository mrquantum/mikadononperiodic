#include <iostream>
#include <ctime>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "writefunctions.h"
#include "structs.h"

using namespace std;
using namespace Eigen;

void  combineElementsOnMikado(vector<connected> &connection,vector<elonstick> &elements){

 
  vector<connected> con=connection;
  int nr;
  
  nr=((con[0].first>con[0].second)?con[0].first:con[0].second);
  
  int newnr;
  for(int i=1; i<con.size();i++){
   newnr=((con[i].first>con[i].second)?con[i].first:con[i].second);
   nr=((newnr>nr)?newnr:nr);
   cout<<"nr"<<nr<<endl;
  }
  
  
  
  vector<elonstick> Elements(nr+1);
  for(int i=0;i<Elements.size();i++){
    Elements[i].sticki=i; 
    cout<<i<<"i"<<endl;
  }
  
  
  int first,second,connr;
  double s1,s2;
  int type;
  int backgroundspring[2];
  
  while(con.size()>0){
    first=connection[0].first;
    second=connection[0].second;
    s1=connection[0].s1;
    s2=connection[0].s2;
    type=connection[0].type;
    connr=connection[0].nrCon;
    cout<<"))()("<<connr<<endl;
    
    if(type==0){
      Elements[first].nr.push_back(connr);
      Elements[first].S.push_back(s1);
      Elements[first].type.push_back(type);
      
      Elements[second].nr.push_back(connr);
      Elements[second].S.push_back(s2);
      Elements[second].type.push_back(type);
    
      //con.erase(con.begin());
      for(int j=0;j<con.size();j++){
	if(con[j].first==second && con[j].second==first){
	  con.erase(con.begin()+j); 
	}
      }
      con.erase(con.begin());
    }
    
    if(type==1){
      Elements[first].nr.push_back(connr);
      Elements[first].S.push_back(s1);
      Elements[first].type.push_back(type);
      
      con.erase(con.begin());
      
    }
    
    
  }
    
    
   elements=Elements; 
    
    
  
  
  
}