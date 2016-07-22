#include <iostream>
#include <ctime>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <array>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "writefunctions.h"
#include "structs.h"

using namespace std;
using namespace Eigen;

void combineElementsOnMikado(vector<connected> &connection,vector<elonstick> &elements){

  vector<connected> con=connection;
  int nr;
  nr=((con[0].first>con[0].second)?con[0].first:con[0].second);
  int newnr;
  //Make a vector of elements with the lenght of the max elements nr.
  
  for(int i=1; i<con.size();i++){
   newnr=((con[i].first>con[i].second)?con[i].first:con[i].second);
   nr=((newnr>nr)?newnr:nr);
  }
  vector<elonstick> Elements(nr+1);
  for(int i=0;i<Elements.size();i++){
    Elements[i].sticki=i; 
  }

  //ads some elements double? Check tomorrow
  int first,second,connr;
  double s1,s2;
  int type;
  int backgroundspring[2];
  int bspringn1, bspringn2;
  vector<int> bdeleted;
  
 // cout<<"first  second  s1      s2      type    nr      spr1    spr2"<<endl;
  int nrELE=0;
  while(con.size()>0){
    first=con[0].first;
    second=con[0].second;
    s1=con[0].s1;
    s2=con[0].s2;
    type=con[0].type;
    connr=con[0].nrCon;
    bspringn1=con[0].backgroundspring[0];
    bspringn2=con[0].backgroundspring[1];
    
   // cout<<first<<"\t"<<second<<"\t"<<s1<<"\t"<<s2<<"\t"<<type<<"\t"<<connr<<"\t"<<bspringn1<<"\t"<<bspringn2<<endl;
    int backgroundspring[2];
    
    if(type==0){
      std::array<int,2> minuslist;
      minuslist[0] = -1;
      minuslist[1] = -1;
      Elements[first].nr.push_back(nrELE);
      Elements[first].S.push_back(s1);
      Elements[first].type.push_back(type);
      Elements[first].backgroundspringn.push_back(minuslist);
      
      Elements[second].nr.push_back(nrELE);
      Elements[second].S.push_back(s2);
      Elements[second].type.push_back(type);
      Elements[second].backgroundspringn.push_back(minuslist);

      nrELE++;
      con.erase(con.begin());
      
	for(int j=0;j<con.size();j++){
	  if((con[j].first==second && con[j].second==first) ||(con[j].first==first && con[j].second==second)){
	      bdeleted.push_back(j); //make a list with elements to be deleted later on
	  }

	}
      while(bdeleted.size()>0){ //list of elements 2 be deleted from the con-vector from back to front bacause it does
	//not mess w.t. indexation of the array
	con.erase(con.begin()+bdeleted[bdeleted.size()-1]);
	bdeleted.pop_back();
      }
    }
    
       
    if(type==1){
      std::array<int,2> minuslist;
      Elements[first].nr.push_back(nrELE);
      Elements[first].S.push_back(s1);
      Elements[first].type.push_back(type);
      minuslist[0] = bspringn1;
      minuslist[1] = bspringn2;
      Elements[first].backgroundspringn.push_back(minuslist);
      nrELE++;
      con.erase(con.begin());
    }
  }

   elements=Elements; 
    

  
}
