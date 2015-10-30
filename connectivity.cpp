#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include <vector>

using namespace std;

struct connect{
    int nodenr;
    std::vector<int> connectedto;
};


double Connectivity(vector<spring> &springlist)
{   
    
    //Check the min, max elements of the springlist in order to make 
    //a list with nodes below
    int min,max;
    if(springlist[0].one<springlist[0].two){
        min=springlist[0].one;
        max=springlist[0].two;
    }else{
        min=springlist[0].one;
        max=springlist[0].two;
    }
    for(int i=1;        i<springlist.size();   i++){
       if(springlist[i].one<min){
           min=springlist[i].one;
       }
       if(springlist[i].two<min){
           min=springlist[i].two;
       }
       if(springlist[i].one>max){
           max=springlist[i].one;
       } 
       if(springlist[i].two>max){
           max=springlist[i].two;
       }
    }
    
    vector<connect> convec;
    connect node;
    for(int i=min       ;i<max+1;       i++){
        node.nodenr=i;
        convec.push_back(node);       
    }
    //fill the connection vector w. the connections
    for(int i=0;        i<springlist.size();    i++){
        convec[springlist[i].one].connectedto.push_back(springlist[i].two);
        convec[springlist[i].two].connectedto.push_back(springlist[i].one);
    }
    
    int totNRconnections=0;
    for(int i=0;        i<convec.size();        i++){
        totNRconnections+=convec[i].connectedto.size();
    }
    
    double totNRconnectionsd=static_cast<double>(totNRconnections);
    
    cout<<"min "<<min<<" max  "<<max<<"  size  "<<convec.size()<<endl;
    return totNRconnectionsd/convec.size();
    
    
}