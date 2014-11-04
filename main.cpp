#include <iostream>
#include <ctime>
#include "random.h"
#include "makemikadonetwork.h"
#include "EnergyandGradients.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <math.h>
#include <functional>

using namespace std;
using namespace Eigen;

const double pi=4.0*atan(1.0);
void getRandom(vector<double> &v)
{
  for (int i = 0; i < v.size(); ++i)
  {
    v[i] = randf();
    cout << v[i] << endl;
  }
}

bool operator<(const stick& first, const stick& second){
  if(first.nr>=second.nr){
    return false;}
    else{ return true;
    }
   }
double inbox(double x,double boxsize){
  if(x<0){
    x=x+boxsize;
  }
  if(x>boxsize){
  x=x-boxsize;
  }
return x;
}

  
int main (int argc,char **argv)
{
  
if(argc>1){
  int SEED=stoi(argv[1]);
  my_random::set_seed(SEED);
}  
  
  
int NumberMikado=150;
double LStick=.20; //Stick Length
double k1=100;
double k2=1;
double rlenshort=.01;
double rlenlong=.001;

vector<stick> mikado=make_sticks(NumberMikado);
vector<stick> mikorig=mikado; //The original set of sticks
//Find the lr wall intercepts and add a ghost to the mikado's
vector<stick> GhostLR = make_ghost_lr(mikado, LStick, mikado.size()); 
    mikado.insert(mikado.end(),GhostLR.begin(),GhostLR.end()); //add the newly found sticks to the existing sticks 
vector<stick> GhostUD=make_ghost_ud(mikado,LStick,mikado.size()); 
if(GhostUD.size()!=0){
vector<stick> GhostLR2=make_ghost_lr(GhostUD,LStick,GhostUD.size());
    mikado.insert(mikado.end(),GhostUD.begin(),GhostUD.end());
if(GhostLR2.size()!=0){
    mikado.insert(mikado.end(),GhostLR2.begin(),GhostLR2.end());
  }
}
std::sort(mikado.begin(),mikado.end());

// for(int i=0;i<mikado.size();i++){
//     cout<<mikado[i].nr<<"\t"<<mikado[i].x<<"\t"<<mikado[i].y<<"\t"<<mikado[i].th<<"\t"<<mikado[i].wlr<<mikado[i].wud<<endl;
// }


vector<connected> Connection=make_connections(mikado,LStick); //Make Connections
vector<connected> Connection2(1); //sives the double elements from the connections
 Connection2[0]=Connection[0];
 for(int i=0;i<Connection.size();i++){
     int flag=1;
    for(int j=0;j<Connection2.size();j++){
        if(Connection[i].first==Connection2[j].first && Connection[i].second==Connection2[j].second){
            flag=0; 
            break;
        }
    }
    if(flag==1){
      Connection2.push_back(Connection[i]);
    } 
}
Connection=Connection2;

//Check for pairs in the connected struct if at least a pair exists (2 points on same mikado)
//then connection[j].recur=1; Else recur =0.
for(int i=0;i<Connection.size()-1;i++){
  for(int j=i+1;j<Connection.size();j++){
    if(Connection[i].first==Connection[j].first){
      Connection[i].recur=1;  Connection[j].recur=1;
    }
   }
 }

//Here we create the nodes, and the springs from the conncection structure that has already
//been made above. 
vector<elonstick> ELONSTICK=sortELEMENTSperMIKADO(Connection);
vector<spring> springlist(0);
vector<node> nodes(0);

        //print here the sticknrs and the coordinates on it. 
        for(int i=0;i<ELONSTICK.size();i++){
            vector<int> nrs=ELONSTICK[i].nr;
            vector<double> pos=ELONSTICK[i].S;
            cout<<"on stick "<<ELONSTICK[i].sticki<<endl;
            for(int j=0;j<nrs.size();j++){
                cout<<nrs[j]<<"\t";
            }
            cout<<endl;
            for(int j=0;j<pos.size();j++){
                cout<<pos[j]<<"\t";
            }
            cout<<endl;
        }

        
 vector<int> order;
 order.push_back(0);
 for(int i=0;i<ELONSTICK.size();i++){
    vector<int> numbervec=ELONSTICK[i].nr;
    
    for(int j=0;j<numbervec.size();j++){
        int flag=0;
        for(int k=0;k<order.size();k++){
            if(numbervec[j]==order[k]) flag=1;
        }
        if(flag==0) order.push_back(numbervec[j]);
    }
 }
     
 std::sort(order.begin(),order.end());
 for(int i=0;i<order.size();i++){
     cout<<"order  "<<i<<"\t"<<order[i]<<endl;
 }
     
 for(int i=0;i<ELONSTICK.size();i++){
  for(int j=0;j<ELONSTICK[i].nr.size();j++){
      for(int k=0;k<order.size();k++){
            if(ELONSTICK[i].nr[j]==order[k]) ELONSTICK[i].nr[j]=k;
      }
  }
}
 
     
     
     
        for(int i=0;i<ELONSTICK.size();i++){
            vector<int> nrs=ELONSTICK[i].nr;
            vector<double> pos=ELONSTICK[i].S;
            cout<<"on stick "<<ELONSTICK[i].sticki<<endl;
            for(int j=0;j<nrs.size();j++){
                cout<<nrs[j]<<"\t";
            }
            cout<<endl;
            for(int j=0;j<pos.size();j++){
                cout<<pos[j]<<"\t";
            }
            cout<<endl;
        }

SpringsAndNodes(ELONSTICK,mikorig,springlist,
                nodes,rlenshort,rlenlong,k1,k2); //Make the springs and Nodes. Input springlist and nodes are (empty vectors)

//Remove all double info
vector<node> singleNodes; 
for(int i=0;i<nodes.size();i++){
  if(nodes[i].number!=nodes[i+1].number){
    node unique=nodes[i];   
    singleNodes.push_back(unique);  
  }
}

cout<<"#############Newnewnodes############"<<endl;
for(int i=0;i<singleNodes.size();i++){
        cout<<singleNodes[i].number<<"\t"<<singleNodes[i].x<<"\t"<<singleNodes[i].y<<endl;
}
cout<<"+++++++++newnewsprings+++++++"<<endl;
for(int i=0;i<springlist.size();i++){
 cout<<springlist[i].one<<"\t"<<springlist[i].two<<"\t wlr:"<<springlist[i].wlr<<"\t wud"<<springlist[i].wud
    <<"\t on stick:\t"<<springlist[i].sticki<<endl;   
}


// //**************MAKE HERE THE PAIR OF SPRINGS

//springpairs is an vector, coding a triplet in each row: 
//springpairs[i][1]= ith triplet first spring
//springpairs[i][2] ith triplet second spring. 
//springpairs[i][3] ith triplet mikado nr.
//the labels for the springs are the indexnumbers of coding the
//springs in the springslist.

vector<vector<int>> springpairs(0);
for(int i=0;i<springlist.size()-1;i++){
    if(springlist[i].sticki==springlist[i+1].sticki){
        vector<int> pair(3);
        pair[0]=i;
        pair[1]=i+1;
        pair[2]=springlist[i].sticki;
        springpairs.push_back(pair);
    }   
}

//****************

//The xy positions
 VectorXd X(singleNodes.size()),Y(singleNodes.size());
 for(int i=0;i<singleNodes.size();i++){
 X(i)=singleNodes[i].x; 
 Y(i)=(singleNodes[i].y);
 }
 for(int i=0;i<singleNodes.size();i++){
 X(i)=inbox(X(i),1.0);
 Y(i)=inbox(Y(i),1.0);
 }
 
 vector<anchor> Anchor;
 anchor Angtemp;
 Angtemp.label=0;
 Angtemp.xpos=.5;
 Angtemp.ypos=.5;
 Anchor.push_back(Angtemp);
 
 VectorXd XY(2*X.size());
 XY<<X,
     Y;
 
int num=XY.size()/2;     
double Energy=Energynetwork(springlist,XY,Anchor);
double EBEND=Ebend(springpairs,springlist,XY);

// int spring7=springpairs[3][0];
// int spring8=springpairs[3][1];
// cout<<"those are springs \t"<<spring7<<"\t and\t"<<spring8<<endl;
// cout<<"on spring \t :"<<spring7<<"\t live coords \t"<<springlist[spring7].one<<"\t and \t"<<springlist[spring7].two<<endl;
// cout<<"on spring \t :"<<spring8<<"\t live coords \t"<<springlist[spring8].one<<"\t and \t"<<springlist[spring8].two<<endl;

// double x12,x13,x14,y12,y13,y14;
// x12=XY(springlist[spring7].one);
// x13=XY(springlist[spring7].two);
// x14=XY(springlist[spring8].two);
// y12=XY(springlist[spring7].one+num);
// y13=XY(springlist[spring7].two+num);
// y14=XY(springlist[spring8].two+num);

// cout<<"x12 \t y12 \t"<<x12<<"\t"<<y12<<endl;
// cout<<"x13 \t y13 \t"<<x13<<"\t"<<y13<<endl;
// cout<<"x14 \t y14 \t"<<x14<<"\t"<<y14<<endl;

// Vector2d v1,v2;
// v1<<x13-x12,y13-y12;
// v2<<x14-x13,y14-y13;

// double t=v1.dot(v2);
// cout<<"teller="<<t<<endl;
// double n=sqrt(v1.dot(v1)*v2.dot(v2));
// double co=t/n;
// cout<<"cos th="<<co<<endl;
// double a=acos(co);
// cout<<"angle="<<a<<endl;
// 
// // 
// //*********************************************************************************
// 
// //Here comes the conjugate gradient
// VectorXd gradE(XY.size());
// VectorXd XYn(XY.size());
// VectorXd gradEn(gradE.size());
// VectorXd s0(gradE.size());
// VectorXd sn(s0.size());
// vector<spring> newsprings(springlist.size());
// double betan;
// 
// 
// 
// gradE=Gradient(springlist,XY,Anchor);
// s0=-gradE;
// 
// 
// ofstream XYfile("conjpoints.txt");
// ofstream springfile("springboundaries.txt");
// int Nit=20;
// for(int i=0;i<Nit;i++)
// {
//   for(int k=0;k<springlist.size();k++){ //loop over all springs to find new values for the borders
//     spring tempspring;
//     tempspring.one=springlist[k].one;
//     tempspring.two=springlist[k].two;
//     tempspring.wlr=0;
//     tempspring.wud=0;
//     double x1=XY(tempspring.one); 
//     double x2=XY(tempspring.two);
//     double y1=XY(tempspring.one+XY.size()/2); 
//     double y2=XY(tempspring.two+XY.size()/2);
//     
//     if((x1>x2)&& abs(x1-x2)>.5){
//       tempspring.wlr=1;
//       }
//     else if((x2>x1)&&abs(x1-x2)>.5){
// 	tempspring.wlr=-1;
//       }
//     if((y1>y2)&& abs(y1-y2)>.5){
//       tempspring.wud=1;
//       }
//     else if((y2>y1)&&abs(y1-y2)>.5){
// 	tempspring.wud=-1;
//       }
//     newsprings[k]=tempspring;
//   }
//   
//   for(int j=0;j<newsprings.size();j++){ //write the newsprings to a file [row wlr ---------- y wud ------]
//   springfile<<newsprings[j].wlr<<"\t"; 
//   }
//   for(int j=0;j<newsprings.size()-1;j++){
//   springfile<<newsprings[j].wud<<"\t";
//   }
//   springfile<<newsprings[newsprings.size()].wud<<endl;
//   
//   for(int j=0;j<XY.size();j++) //write the XY-data to txt
//   {
//    XYfile<<XY(j)<<"\t";
//   }
//   XYfile<<endl;
// 
//   
// //This is the secant method
// double an2=0.01;
// double an1=0;
// double an;
// double tol=.0001;
// do{ 
// an=an1-dEda(XY+an1*s0,Anchor,s0,springlist)*(an1-an2)/(dEda(XY+an1*s0,Anchor,s0,springlist)-
//   dEda(XY+an2*s0,Anchor,s0,springlist));
// an2=an1;
// an1=an;
// }while(abs(an-an1)>tol);
// 
// //Update variables
// XYn=XY+an*s0;
// gradEn=Gradient(springlist,XYn,Anchor);
// betan=gradEn.dot(gradEn)/(gradE.dot(gradE));
// sn=-gradEn+betan*s0;
// s0=sn;
// gradE=gradEn;
// XY=XYn;
// }
// 
// XYfile.close();
// springfile.close();
// 
// FILE *fp = fopen("mikado.txt","w");
//   for(int i=0;i<m.size();i++){
//     fprintf(fp,"%d \t %1.8f \t %1.8f \t %1.8f \t %d \t %d\n",m[i].nr,m[i].x,m[i].y,m[i].th,m[i].wlr,m[i].wud);
//   }
//   fclose(fp);
// FILE *fp2=fopen("nodes.txt","w");
//  for(int i=0;i<singleNodes.size();i++){
//  fprintf(fp2,"%d \t %1.8f \t %1.8f \n",singleNodes[i].number,X(i),Y(i));
//  }
//  fclose(fp2);
FILE *fp3 = fopen("springs.txt","w");
for(int i=0;i<springlist.size();i++){
  fprintf(fp3,"%d \t %d \t %d \t %d \lt %1.8f \t %1.8f \t %d\n",springlist[i].one,springlist[i].two,springlist[i].wlr,springlist[i].wud,springlist[i].rlen, springlist[i].k , springlist[i].sticki);
 }
 fclose(fp3);
FILE *fp4=fopen("mikado1.txt","w");
 for(int i=0;i<mikorig.size();i++){
  fprintf(fp4,"%1.8f \t %1.8f \t %1.8f\n",mikorig[i].x,mikorig[i].y,mikorig[i].th);
 }
 fclose(fp4);

return 0;
 }
