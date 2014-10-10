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
  
int main ()
{
int NumberMikado=100;
double LStick=.25; //Stick Length
double k=1;
double L=.01;

vector<stick> m=make_sticks(NumberMikado);
vector<stick> mikorig=m; //The original set of sticks
//Find the lr wall intercepts and add a ghost to the mikado's
vector<stick> GhostLR = make_ghost_lr(m, LStick, m.size()); 
    m.insert(m.end(),GhostLR.begin(),GhostLR.end()); //add the newly found sticks to the existing sticks 
vector<stick> GhostUD=make_ghost_ud(m,LStick,m.size()); 
if(GhostUD.size()!=0){
vector<stick> GhostLR2=make_ghost_lr(GhostUD,LStick,GhostUD.size());
    m.insert(m.end(),GhostUD.begin(),GhostUD.end());
if(GhostLR2.size()!=0){
    m.insert(m.end(),GhostLR2.begin(),GhostLR2.end());
  }
}
std::sort(m.begin(),m.end());


vector<connected> Connection=make_connections(m,LStick); //Make Connections
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


//Misschien gooi ik deze weg! 
// Check for pairs in the connected struct if at least a pair exists (2 points on same mikado)
// then connection[j].recur=1; Else recur =0.
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
vector<spring> springlist;
vector<node> nodes;
SpringsAndNodes(ELONSTICK,mikorig,springlist,nodes); //Make the springs and Nodes. Input springlist and nodes are (empty vectors)

//Remove all double info
vector<node> singleNodes; 
for(int i=0;i<nodes.size();i++){
  if(nodes[i].number!=nodes[i-1].number){
    node unique=nodes[i];   
    singleNodes.push_back(unique);  
  }
}

//Now put the nodes and spring elements in a ascending logical ordering from 0 to #numberofsprings
for(int i=0;i<singleNodes.size(); i++){
  if(singleNodes[i].number!=i){
    for(int j=0;j<springlist.size();j++){
      if(springlist[j].one==singleNodes[i].number){
	springlist[j].one=i;
      }
      if(springlist[j].two==singleNodes[i].number){
	springlist[j].two=i;
     }
   }
   singleNodes[i].number=i;
 }
}

//The xy positions
VectorXd X(singleNodes.size()),Y(singleNodes.size());
for(int i=0;i<singleNodes.size();i++){
X(i)=singleNodes[i].x; //dit is niet wat ik wil. Schrijf even een func. 
Y(i)=(singleNodes[i].y);
}
for(int i=0;i<singleNodes.size();i++){
X(i)=inbox(X(i),1.0);
Y(i)=inbox(Y(i),1.0);
}


VectorXd XY(2*X.size());
XY<<X,
    Y;
double Energy=Energynetwork(springlist,XY,k,L);

//Initiate the gradient.
VectorXd gradE(2*X.size());
VectorXd XYn(XY.size());
VectorXd gradEn(gradE.size());
VectorXd s0(gradE.size());
VectorXd sn(s0.size());
double betan;

gradE=Gradient(springlist,XY,1,.01);
s0=-gradE;

//Here comes the conjugate gradient

ofstream XYfile("test.txt");

for(int i=0;i<10;i++)
{
  for(int k=0;k<XY.size();k++)
  {
    XYfile<<XY(k)<<"\t";
  }
  XYfile<<endl;
  
//This is the secand method
double an2=0.01;
double an1=0;
double an;
double tol=.0001;
do{
an=an1-dEda(XY+an1*s0,s0,springlist,k,L)*(an1-an2)/(dEda(XY+an1*s0,s0,springlist,k,L)-dEda(XY+an2*s0,s0,springlist,k,L));
an2=an1;
an1=an;
}while(abs(an-an1)>tol);

XYn=XY+an*s0;
gradEn=Gradient(springlist,XYn,k,L);
betan=gradEn.dot(gradEn)/(gradE.dot(gradE));
sn=-gradEn+betan*s0;
s0=sn;
gradE=gradEn;
XY=XYn;

  for(int j=0;j<XY.size();j++){
    cout<<XY[j]<<endl; 
  }
}

XYfile.close();





FILE *fp = fopen("mikado.txt","w");
  for(int i=0;i<m.size();i++){
    fprintf(fp,"%d \t %1.8f \t %1.8f \t %1.8f \t %d \t %d\n",m[i].nr,m[i].x,m[i].y,m[i].th,m[i].wlr,m[i].wud);
  }
  fclose(fp);
 
 FILE *fp4=fopen("mikado1.txt","w");
 for(int i=0;i<mikorig.size();i++){
  fprintf(fp,"%1.8f \t %1.8f \t %1.8f\n",mikorig[i].x,mikorig[i].y,mikorig[i].th);
 }
 fclose(fp4);
 
 FILE *fp2=fopen("nodes.txt","w");
 for(int i=0;i<singleNodes.size();i++){
 fprintf(fp2,"%d \t %1.8f \t %1.8f \n",singleNodes[i].number,X(i),Y(i));
 }
 fclose(fp2);
 
 
 FILE *fp3 = fopen("springs.txt","w");
for(int i=0;i<springlist.size();i++){
  fprintf(fp3,"%d \t %d \t %d \t %d\n",springlist[i].one,springlist[i].two,springlist[i].wlr,springlist[i].wud);
 }
 fclose(fp3);

    return 0;
 }