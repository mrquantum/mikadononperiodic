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
  for (std::size_t i = 0; i < v.size(); ++i)
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
  
  
int NumberMikado=300;
double LStick=.3; //Stick Length
double k1=1;
double k2=1;
//double kappa=0.000001;
double kappa=0.00001;
double rlenshort=.0001;
double rlenlong=.01;

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


vector<connected> Connection=make_connections(mikado,LStick); //Make Connections
vector<connected> Connection2(1); //sives the double elements from the connections
Connection2[0]=Connection[0];
for(std::size_t i=0;i<Connection.size();i++){
    int flag=1;
    for(std::size_t j=0;j<Connection2.size();j++){
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
for(std::size_t i=0;i<Connection.size()-1;i++){
  for(std::size_t j=i+1;j<Connection.size();j++){
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
        
 vector<int> order;
 order.push_back(0);
 for(std::size_t i=0;i<ELONSTICK.size();i++){
    vector<int> numbervec=ELONSTICK[i].nr;
    for(std::size_t j=0;j<numbervec.size();j++){
        int flag=0;
        for(std::size_t k=0;k<order.size();k++){
            if(numbervec[j]==order[k]) flag=1;
        }
        if(flag==0) order.push_back(numbervec[j]);
    }
 }
     
 std::sort(order.begin(),order.end());
     
 for(std::size_t i=0;i<ELONSTICK.size();i++){
  for(std::size_t j=0;j<ELONSTICK[i].nr.size();j++){
      for(std::size_t k=0;k<order.size();k++){
            if(ELONSTICK[i].nr[j]==order[k]) ELONSTICK[i].nr[j]=k;
      }
  }
}

//Make the springs and Nodes. Input springlist and nodes are (empty vectors)
SpringsAndNodes(ELONSTICK,mikorig,springlist,nodes,rlenshort,rlenlong,k1,k2); 

//Remove all double info
vector<node> singleNodes; 
for(std::size_t i=0;i<nodes.size();i++){
  if(nodes[i].number!=nodes[i+1].number){
    node unique=nodes[i];   
    singleNodes.push_back(unique);  
  }
}



/*
**************MAKE HERE THE PAIR OF SPRINGS
springpairs is an vector, coding a triplet in each row: 
springpairs[i][1]= ith triplet first spring
springpairs[i][2] ith triplet second spring. 
springpairs[i][3] ith triplet mikado nr.
the labels for the springs are the indexnumbers of coding the
springs in the springslist.
*/

vector<vector<int>> springpairs(0);
for(std::size_t i=0;i<springlist.size()-1;i++){
    if(springlist[i].sticki==springlist[i+1].sticki){
        vector<int> pair(3);
        pair[0]=i;
        pair[1]=i+1;
        pair[2]=springlist[i].sticki;
        springpairs.push_back(pair);
    }   
}


//The xy positions
 VectorXd X(singleNodes.size()),Y(singleNodes.size());
 for(std::size_t i=0;i<singleNodes.size();i++){
 X(i)=singleNodes[i].x; 
 Y(i)=(singleNodes[i].y);
 }
 for(std::size_t i=0;i<singleNodes.size();i++){
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
double EBEND=Ebend(springpairs,springlist,XY,kappa);

//Here comes the conjugate gradient
 VectorXd gradE(XY.size());
 VectorXd XYn(XY.size());
 VectorXd gradEn(gradE.size());
 VectorXd s0(gradE.size());
 VectorXd sn(s0.size());
 vector<spring> newsprings(springlist.size());
 double betan;
 gradE=Gradient(springlist,XY,Anchor)+gradEbend(springpairs,springlist,XY,kappa);
 s0=-gradE;
 
 ofstream XYfile("conjpoints.txt");
 ofstream EFile("Energy.txt");
 EFile<<Energy<<"\t"<<EBEND<<endl;
 
  //The loop of the conj grad method
 int Nit=50;
 for(int i=0;i<Nit;i++)
 {

for(int j=0;j<XY.size();j++){ //write the XY-data to txt
    XYfile<<XY(j)<<"\t";
   }
   XYfile<<endl;
 
   
 //This is the secant method
 double an2=0.0;
 double an1=0.000000001;
 double an;
 double tol=0.0000000001;
 int q=0; 
 double dEda2,dEda1;
 
dEda2=dEda(XY+an2*s0,Anchor,s0,springlist,springpairs,kappa);

 do{ 
    dEda1=dEda(XY+an1*s0,Anchor,s0,springlist,springpairs,kappa);
    an=an1-dEda1*(an1-an2)/(dEda1-dEda2);
    an2=an1;
    an1=an;
    dEda2=dEda1;   
    q++;
 }while(q<50 && abs(an2-an1)>tol);
 cout<<i<<"\t"<<q<<endl;
 //Update variables
 double adeptsteps;
 //adeptsteps=sqrt(1+s0.dot(s0));
 //adeptsteps=1.0/adeptsteps;
 adeptsteps=1.0;
 XYn=XY+an*adeptsteps*s0;
 gradEn=Gradient(springlist,XYn,Anchor)+gradEbend(springpairs,springlist,XYn,kappa);
 //betan=gradEn.dot(gradEn)/(gradE.dot(gradE));
 betan=(gradEn-gradE).dot(gradEn)/gradE.dot(gradE);
 
 
 //Rest the gradient
 if(i%10==0){ 
     sn=-gradEn;}
 else{
    sn=-gradEn+betan*s0;
 }
s0=sn;
gradE=gradEn;
XY=XYn;
Energy=Energynetwork(springlist,XY,Anchor);
EBEND=Ebend(springpairs,springlist,XY,kappa);    
EFile<<Energy<<"\t"<<EBEND<<endl;

}

 XYfile.close();
 EFile.close();
 
 
FILE *fp = fopen("mikado.txt","w");
  for(std::size_t i=0;i<mikado.size();i++){
    fprintf(fp,"%d \t %1.8f \t %1.8f \t %1.8f \t %d \t %d\n",mikado[i].nr,mikado[i].x,mikado[i].y,
            mikado[i].th,mikado[i].wlr,mikado[i].wud);
  }
  fclose(fp);
FILE *fp2=fopen("nodes.txt","w");
 for(std::size_t i=0;i<singleNodes.size();i++){
 fprintf(fp2,"%d \t %1.8f \t %1.8f \n",singleNodes[i].number,X(i),Y(i));
 }
 fclose(fp2);
FILE *fp3 = fopen("springs.txt","w");
for(std::size_t i=0;i<springlist.size();i++){
  fprintf(fp3,"%d \t %d \t %d \t %d \t %1.8f \t %1.8f \t %d\n",springlist[i].one,springlist[i].two,
          springlist[i].wlr,springlist[i].wud,springlist[i].rlen, springlist[i].k , springlist[i].sticki);
 }
fclose(fp3);

FILE *fp4=fopen("mikado1.txt","w");
 for(std::size_t i=0;i<mikorig.size();i++){
  fprintf(fp4,"%1.8f \t %1.8f \t %1.8f\n",mikorig[i].x,mikorig[i].y,mikorig[i].th);
 }
 fclose(fp4);

return 0;
 }
