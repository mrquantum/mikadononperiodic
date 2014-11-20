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
  
int main (int argc,char **argv)
{
  
if(argc>1){
  int SEED=stoi(argv[1]);
  my_random::set_seed(SEED);
}  
  
//System parameters  
int Nit=150;  
double tolE=.001;
int NumberMikado=200;
double LStick=.4; //Stick Length
double k1=.05;
double k2=.1;
double kappa=0.000006;
double rlenshort=.0001;
double rlenlong=.01;

vector<stick> mikado;
vector<stick> mikorig;
makeSticks(mikado,mikorig,NumberMikado,LStick);



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
makeSpringsAndNodes(ELONSTICK,mikorig,springlist,nodes,rlenshort,rlenlong,k1,k2); 

//Remove all double info
vector<node> singleNodes; 
for(std::size_t i=0;i<nodes.size();i++){
  if(nodes[i].number!=nodes[i+1].number){
    node unique=nodes[i];   
    singleNodes.push_back(unique);  
  }
}


//**************MAKE HERE THE PAIR OF SPRINGS
//springpairs is an vector, coding a triplet in each row: 
//springpairs[i][1]= ith triplet first spring
//springpairs[i][2] ith triplet second spring. 
//springpairs[i][3] ith triplet mikado nr.
//the labels for the springs are the indexnumbers of coding the
//springs in the springslist.
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
    Y(i)=singleNodes[i].y;
 }
 for(std::size_t i=0;i<singleNodes.size();i++){
    X(i)=inbox(X(i),1.0);
    Y(i)=inbox(Y(i),1.0);
 }
  
 VectorXd XY(2*X.size());
 XY<<X,Y;
 
 
int num=XY.size()/2;     
double Energy=Energynetwork(springlist,XY);
double EBEND=Ebend(springpairs,springlist,XY,kappa);


//Here comes the conjugate gradient
    VectorXd gradE(XY.size());
    VectorXd XYn(XY.size());
    VectorXd gradEn(gradE.size());
    VectorXd s0(gradE.size());
    VectorXd sn(s0.size());
    VectorXd dXY(XY.size());
    double maxdispl=.01;
    double betan;
    gradE=Gradient(springlist,XY)+
          gradEbend(springpairs,springlist,XY,kappa);
    s0=-gradE;
 
 
    ofstream XYfile("conjpoints.txt");
    ofstream EFile("Energy.txt");
    EFile<<Energy<<"\t"<<EBEND<<endl;
 
    //The loop of the conj grad method
    int conjsteps=0;
    do{
        conjsteps++;

        for(int j=0;j<XY.size();j++){ //write the XY-data to txt
            XYfile<<XY(j)<<"\t";
        }
        XYfile<<endl;
 
        double a1,a2,root;
        a1=0;
        a2=.1;
        //doBracketfind(a1,a2,XY,s0,springlist,springpairs,kappa);
        //cout<<"The bracket:  "<<a1<<"  "<<a2<<endl;
        //doBisection(a1,a2,root,XY,s0,springlist,springpairs,kappa);
        //doFalsePosition(a1,a2,root,XY,s0,springlist,springpairs,kappa);
        doSecant(root,XY,s0,springlist,springpairs,kappa);
        cout<<conjsteps<<"\t"<<root<<endl;
        double an=root;  
  
        //Update variables
        double adeptsteps;
        //adeptsteps=sqrt(1+sqrt(s0.dot(s0)));
        //adeptsteps=1.0/adeptsteps;
        adeptsteps=1.0;
        
        dXY=an*adeptsteps*s0;
        for(int q=0;q<dXY.size();q++){
            if(dXY(q)>maxdispl) dXY(q)=maxdispl*sgn(dXY(q));
        }
        XYn=XY+dXY;
        //XYn=XY+an*adeptsteps*s0;
        gradEn=Gradient(springlist,XYn)+gradEbend(springpairs,springlist,XYn,kappa);
        //betan=gradEn.dot(gradEn)/(gradE.dot(gradE));
        betan=(gradEn-gradE).dot(gradEn)/gradE.dot(gradE);
 
 
        if(conjsteps%10==0){ 
            sn=-gradEn;}
        else{
            sn=-gradEn+betan*s0;
        }
        s0=sn;
        gradE=gradEn;
        XY=XYn;
        Energy=Energynetwork(springlist,XY);
        EBEND=Ebend(springpairs,springlist,XY,kappa);    
        EFile<<Energy<<"\t"<<EBEND<<endl;

    }while(conjsteps<Nit && sqrt(gradE.dot(gradE))>tolE);

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
