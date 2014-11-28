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
//my_random::set_seed(0);
//System parameters  
int Nit=500;  
double tolE=.001;
int NumberMikado=100;
double LStick=.4; //Stick Length
double k1=.05;
double k2=.1;
double kappa=.000007;
double rlenshort=.0001;
double rlenlong=.01;

vector<stick> mikado;
vector<stick> mikorig;
makeSticks(mikado,mikorig,NumberMikado,LStick);
vector<connected> Connection;
makeConnections(Connection,mikado,LStick);
//Here we create the nodes, and the springs from the conncection structure that has already
//been made above. 
vector<elonstick> ELONSTICK;
sortELEMENTSperMIKADO(ELONSTICK,Connection);
vector<int> order;
orderElonstick(order,ELONSTICK); 

//Make the springs and Nodes. Input springlist and nodes are (empty vectors)
vector<spring> springlist(0);
vector<node> nodes(0);
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
vector<vector<int>> springpairs(0);
makeSpringpairs(springpairs,springlist);


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
 
 
double ESTRETCH=Energynetwork(springlist,XY);
double EBEND=Ebend(springpairs,springlist,XY,kappa);
double ETOT=ESTRETCH+EBEND;

//Here comes the conjugate gradient
    VectorXd gradE(XY.size());
    VectorXd XYn(XY.size());
    VectorXd XYcopy(XY.size());
    VectorXd gradEn(gradE.size());
    VectorXd s0(gradE.size());
    double betan;
    double Estr, Eben, Etot, lenGrad;

    gradE=Gradient(springlist,XY)+gradEbend(springpairs,springlist,XY,kappa);
    s0=-gradE; 
     
    ofstream XYfile("conjpoints.txt");
    ofstream EFile("Energy.txt");
    ofstream dEdafile("dEda.txt");
    ofstream Rootalpha("rootalpha.txt");

    vector<double> alpha(200);
    for(int i=0;i<200;i++){
        alpha[i]=-.1+0.001*i;
        dEdafile<<alpha[i]<<"\t";    
    }
    dEdafile<<endl;
    
    EFile<<ESTRETCH<<"\t"<<EBEND<<"\t"<<ETOT<<"\t"<<0<<endl;
 
    //The loop of the conj grad method
    int conjsteps=0;
    double root1;
    double ETOT0;
    ESTRETCH=Energynetwork(springlist,XY);
    EBEND=Ebend(springpairs,springlist,XY,kappa);
    ETOT0=ESTRETCH+EBEND;
    
    
    
 
    //loop of the cg-method
    do{
       for(int j=0;j<XY.size();j++){ //write the XY-data to txt
            XYfile<<XY(j)<<"\t";
        } XYfile<<endl;
      
      conjsteps++;
      cout<<conjsteps<<endl;
    //doSteepestDescent(XY,s0,gradE,springlist,springpairs,root,kappa);
    
      for(int k=0;k<200;k++){
      dEdafile<<dEda(XY+alpha[k]*s0,s0,springlist,springpairs,kappa)<<"\t";  
    }
    dEdafile<<endl;
      
      
      doConjStep(XY,s0,gradE,springlist,springpairs,root1,kappa,conjsteps);         
      Rootalpha<<root1<<"\t";
      
      
      ESTRETCH=Energynetwork(springlist,XY);
      EBEND=Ebend(springpairs,springlist,XY,kappa);    
      ETOT=ESTRETCH+EBEND;    
      lenGrad=sqrt(gradE.dot(gradE));
    EFile<<ESTRETCH<<"\t"<<EBEND<<"\t"<<ETOT<<"\t"<<lenGrad<<endl; //Write the Energy to a txt-file.
    }while(conjsteps<Nit && lenGrad>tolE);

    XYfile.close();
    EFile.close();
    dEdafile.close();
    Rootalpha.close();
    
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

 
 