#include <iostream>
#include <ctime>
#include "random.h"
#include "makemikadonetwork.h"
#include "EnergyandGradients.h"
#include "minimizers.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Sparse>

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <math.h>
#include <functional>
#include "importparam.h"
#include "prestress.h"

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
  
int get_next_script_tag(FILE *in, char *buf)
{
        if(fscanf(in, "%*[^<]") != 0) {
                printf("White space removal error in get_next_script_tag(). Something odd has happened: this error should never occur...\n");
        }
        if(fscanf(in, "<%255[^>]>", buf) != 1) {
                printf("Error reading tag in script file.\n");
                return 0;
        }

        return 1;
} 
    
int init(param &Mikadoparameters ,const char *FILENAME)
{
    FILE *paramfile;
    paramfile=fopen(FILENAME,"r");
    if(!paramfile){
     cout<<"Error opening paramfile"<<endl;
     return -1;     
    }
    const int buffersize=101;
    char buf[buffersize];
//     if(get_next_script_tag(paramfile,buf)==0){
//         cout<<"Error reading script file"<<endl;   
//         fclose(paramfile);
//         return -1;
//    }

for(;;) {
        if(get_next_script_tag(paramfile, buf) == 0) {
                printf("\tError when reading tag in scriptfile\n");
                fclose(paramfile);
                return -1;
        }
//         if(feof(script_file)) {
//                 FFEA_error_text();
//                 printf("\tReached end of file before end of <param> block\n");
//                 fclose(script_file);
//                 return FFEA_ERROR;
//         }

        if(strcmp(buf, "END") == 0) {
                break;
        }

        if(Mikadoparameters.parse_param_assignment(buf) == -1) {
                printf("Parameter assignment failed.\n");
                fclose(paramfile);
                return -1;
        }
}    
 return 0;   
}
 
int main (int argc,char **argv)
{
  
if(argc>1){
  int SEED=stoi(argv[1]);
  my_random::set_seed(SEED);
}  
//my_random::set_seed(0);
//System parameters  

param Mikadoparameters;
//char parameterfile = "params.txt";

init(Mikadoparameters,"params.txt");

double ESTRETCH,EBEND,ETOT;
vector<int> order;
vector<stick> mikado;
vector<stick> mikorig;
vector<connected> Connection;
vector<elonstick> ELONSTICK;
vector<spring> springlist(0);
vector<node> nodes(0);
vector<node> singleNodes; 
vector<vector<int>> springpairs(0);
double lenGrad;
//vector<double> alpha(100);

int Nit=Mikadoparameters.Nit;  
double tolE=Mikadoparameters.tolE;
int NumberMikado=Mikadoparameters.NumberMikado;
double LStick=Mikadoparameters.LStick; //Stick Length
double k1=Mikadoparameters.k1;
double k2=Mikadoparameters.k2;
double kappa=Mikadoparameters.kappa;
double rlenshort=Mikadoparameters.rlenshort;
double rlenlong=Mikadoparameters.rlenlong;

ofstream mikadofile("mikado.txt"); 
ofstream nodefile("nodes.txt");
ofstream springfile("springs.txt");
ofstream XYfile("conjpoints.txt");
ofstream EFile("Energy.txt");
ofstream shearcoordinates("shearcoordinates.txt");
ofstream shearenergy("shearenergy.txt");
ofstream dEdafile("dEda.txt");

makeSticks(mikado,mikorig,NumberMikado,LStick);
    //write sticks to mikado.txt
for(int i=0;i<mikado.size();i++){
    mikadofile<<mikado[i].nr<<"\t"<<mikado[i].x<<"\t"<<mikado[i].y<<"\t"<<mikado[i].th<<"\t"<<mikado[i].wlr<<
    mikado[i].wud<<endl;
} mikadofile.close();

makeConnections(Connection,mikado,LStick);
     //Here we create the nodes, and the springs from the conncection structure that has already
     //been made above. 
sortELEMENTSperMIKADO(ELONSTICK,Connection);
orderElonstick(order,ELONSTICK); 
 
     //Make the springs and Nodes. Input springlist and nodes are (empty vectors)
makeSpringsAndNodes(ELONSTICK,mikorig,springlist,nodes,rlenshort,rlenlong,k1,k2); 
    
//write sticks to springs.txt
for(int i=0;i<springlist.size();i++){
    springfile<<springlist[i].one<<"\t"
              <<springlist[i].two<<"\t"
              <<springlist[i].wlr<<"\t"
              <<springlist[i].wud<<"\t"
              <<springlist[i].rlen<<"\t"
              <<springlist[i].k<<"\t"
              <<springlist[i].sticki<<endl;
}               springfile.close();

     //Remove all double info
for(std::size_t i=0;i<nodes.size();i++){
   if(nodes[i].number!=nodes[i+1].number){
     node unique=nodes[i];   
     singleNodes.push_back(unique);  
   }
}
//**************MAKE HERE THE PAIR OF SPRINGS
    makeSpringpairs(springpairs,springlist);
        for(int i=0;i<singleNodes.size();i++){
        }

    VectorXd X(singleNodes.size()),Y(singleNodes.size());
    VectorXd XY(2*X.size());
    VectorXd gradE(XY.size());
    VectorXd XYn(XY.size());
    VectorXd XYcopy(XY.size());
    VectorXd gradEn(gradE.size());
    VectorXd s0(gradE.size());


    for(int i=0;i<singleNodes.size();i++){
        nodefile<<singleNodes[i].number<<"\t"<<singleNodes[i].x<<"\t"<<singleNodes[i].y<<endl;
    } nodefile.close();


//The xy positions
    for(std::size_t i=0;i<singleNodes.size();i++){
        X(i)=singleNodes[i].x; 
        Y(i)=singleNodes[i].y;
    }
    for(std::size_t i=0;i<singleNodes.size();i++){
        X(i)=inbox(X(i),1.0);
        Y(i)=inbox(Y(i),1.0);
    }

    XY<<X,Y;
    double g11=1.0;
    double g12=0.0;
    double g22=1.0;
    int relax=0;
    gradE=Gradient(springlist,XY,g11,g12,g22)+gradEbendn(springpairs,springlist,XY,g11,g12,g22,kappa);
    s0=-gradE; 
    double lengradBegin=sqrt(s0.dot(s0));
    double lengradconjstep;
    cout<<lengradBegin<<endl;
    
//     vector<double> alpha(2000);
//     for(int i=0; i<alpha.size();i++){
//      alpha[i]=-1+i*.002;        
//     dEdafile<<alpha[i]<<"\t";
//     }
//     dEdafile<<endl;
    
//Initial relaxation
    
    do{
//         for(int q=0;q<alpha.size();q++){
//         dEdafile<<dEda(XY+alpha[q]*s0,s0,springlist,springpairs,kappa,g11,g12,g22)<<"\t";        
//         }
//         dEdafile<<endl;
        doConjStep(XY,s0,gradE,springlist,springpairs,kappa,relax,g11,g12,g22);
        relax++;
        lengradconjstep=sqrt(gradE.dot(gradE));
        cout<<relax<<"  "<<lengradconjstep/lengradBegin<<endl;
    }while(relax<Nit && lengradBegin<2000000*lengradconjstep);

    double P=prestress(springlist,springpairs,XY,kappa,Nit);
    cout<<"THE PRESTRESS IS "<<P<<endl;




// //Shearproperties
    double deltaAngle=.005;
    double angle=0; 
//This is the shearingloop
    for(int k=0;k<175;k++){ 
        double g11=1;
        double g12=tan(angle);
        double g21=g12;
        double g22=1+tan(angle)*tan(angle);

//Here comes the conjugate gradient
        gradE=Gradient(springlist,XY,g11,g12,g22)+gradEbendn(springpairs,springlist,XY,g11,g12,g22,kappa);
        s0=-gradE;  

        int conjsteps=0;
        double root1=0;
        double ESTRETCH=Energynetwork(springlist,XY,g11,g12,g22);
        double EBEND=Ebend(springpairs,springlist,XY,g11,g12,g22,kappa);
        double ETOT=ESTRETCH+EBEND; 
        EFile<<ESTRETCH<<"\t"<<EBEND<<"\t"<<ETOT<<"\t"<<0<<endl;

//loop of the cg-method
        lengradBegin=sqrt(gradE.dot(gradE));
        do{
        //for(int j=0;j<XY.size();j++){ //write the XY-data to txt
        //XYfile<<XY(j)<<"\t";
        //} XYfile<<endl;
        conjsteps++;
        cout<<"\r"<<conjsteps<<"   "<<flush;
        doConjStep(XY,s0,gradE,springlist,springpairs,kappa,conjsteps,g11,g12,g22);         
        ESTRETCH=Energynetwork(springlist,XY,g11,g12,g22);
        EBEND=Ebend(springpairs,springlist,XY,g11,g12,g22,kappa);    
        ETOT=ESTRETCH+EBEND;    
        lengradconjstep=sqrt(gradE.dot(gradE));
        //Write the Energy to a txt-file.
        EFile<<ESTRETCH<<"\t"<<EBEND<<"\t"<<ETOT<<"\t"<<lenGrad<<endl; 
    }while(conjsteps<Nit && lengradBegin<2000000*lengradconjstep);

    for(int ii=0;ii<XY.size();ii++){
        shearcoordinates<<XY(ii)<<"\t";
    }
    shearenergy<<angle<<"\t"<<ETOT<<endl;

        if(k<25||(k>=75&&k<125)){
    angle=angle+deltaAngle;
    }
    if((k>=25&&k<75)||k>=125){
        angle=angle-deltaAngle;
    }
    angle=angle+deltaAngle;
    shearcoordinates<<endl;
}




XYfile.close();
EFile.close();
shearcoordinates.close();
shearenergy.close();
//dEdafile.close();


return 0;
}

 
 