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
#include "BendingEnergy.h"
#include "BendingGrad.h"

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

makeSticks(mikado,mikorig,NumberMikado,LStick);

//write sticks to mikado.txt
for(int i=0;i<mikado.size();i++){
    mikadofile<<mikado[i].nr<<"\t"<<mikado[i].x<<"\t"<<mikado[i].y<<"\t"<<mikado[i].th<<"\t"<<mikado[i].wlr<<
    mikado[i].wud<<endl;
} mikadofile.close();

makeConnections(Connection,mikado,LStick); //Here we create the nodes, and the springs from the conncection structure that has already
                                           //been made above. 
sortELEMENTSperMIKADO(ELONSTICK,Connection);
orderElonstick(order,ELONSTICK); 
makeSpringsAndNodes(ELONSTICK,mikorig,springlist,nodes,rlenshort,rlenlong,k1,k2);//Make the springs and Nodes. 
                                                                                //Input springlist and nodes are (empty vectors)
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

//MAKE HERE THE PAIR OF SPRINGS
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

//Write the original nodes to .../nodes.txt
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

//Shearing steps
double angle=0.0; 
double deltaboxdx=0.01;
double boxdx=0;

for(int k=0;k<100;k++){
double g11=1;
double g12=boxdx;
double g22=1+boxdx*boxdx;
       angle=atan(boxdx);

 double EBEND=Ebending(springpairs,springlist,XY,kappa,g11,g12,g22);
 double ESTRETCH=Energynetwork(springlist,XY,g11,g12,g22);
double ETOT=ESTRETCH+EBEND;
 
//Here comes the conjugate gradient
gradE=HarmonicGradient(springlist,XY,g11,g12,g22)+BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22);
s0=-gradE;
 
//     //Here be some testing code
// if(k==2){
//     VectorXd testvector(XY.size());
//     for (int i=0;i<testvector.size();i++) {
//         testvector[i]=0;
//     }
//     cout << "Now comparing delta E direct with delta E from gradient\n";
//     for (int i=0;i<testvector.size();i++) {
//             testvector[i]=1e-5;
//             cout << ESTRETCH << "\t" << Energynetwork(springlist,XY+testvector,g11,g12,g22)-ESTRETCH << "\t";
//             cout << gradE.dot(testvector) << "\t";
//             cout << EBEND << "\t" << Ebending(springpairs,springlist,XY+testvector,kappa,g11,g12,g22)-EBEND << "\t";
//             cout << (BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22)).dot(testvector) << endl;
//             testvector[i]=0;
//     }
// }
//     //end of testing code

EFile<<ESTRETCH<<"\t"<<EBEND<<"\t"<<ETOT<<"\t"<<0<<endl;
ESTRETCH=Energynetwork(springlist,XY,g11,g12,g22);
EBEND=Ebend(springpairs,springlist,XY,g11,g12,g22,kappa);  

//loop of the cg-method
int conjsteps=0;
do{
//     for(int j=0;j<XY.size();j++){ //write the XY-data to txt
//       XYfile<<XY(j)<<"\t";
//     } XYfile<<endl;
    conjsteps++;
    doConjStep(XY,s0,gradE,springlist,springpairs,kappa,conjsteps,g11,g12,g22);
    ESTRETCH=Energynetwork(springlist,XY,g11,g12,g22);
    EBEND=Ebend(springpairs,springlist,XY,g11,g12,g22,kappa);
    ETOT=ESTRETCH+EBEND;
    lenGrad=sqrt(gradE.dot(gradE));
    EFile<<ESTRETCH<<"\t"<<EBEND<<"\t"<<ETOT<<"\t"<<lenGrad<<endl; //Write the Energy to a txt-file.
}while(conjsteps<Nit && lenGrad>tolE);



for(int ii=0;ii<XY.size();ii++){
    shearcoordinates<<XY(ii)<<"\t";
}

//shearenergy<<ETOT<<"\t"<<angle<<endl;
shearenergy<<boxdx<<"\t"<<ETOT<<endl; //Write the deformation and the Energy to a txtfile
shearcoordinates<<endl;

if(k<100){
boxdx=boxdx+deltaboxdx;
}
if(k>100){
    boxdx=boxdx-deltaboxdx;
}
cout<<"This was a shearstep"<<endl;
}


XYfile.close();
EFile.close();
shearcoordinates.close();
shearenergy.close();


return 0;
}



