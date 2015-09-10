#include <iostream>
#include <ctime>
#include "random.h"
#include "makemikadonetwork.h"
#include "EnergyandGradients.h"
#include "minimizers.h"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <math.h>
#include <functional>
#include "importparam.h"
#include "BendingEnergy.h"
#include "BendingGrad.h"
#include "clusters.h"
#include <stdio.h>
#include <stdlib.h>
//#include "exportfiles.h"
#include "writefunctions.h"
#include "shearleader.h"

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
    int SEED;
    if(argc==0) SEED=0;
    if(argc>1){
        SEED=stoi(argv[1]);
        my_random::set_seed(SEED);
    }
//my_random::set_seed(0);
//System parameters  

    param Mikadoparameters;
    init(Mikadoparameters,"params.txt");

    vector<int> order;
    vector<stick> mikado(0);
    vector<stick> mikorig(0);
    vector<connected> Connection(0);
    vector<elonstick> ELONSTICK;
    vector<spring> springlist(0);
    vector<node> nodes(0);
    vector<node> singleNodes; 
    vector<vector<int>> springpairs(0);
    
    
    //a list which txt-files to produce
    ofstream mikadofile("mikado.txt"); 
    ofstream nodefile("nodes.txt");
    ofstream springfile("springs.txt");
    ofstream XYfile("conjpoints.txt");
    ofstream shearcoordinates("shearcoordinates.txt");
    ofstream shearenergy("shearenergy.txt");
    ofstream cluster("clusters.txt");
    ofstream clusterdata("clusterdata.txt", ios_base::app | ios_base::out);
    ofstream anglefile("angles.txt");
    char s[80];
    ofstream clusterdistribution(s);
    
    
    
    double lenGrad;
    int Nit=Mikadoparameters.Nit;  
    double tolGradE=Mikadoparameters.tolGradE;
    int NumberMikado=Mikadoparameters.NumberMikado;
    double LStick=Mikadoparameters.LStick; //Stick Length
    double k1=Mikadoparameters.k1;
    double k2=Mikadoparameters.k2;
    double kappa=Mikadoparameters.kappa;
    double rlenshort=Mikadoparameters.rlenshort;
    double rlenlong=Mikadoparameters.rlenlong;
    double stretchf=Mikadoparameters.stretchf;
    double deltaboxdx=Mikadoparameters.StepSize;
    int NumberStepsRight=Mikadoparameters.NumberStepsRight;
    int NumberStepsLeft=Mikadoparameters.NumberStepsLeft;
    
    int springone,springtwo,node1,node2,node3;

    cout<<"The number of rightssteps= "<<NumberStepsRight<<endl;
    cout<<"The number of leftsteps= "<<NumberStepsLeft<<endl;
 
    makeSticks(mikado,mikorig,NumberMikado,LStick);
    Write_Mikado_2txt(mikadofile,mikado);
    makeConnections(Connection,mikado,LStick); 
    sortELEMENTSperMIKADO(ELONSTICK,Connection);
    orderElonstick(order,ELONSTICK); 
    makeSpringsAndNodes(ELONSTICK,mikorig,springlist,nodes,rlenshort,rlenlong,k1,k2,stretchf);//Make the springs and Nodes. 
    Write_Springs_2txt(springfile,springlist);

    //make a table with sticks that are connected
    vector<int> NEWROW(2);
    vector<vector<int>> ConnectSticks;
        for(int i=0;i<Connection.size();i++){
            NEWROW[0]=Connection[i].first; 
            NEWROW[1]=Connection[i].second;
            ConnectSticks.push_back(NEWROW);
        }

    //Clusterv is a vector w. on its entries the clusters 
    vector<vector<int>> conmatr=connectivitymatrix(ConnectSticks,NumberMikado);
    vector<vector<int>> Clusterv=clusters(conmatr);
    //Make a cluster size distribution from the clusters
    vector<vector<int>> numberdistribution=Numberdistribution(Clusterv,NumberMikado);

    Write_Clusterdistribution_2txt(clusterdistribution,numberdistribution);
    Write_Clusterdistribution_2txt(cluster,Clusterv);
    Write_Clusterdata_2txt(clusterdata,NumberMikado,SEED,LStick,Clusterv);

    //Remove all double info
    for(std::size_t i=0;i<nodes.size();i++){
        if(nodes[i].number!=nodes[i+1].number){
            node unique=nodes[i];   
            singleNodes.push_back(unique);  
        }
    }

    //MAKE HERE THE PAIR OF SPRINGS
    makeSpringpairs(springpairs,springlist);
for(int i=0;i<springpairs.size();i++){
    springone=springpairs[i][0];
    springtwo=springpairs[i][1];
    node1=springlist[springone].one;
    node2=springlist[springone].two;
    node3=springlist[springtwo].two;
    
  Write_Angles_2txt(anglefile,node1,node2,node3);
    
}


VectorXd X(singleNodes.size()),Y(singleNodes.size());
VectorXd XY(2*X.size());
VectorXd gradE(XY.size());
VectorXd XYn(XY.size());
VectorXd XYcopy(XY.size());
VectorXd gradEn(gradE.size());
VectorXd s0(gradE.size());


    Write_Nodes_2txt(nodefile,singleNodes);

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

    double boxdx=0;
    double g11,g12,g22;
    double EBEND,ESTRETCH,ETOT;
    
    for(int k=0;k<(NumberStepsRight+NumberStepsLeft);k++){
         g11=1.0;
         g12=boxdx;
         g22=1.0+boxdx*boxdx;

         EBEND=EbendingC(springpairs,springlist,XY,kappa,g11,g12,g22);
         ESTRETCH=Energynetwork(springlist,XY,g11,g12,g22);
         ETOT=ESTRETCH+EBEND;
	 //Here comes the conjugate gradient
        gradE=HarmonicGradient(springlist,XY,g11,g12,g22)+BendingGrad(springpairs,springlist,XY,kappa,g11,g12,g22);
        s0=-gradE;

        //loop of the cg-method
        int conjsteps=0;
        do{
            conjsteps++;
            doConjStep(XY,s0,gradE,springlist,springpairs,kappa,conjsteps,g11,g12,g22);
            ESTRETCH=Energynetwork(springlist,XY,g11,g12,g22);
            EBEND=EbendingC(springpairs,springlist,XY,kappa,g11,g12,g22);  
            ETOT=ESTRETCH+EBEND;
            lenGrad=sqrt(gradE.dot(gradE));
        }while(conjsteps<Nit && lenGrad>tolGradE);

        Write_ShearCoordinates_2txt(shearcoordinates,XY);
	Write_ShearEnergy_2txt(shearenergy,boxdx,ETOT,ESTRETCH,EBEND,lenGrad,conjsteps);

        //Perform the sheardeformation for the next step.
        if(k<NumberStepsRight){
        boxdx=boxdx+deltaboxdx;
        }
        if(k>NumberStepsRight){ 
            boxdx=boxdx-deltaboxdx;
        }
        cout<<"\rThis was a shearstep     "<<k<<flush;
    }


    XYfile.close();
    shearcoordinates.close();
    shearenergy.close();
    cout<<endl;

    return 0;
}



