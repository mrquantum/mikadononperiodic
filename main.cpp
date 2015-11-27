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
#include "writefunctions.h"
#include "shearleader.h"
#include "randomnetwork.h"
#include "makeSpringNetworks.h"
#include "structs.h"
#include "cutsprings.h"
#include "connectivity.h"
#include "cgmethod.h"

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
    vector<vector<int>> conmatr;
    vector<vector<int>> Clusterv;
    vector<vector<int>> numberdistribution;

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

    //double lenGrad;
    int Nit=Mikadoparameters.Nit;
    double tolGradE=Mikadoparameters.tolGradE;
    int NumberMikado=Mikadoparameters.NumberMikado;
    double LStick=Mikadoparameters.LStick; //Stick Length
    double k1=Mikadoparameters.k1;
    double k2=Mikadoparameters.k2;
    int bendingOn=Mikadoparameters.bendingOn;
    double kappa=Mikadoparameters.kappa;
    double rlenshort=Mikadoparameters.rlenshort;
    double rlenlong=Mikadoparameters.rlenlong;
    double stretchf=Mikadoparameters.stretchf;
    double deltaboxdx=Mikadoparameters.StepSize;
    int NumberStepsRight=Mikadoparameters.NumberStepsRight;
    int NumberStepsLeft=Mikadoparameters.NumberStepsLeft;
    int backGroundOn=Mikadoparameters.backGroundOn;
    int backGroundType=Mikadoparameters.backGroundType;
    int mode=Mikadoparameters.mode;
    double Z_aim=Mikadoparameters.Z_aim;
    
//     cout<<"The number of rightssteps= "<<NumberStepsRight<<endl;
//     cout<<"The number of leftsteps= "<<NumberStepsLeft<<endl;

    
        vector<spring> background(0);
        VectorXd XYb;
        VectorXd XY;
        VectorXd gradE(XY.size());
        VectorXd gradEn(gradE.size());
        VectorXd s0(gradE.size());
        double Z; //Connectivity
        int Numberf=18;
        //initiate the backgroundnetwork
        if(backGroundOn==1){
            XYb=makeSquareNetwork(Numberf,background);
            cout<<"SQUARE"<<endl;
        } else if(backGroundOn==2){             
            XYb=triangularnetwork(Numberf,background);
            cout<<"TRIANGLE"<<endl;
        }


    if(mode==0){
         if(NumberMikado>0){
                    XY=initiateRandomNetwork(springlist,springpairs,mikado,mikorig,ELONSTICK,Connection,nodes,
                                            singleNodes,conmatr,Clusterv,numberdistribution,NumberMikado,background,XYb,
                                            SEED,LStick,rlenshort,rlenlong,k1,k2,stretchf,
                                            springfile,anglefile,mikadofile,clusterdistribution,cluster,
                                            clusterdata,nodefile);
                } else{
                    XY=XYb;
                    springlist=background;
                }
                Z=Connectivity(springlist);
    }
    
    
    ofstream zfile("zfile.txt",ios_base::app | ios_base::out);
    if(mode==1){
        int number=50;
        //initiate the mikadonetwork
        do{
            springlist.clear();
            springpairs.clear();
            mikado.clear();
            mikorig.clear();
            ELONSTICK.clear();
            Connection.clear();
            nodes.clear();
            singleNodes.clear();
            conmatr.clear();
            NumberMikado=number;
                if(NumberMikado>0){
                    XY=initiateRandomNetwork(springlist,springpairs,mikado,mikorig,ELONSTICK,Connection,nodes,
                                            singleNodes,conmatr,Clusterv,numberdistribution,NumberMikado,background,XYb,
                                            SEED,LStick,rlenshort,rlenlong,k1,k2,stretchf,
                                            springfile,anglefile,mikadofile,clusterdistribution,cluster,
                                            clusterdata,nodefile);
                } else{
                    XY=XYb;
                    springlist=background;
                }
                
                Z=Connectivity(springlist);
                zfile<<NumberMikado<<"\t"<<Z<<endl;
                number++;
        }while(Z<Z_aim);
    }
    
        networkinfo info;
        info.g11=1.0;
        info.g12=0.0;
        info.g22=1.0;
        info.springlist=springlist;
        info.size=XY.size();
  
        VectorXd Xi(XY.size());
        HarmonicGradientn(XY.data(),Xi.data(),info);
        VectorXd Xiold=HarmonicGradient(springlist,XY,1.0,0.0,1.0);
        
        //test for the new CG method
        int iter;
        double fret;
        cout<<EnergyNetworkn(XY.data(),info)<<" FUCK NOH"<<endl;
        frprmn(XY.data(),XY.size(),1.0e-10,&iter,&fret,EnergyNetworkn,HarmonicGradientn,info);
        cout<<EnergyNetworkn(XY.data(),info)<<" FUCK YEAH"<<endl;
        //end test
        cout<<iter<<"       steps"<<endl;

    Write_Springs_2txt(springfile,springlist);

    //Shearing
    vector<vector<int>> springp(0);
    shearsteps(deltaboxdx,NumberStepsRight,NumberStepsLeft,springlist,
            springp,XY,bendingOn,kappa,Nit,tolGradE,shearcoordinates,shearenergy);
    XYfile.close();
    shearcoordinates.close();
    shearenergy.close();
    zfile.close();
    cout<<endl;
    return 0;
}



