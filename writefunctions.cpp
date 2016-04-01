#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"
#include <fstream>
#include "makemikadonetwork.h"
#include "clusters.h"
using namespace Eigen;
using namespace std;

void Write_Mikado_2txt(ofstream &mikadofile,vector<stick> &mikado){
    for(int i=0;i<mikado.size();i++){
        mikadofile<<mikado[i].nr<<"\t"<<mikado[i].x<<"\t"<<mikado[i].y<<"\t"<<mikado[i].th<<"\t"<<mikado[i].wlr<<"\t"<<
        mikado[i].wud<<endl;
    } 
    mikadofile.close();
};


void Write_Springs_2txt(ofstream &springfile,vector<spring> &springlist){
    for(int i=0;i<springlist.size();i++){
        springfile<<springlist[i].one<<"\t"
                  <<springlist[i].two<<"\t"
                <<springlist[i].wlr<<"\t"
                <<springlist[i].wud<<"\t"
                <<springlist[i].rlen<<"\t"
                <<springlist[i].k<<"\t"
                <<springlist[i].sticki<<endl;       
    }
    springfile.close();
};


void Write_Clusterdistribution_2txt(ofstream &clusterdistribution,vector<vector<int>> &numberdistribution){
    for(int i=0;i<numberdistribution.size();i++){
            clusterdistribution<<numberdistribution[i][0]<<"    "<<numberdistribution[i][1]<<endl;
        }
};


void Write_Clusterdata_2txt(ofstream &clusterdata,int SEED,int NumberMikado,double LStick,vector<vector<int>> &C){
    clusterdata<<SEED<<"    "<<NumberMikado<<"      "<<LStick<<"    "<<C.size()<<endl;
    clusterdata.close();   
};

void Write_Nodes_2txt(ofstream &nodefile,vector<node> &singleNodes){
    for(int i=0;i<singleNodes.size();i++){
        nodefile<<singleNodes[i].number<<"\t"<<singleNodes[i].x<<"\t"<<singleNodes[i].y<<endl;
    } 
    nodefile.close();
};

void Write_ShearCoordinates_2txt(ofstream &shearcoordinates,VectorXd &XY){
    for(int ii=0;ii<XY.size();ii++){
	shearcoordinates<<XY(ii)<<"\t";
         }
         shearcoordinates<<endl; 

};

void Write_ShearEnergy_2txt(ofstream &shearenergy,double boxdx,double ETOT,double ESTRETCH, double EBEND, double lenGrad, int conjsteps){
  shearenergy<<boxdx<<"\t"<<ETOT<<"\t"<<ESTRETCH<<"\t"<<EBEND<<"\t"<<lenGrad<<"\t"<<conjsteps<<endl;
  
};

void Write_Angles_2txt(ofstream &anglefile,int node1,int node2,int node3){
      anglefile<<node1<<"\t"<<node2<<"\t"<<node3<<endl;  
};

void Write_IndividualClusters_2txt(ofstream &cluster,vector<vector<int>> Clusterv){
    for(int i=0;i<Clusterv.size();i++){
        for(int j=0;j<Clusterv[i].size();j++){
            cluster<<Clusterv[i][j]<<",";
        }
        cluster<<endl;
    }
    cluster.close();
}


void Write_Stresstensor(ofstream &stresstens, double strain, stresstensor sigma){
    stresstens<<strain<<"\t"<<sigma.sxx<<"\t"<<sigma.sxy<<"\t"<<sigma.syy<<"\t"<<0.5*(sigma.sxx+sigma.syy)<<endl;
}
