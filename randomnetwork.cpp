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
#include "makeSpringNetworks.h"
#include "structs.h"


using namespace std;
using namespace Eigen;

VectorXd initiateRandomNetwork(vector<spring> &springlist,
                                vector<vector<int>> &springpairs,
                                vector<stick> &mikado,
                                vector<stick> &mikorig,
                                vector<elonstick> &ELONSTICK,
                                vector<connected> &Connection,
                                vector<node> &nodes,
                                vector<node> &singleNodes,
                                vector<vector<int>> &conmatr,
                                vector<vector<int>> &Clusterv,
                                vector<vector<int>> &numberdistribution,
                                int NumberMikado,
                                vector<spring> &background,
                                VectorXd &XYb,
                                int SEED,
                                double LStick,
                                double rlenshort,
                                double rlenlong,
                                double k1,
                                double k2,
                                double stretchf,
                                ofstream &springfile,
                                ofstream &anglefile,
                                ofstream &mikadofile,
                                ofstream &clusterdistribution,
                                ofstream &cluster,
                                ofstream &clusterdata,
                                ofstream &nodefile                                
                               )
                               

{
    vector<int> order;
    makeSticks(mikado,mikorig,NumberMikado,LStick);
    Write_Mikado_2txt(mikadofile,mikado);
    makeConnections(Connection,mikado,LStick,background,XYb); 
    sortELEMENTSperMIKADO(ELONSTICK,Connection);
    orderElonstick(order,ELONSTICK); 
    makeSpringsAndNodes(ELONSTICK,mikorig,springlist,nodes,rlenshort,rlenlong,k1,k2,stretchf);//Make the springs and Nodes. 
    //Write_Springs_2txt(springfile,springlist);

    //make a table with sticks that are connected
    vector<int> NEWROW(2);
    vector<vector<int>> ConnectSticks;
        for(int i=0;i<Connection.size();i++){
            NEWROW[0]=Connection[i].first; 
            NEWROW[1]=Connection[i].second;
            ConnectSticks.push_back(NEWROW);
        }

    //Clusterv is a vector w. on its entries the clusters 
    conmatr=connectivitymatrix(ConnectSticks,NumberMikado);
    Clusterv=clusters(conmatr);
    //Make a cluster size distribution from the clusters
    numberdistribution=Numberdistribution(Clusterv,NumberMikado);

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

    
    VectorXd X(singleNodes.size()),Y(singleNodes.size());
    VectorXd XY(2*X.size());
    VectorXd XYn(XY.size());
    VectorXd XYcopy(XY.size());
   
   int springone,springtwo;
   int node1,node2,node3;
    
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

return XY;
}