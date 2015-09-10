#include <iostream>
#include <ctime>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <math.h>
#include "makemikadonetwork.h"
#include "clusters.h"

using namespace std;
using namespace Eigen;




vector<vector<int>> connectivitymatrix(vector<spring> springs)
{
    vector<int> ones;
    vector<int> seconds;
    for(int i=0;i<springs.size();i++){
        //cout<<springs[i].one<<"  "<<springs[i].two<<endl;
        ones.push_back(springs[i].one);
        seconds.push_back(springs[i].two);
    }
    int numbersprings=ones.size()-1;
    vector<int> sorted=ones;
    
     //check how many nodes are there by sorting the lists one and two and taking
    //the max. element
    sort(sorted.begin(),sorted.end());
    int max=sorted[numbersprings];
    sorted.clear();
    sorted=seconds;
    sort(sorted.begin(),sorted.end());
    max<sorted[numbersprings]? max=sorted[numbersprings] :1 ;
    max++; //because start counting at 0
    
    vector<vector<int>> M(max,vector<int> (max));
    for(int i=0;i<max;i++){
        for(int j=0;j<max;j++){
            M[i][j]=0;
            if(i==j) M[i][j]=1;
        }
    }
    //cout<<"M  "<<M.size()<<"  "<<M[1].size()<<"  "<<max<<endl;

    int k,l;
    for(int i=0;i<ones.size();i++){
        k=ones[i];
        l=seconds[i];
        M[k][l]=1;
        M[l][k]=1;
    }
    
    return M;
    
};


vector<vector<int>> connectivitymatrix(vector<vector<int>> connectedSticks,int numbermikado)
{
    vector<int> ones;
    vector<int> seconds;
    for(int i=0; i<connectedSticks.size();i++){
    ones.push_back(connectedSticks[i][0]);
    seconds.push_back(connectedSticks[i][1]);
    }
    //Make an identitymatrix;
    vector<vector<int>> M(numbermikado,vector<int> (numbermikado));
    for(int i=0;i<numbermikado;i++){
        for(int j=0;j<numbermikado;j++){
            M[i][j]=0;
            if(i==j) M[i][j]=1;
        }
    }

    //Fill the rest.
    int k,l;
    for(int i=0;i<ones.size();i++){
        k=ones[i];
        l=seconds[i];
        M[k][l]=1;
        M[l][k]=1;
    }    
    
    return M;
    
}

void addelement(int element, vector<int> &v){
    //adds unique elements to a list. 
    int flag=0;
    int no_ele=0;
    if(v.size()==0){
        v.push_back(element);
        return;
    }
    
    for(int i=0;i<v.size();i++){
        if(v[i]==element){
            flag=1;
            break;
        }
    }
    if(flag==0) v.push_back(element);
};

vector< vector< int> > clusters(vector<vector<int> > M)
{
    vector< vector< int> > Clusters;
    vector<int> stack(0);
    vector<int> out(0);
    int seed;
    int numbernodes=M.size();
    //Here we create the connectivitymatrix.
    vector<vector<site>> m(numbernodes,vector<site> (numbernodes) );
    for(int i=0;i<numbernodes;i++){
        for(int j=0;j<numbernodes;j++){
            m[i][j].number=M[i][j];
            m[i][j].beenthere=0;
        }
    }
    for(int n=0;n<numbernodes;n++){
        if(m[n][n].beenthere==0){ //Loops over the diagonal, to check if the n^th node is already in a 
                                  //cluster if not, begin a new cluster
            m[n][n].beenthere=1;
            seed=n;
            stack.clear();
            stack.push_back(seed);
            out.clear();
            addelement(seed,out); //ads seed to out if seed wasn't already there.
            //out.push_back(seed);
            m[seed][seed].beenthere=1; 
            
            do{ //Loop over ONE cluster
                for(int k=0;k<numbernodes;k++){
                    if((m[seed][k].number==1 && m[seed][k].beenthere==0 && m[k][seed].beenthere==0) && seed!=k){

                        addelement(k,out);
                        stack.push_back(k);
                        m[seed][k].beenthere=1;
                        m[k][seed].beenthere=1;
                        m[k][k].beenthere=1;
                        //cout<<"check2"<<endl;
                        // Not really nescecary for counting in the cluster, but it might be handy to check that node k is
                        // already part of a cluster, so we don't have to take this one into account for looking for more 
                        //clusters.

                        seed=k;
                        break;
                    }else if(k==numbernodes-1){
                        stack.pop_back();
                        //cout<<"check"<<endl;
                        if(stack.size()==0) break;
                        seed=stack[stack.size()-1];
                    }
                }
            }while(stack.size()>0);
            sort(out.begin(),out.end());
        Clusters.push_back(out);
        }
    }

return Clusters;    
}

vector<vector<int>> Numberdistribution(vector<vector<int>> Clusters,int NSticks)
{
    vector<vector<int>> numberdistribution(NSticks,vector<int> (2));
    //initiate the array with 'bins'
    for(int i=0;i<NSticks;i++){
        numberdistribution[i][0]=i+1;
        numberdistribution[i][1]=0;
    }
    
    int size;
    for(int i=0;i<Clusters.size();i++){
        size=Clusters[i].size(); //the size of the i-th cluster
        numberdistribution[size-1][1]++;
    }
    
    return numberdistribution;
    
    
};
