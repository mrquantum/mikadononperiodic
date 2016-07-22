#include <iostream>
#include <ctime>
#include "random.h"
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
#include "combineElementsOnMikado.h"

using namespace Eigen;
using namespace std;
const double pi=4.0*atan(1.0);

//This function overloads the < operator, so that it can
//compare nodes
bool operator<(const node& first,const node& second){
   if(first.number>=second.number){
   return false;}
   else{
     return true;
   }
}

bool operator<(const elonstick& first, const elonstick& second){
  if(first.sticki>=second.sticki){
    return false;}
    else{
      return true;
    }
}

bool operator<(const stick &first, const stick &second)
{
  if(first.nr>=second.nr){
    return false;}
    else{ return true;
    }
   }



//Here we make the initial mikadonetwork
std::vector<stick> make_sticks(int N,double sticklen)
{
  my_random::get_gre(2);
std::vector<stick> m(N);
//stick S;
for(int i=0;i<N;i++){
  m[i].x=randf();
  m[i].y=randf();
  m[i].th=2*pi*randf();
  m[i].nr=i;
  m[i].wlr=0;
  m[i].wud=0;
  m[i].length=sticklen;
  }
return m;
 }
 
//Check for passes through the left and right wall 
vector<stick> make_ghost_lr( const vector<stick> &m, double LStick, int NumberMikado)
{
  std::vector<stick> GhostLR;
  for (int i=0;i<NumberMikado;i++){
    if(m[i].th!=pi/2 || m[i].th!=3*pi/2){ // Check here for the left wall 
      stick ghostRowLR;
      Matrix2d A; //The matrix to solve with
      Vector2d b,b2; // A*st=b
      Vector2d st,st2;
 
      A<<-cos(m[i].th),0,-sin(m[i].th),1; 
      b<<m[i].x,m[i].y;
      b2<<m[i].x-1,m[i].y;
      st=A.lu().solve(b);
      st2=A.lu().solve(b2);
      if((st(0)>0 && st(0)<LStick) && (st(1)>0&&st(1)<1)){ //If mikado passes throug wall make ghost mikado 
        ghostRowLR=m[i];
        ghostRowLR.x=ghostRowLR.x+1;
        ghostRowLR.wlr=ghostRowLR.wlr-1;
        GhostLR.push_back(ghostRowLR); 
      }
      if((st2(0)>0&&st2(0)<LStick)&&(st2(1)>0&&st2(1)<1)){ //Now do the same for the upper wall
        ghostRowLR=m[i];
        ghostRowLR.x=ghostRowLR.x-1;
        ghostRowLR.wlr=ghostRowLR.wlr+1;
        GhostLR.push_back(ghostRowLR);
      }
    }
  }
  return GhostLR;
}

//Check for passes through the down and up wall;
vector<stick> make_ghost_ud(const vector<stick> &m, double LStick, int NumberMikado){
  vector<stick> GhostUD;
  for(int i=0;i<NumberMikado;i++){
    if(m[i].th!=0||m[i].th!=pi){
      stick ghostRowUD;
      Matrix2d A;
      Vector2d b,b2;
      Vector2d st, st2;
      
      A<<-cos(m[i].th),1,-sin(m[i].th),0;
      b<<m[i].x,m[i].y;
      b2<<m[i].x,m[i].y-1;
      st=A.lu().solve(b);
      st2=A.lu().solve(b2);
      if((st(0)>0&&st(0)<LStick)&&(st(1)>0&&st(1)<1)){
	ghostRowUD=m[i];
	ghostRowUD.y=ghostRowUD.y+1;
	ghostRowUD.wud=ghostRowUD.wud-1;
	GhostUD.push_back(ghostRowUD);
	}
      if((st2(0)>0&&st2(0)<LStick)&&(st2(1)>0&&st2(1)<1)){
	ghostRowUD=m[i];
	ghostRowUD.y=ghostRowUD.y-1;
	ghostRowUD.wud=ghostRowUD.wud+1;
	GhostUD.push_back(ghostRowUD);
	}
     }
    
  }
  return GhostUD;
}

//Checks whether the mikado's are connected
void make_connections(vector<connected> &Connection,
                      const vector<stick> &m, 
                      double LStick,
                      const vector<spring> &background,
                      const VectorXd &XYb)
{
  Matrix2d A; //The matrix to solve with
  Vector2d b,st; // A*st=b
  int nr=0;
  connected xtrarow, xtrarow2;

  //Loop over all sticks to find coordinates
    for(int i=0;i<m.size()-1;i++){
        for(int j=i+1;j<m.size();j++){
            if(m[i].th!=m[j].th || m[i].nr!=m[j].nr){
                A<<-cos(m[i].th),cos(m[j].th),-sin(m[i].th),sin(m[j].th);
                b<<m[i].x-m[j].x,m[i].y-m[j].y;
                st=A.lu().solve(b);
                if ((st(0)>0.0 && st(0)<LStick)&&(st(1)>0.0 && st(1)<LStick)){
                    xtrarow.first=m[i].nr;	//[stick i stick j sij sji] 
                    xtrarow.second=m[j].nr;
                    xtrarow.s1=st(0);
                    xtrarow.s2=st(1);
                    xtrarow.nrCon=nr;
                    xtrarow.recur=0;
                    xtrarow.type=0;
                    xtrarow.backgroundspring[0]=-1;
                    xtrarow.backgroundspring[1]=-1;
                    
                    xtrarow2.first=m[j].nr;	//[stick j sticki sji sij]
                    xtrarow2.second=m[i].nr;
                    xtrarow2.s1=st(1);
                    xtrarow2.s2=st(0);
                    xtrarow2.nrCon=nr;
                    xtrarow2.recur=0;
                    xtrarow2.type=0;
                    xtrarow2.backgroundspring[0]=-1;
                    xtrarow2.backgroundspring[1]=-1;
                    Connection.push_back(xtrarow);
                    Connection.push_back(xtrarow2);
                    nr++;
                }
            }
        }
    }

  //Now loop over the background springs
    int one;
    int two;
    int number=XYb.size()/2;
    double xsb,xse,ysb,yse,xmik,ymik;
    double thspr, thmik;
    double lenspring; //the coordinates + parameters of the spring 
    double s,t;
    if(background.size()>0){//check wheather this block needs to be executed
        for(int i=0; i<m.size();i++){
            xmik=m[i].x;
            ymik=m[i].y;
            thmik=m[i].th;

            for(int j=0;j<background.size();j++){
            //calculate the intersections similar to the previous
                one=background[j].one;
                two=background[j].two;
                xsb=XYb(one);
                ysb=XYb(one+number);
                xse=XYb(two)+background[j].wlr;
                yse=XYb(two+number)+background[j].wud;
                thspr=atan2((yse-ysb),(xse-xsb));
                lenspring=sqrt(pow((yse-ysb),2)+pow((xse-xsb),2));
        
                if(fabs(thspr-thmik)>1e-10 && fabs(fabs(thspr-thmik)-pi)>1e-10){
                    A<<cos(thmik),-cos(thspr),sin(thmik),-sin(thspr);
                    b<<xsb-xmik,ysb-ymik;
                    st=A.lu().solve(b);
                    s=st(0);
                    t=st(1);
            
                    if(s>0.0 && t>0.0 && s<LStick && t<lenspring){
                        //then we found a node!   
                        xtrarow.first=m[i].nr;
                        xtrarow.second=-1; //the signature of a collision w spring
                        xtrarow.s1=s;
                        xtrarow.s2=t;
                        xtrarow.nrCon=nr;
                        xtrarow.type=1;
                        xtrarow.backgroundspring[0]=one;
                        xtrarow.backgroundspring[1]=two;
                        xtrarow.recur=0;
                        Connection.push_back(xtrarow);
                        nr++;
                    }
                }
            }
        }
    }
}


void sortELEMENTSperMIKADO(vector<elonstick> &ELONSTICK,vector<connected> &Connection)
{
  combineElementsOnMikado(Connection,ELONSTICK);
  
// now sort extrarow on descending order per stick;
for(int j=0; j<ELONSTICK.size();j++){
    vector<double> distances=ELONSTICK[j].S;
    vector<int> numbers=ELONSTICK[j].nr;
    vector<int> type=ELONSTICK[j].type;
    vector<array<int,2>> backgroundspring=ELONSTICK[j].backgroundspringn;
// now sort extrarow on descending order per stick;
    
    int swapped=0;
      do{
        int k=0;
        swapped=0; //this is the control parameter, checks 1 if elements are swapped
            if(distances.size()>1){
                for(int i=0;i<distances.size()-1;i++){ //loop through the list thill the end-k-1 th element;
                    if(distances[i]>distances[i+1]){ //checks if neighbours are in right order, if not then swap en change swap parameter
                        double d1=distances[i];
                        double d2=distances[i+1];
                        distances[i]=d2;
                        distances[i+1]=d1;
                        int n1=numbers[i]; 
                        int n2=numbers[i+1];
                        numbers[i]=n2; 
                        numbers[i+1]=n1;
                        int t1=type[i];
                        int t2=type[i+1];
                        type[i]=t2;
                        type[i+1]=t1;
                        int backspring1 [2]={backgroundspring[i][0],backgroundspring[i][1]};
                        int backspring2 [2]={backgroundspring[i+1][0],backgroundspring[i+1][1]};
                        backgroundspring[i][0]=backspring2[0];
                        backgroundspring[i][1]=backspring2[1];
                        backgroundspring[i+1][0]=backspring1[0];
                        backgroundspring[i+1][1]=backspring1[1];
                        
                        swapped=1;
                        k++;
                    }
                }
            }
    } while(swapped==1);
        ELONSTICK[j].S=distances; //Put the new data back into the original vectors
        ELONSTICK[j].nr=numbers;
	ELONSTICK[j].type=type;
        ELONSTICK[j].backgroundspringn=backgroundspring;
    }
std::sort(ELONSTICK.begin(),ELONSTICK.end());
}

void orderElonstick(vector<int> &order,vector<elonstick> &ELONSTICK)
{
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
}


//With the points per stick, we can now make nodes and springs.
void makeSpringsAndNodes(const vector<elonstick> &ELONSTICK,const vector<stick> &mikorig, vector<spring> &springlist,
vector<node> &nodes,
double rlenshort, double rlenlong,double k1,double k2,double stretchf,vector<spring> &background,VectorXd &XYb)
{

    int background_size=XYb.size()/2;
    //first make the background networks into 'nodes'
    VectorXd Xb(background_size), Yb(background_size);
    node backgroundnode;
    for(int i=0;i<background_size;i++){
        Xb(i)=XYb(i);
        Yb(i)=XYb(i+background_size);
        backgroundnode.number=i;
        backgroundnode.x=Xb(i);
        backgroundnode.y=Yb(i);
        nodes.push_back(backgroundnode);
    }
    springlist=background;
    //we have now set up the background spring network in terms of the struct node and we the springs
    //are loaded into the springlist as the first N springs. 
    //Now we are going to append the new springs to the network.

    for(int i=0;i<ELONSTICK.size();i++){
        int sticknr=ELONSTICK[i].sticki;
        vector<int> nodesonsticki=ELONSTICK[i].nr;
        vector<double> posonsticki=ELONSTICK[i].S;
        stick CURRENTSTICK=mikorig[sticknr];
        int node1, node2;
        double lenspring;
        double springconstant;
        if(nodesonsticki.size()>1){
            for(int j=0;j<nodesonsticki.size()-1;j++){
                spring newspring;
                node1=nodesonsticki[j]+background_size;
                node2=nodesonsticki[j+1]+background_size;
                lenspring=posonsticki[j+1]-posonsticki[j];
                springconstant=lenspring*k1/CURRENTSTICK.length;
//                 springconstant=0.5;
                double x1, x2, y1, y2;
                x1=CURRENTSTICK.x+posonsticki[j]*cos(CURRENTSTICK.th); //calculate the position of the node
                x2=CURRENTSTICK.x+posonsticki[j+1]*cos(CURRENTSTICK.th);//and the position of the adjacent one
                y1=CURRENTSTICK.y+posonsticki[j]*sin(CURRENTSTICK.th);
                y2=CURRENTSTICK.y+posonsticki[j+1]*sin(CURRENTSTICK.th);
                newspring=makespring(node1,node2,x1,x2,y1,y2,sticknr,springconstant,stretchf);
                node nodetemp1, nodetemp2;
                nodetemp1.number=newspring.one;
                nodetemp1.x=x1-floor(x1);
                nodetemp1.y=y1-floor(y1);
                
                nodetemp2.number=newspring.two;
                nodetemp2.x=x2-floor(x2);
                nodetemp2.y=y2-floor(y2);
                nodes.push_back(nodetemp1);
                nodes.push_back(nodetemp2);
                springlist.push_back(newspring);

                int backnode1,backnode2,midpoint;
                int wlr,wud;
                double xbn1,xbn2,ybn1,ybn2,xm,ym;
                double xm1,xm2,ym1,ym2;
                
                //first and last element
                if(j==0 && ELONSTICK[i].type[0]==1){
                    makeanddeletebondsonbackground(springlist,ELONSTICK,CURRENTSTICK,posonsticki,Xb,Yb,background_size,i,0);
                }
                //check for the other elements
                if(ELONSTICK[i].type[j+1]==1){
                    makeanddeletebondsonbackground(springlist,ELONSTICK,CURRENTSTICK,posonsticki,Xb,Yb,background_size,i,j+1);
                }
            }
        }
    }
    //sort the nodes and if the background absent remove double nodes.
    //in the other case this happens in makeanddeletebondsonbackground
    std::sort(nodes.begin(),nodes.end());
    if(background.size()==0) {
        ordernodes(nodes,springlist);
    }
 
    int mi=0;
    int ma=0;
    for(int i=0;i<springlist.size();i++){
        springlist[i].one>ma ? ma=springlist[i].one : ma;
        springlist[i].two>ma ? ma=springlist[i].two : ma;

    }
}

// void makeSpringsAndNodes(const vector<elonstick> &ELONSTICK,const vector<stick> &mikorig, vector<spring> &springlist,
// vector<node> &nodes,
// double rlenshort, double rlenlong,double k1,double k2,double stretchf,vector<spring> &background,VectorXd &XYb)
// {
// 
//     int background_size=XYb.size()/2;
//     //first make the background networks into 'nodes'
//     VectorXd Xb(background_size), Yb(background_size);
//     node backgroundnode;
//     for(int i=0;i<background_size;i++){
//         Xb(i)=XYb(i);
//         Yb(i)=XYb(i+background_size);
//         backgroundnode.number=i;
//         backgroundnode.x=Xb(i);
//         backgroundnode.y=Yb(i);
//         nodes.push_back(backgroundnode);
//     }
//     springlist=background;
//     //we have now set up the background spring network in terms of the struct node and we the springs
//     //are loaded into the springlist as the first N springs. 
//     //Now we are going to append the new springs to the network.
// 
//     for(int i=0;i<ELONSTICK.size();i++){
//         int sticknr=ELONSTICK[i].sticki;
//         vector<int> nodesonsticki=ELONSTICK[i].nr;
//         vector<double> posonsticki=ELONSTICK[i].S;
//         stick CURRENTSTICK=mikorig[sticknr];
//         int node1, node2;
//         for(int j=0;j<nodesonsticki.size()-1;j++){
//             spring newspring;
//             node1=nodesonsticki[j]+background_size;
//             node2=nodesonsticki[j+1]+background_size;
//  
//             double x1, x2, y1, y2;
//             x1=CURRENTSTICK.x+posonsticki[j]*cos(CURRENTSTICK.th); //calculate the position of the node
//             x2=CURRENTSTICK.x+posonsticki[j+1]*cos(CURRENTSTICK.th);//and the position of the adjacent one
//             y1=CURRENTSTICK.y+posonsticki[j]*sin(CURRENTSTICK.th);
//             y2=CURRENTSTICK.y+posonsticki[j+1]*sin(CURRENTSTICK.th);
//             newspring=makespring(node1,node2,x1,x2,y1,y2,sticknr,1.0,stretchf);
//             node nodetemp1, nodetemp2;
//             nodetemp1.number=newspring.one;
//             nodetemp1.x=x1-floor(x1);
//             nodetemp1.y=y1-floor(y1);
//             nodetemp2.number=newspring.two;
//             nodetemp2.x=x2-floor(x2);
//             nodetemp2.y=y2-floor(y2);
//             nodes.push_back(nodetemp1);
//             nodes.push_back(nodetemp2);
//             springlist.push_back(newspring);
//             
//             int backnode1,backnode2,midpoint;
//             int wlr,wud;
//             double xbn1,xbn2,ybn1,ybn2,xm,ym;
//             double xm1,xm2,ym1,ym2;
//             
//             //check for the first element on the stick.
//             if(j==0 && ELONSTICK[i].type[0]==1){
//                 //spring to be removed has nodes
//                 backnode1=ELONSTICK[i].backgroundspringn[0][0];
//                 backnode2=ELONSTICK[i].backgroundspringn[0][1];
//                 for(int m=0;m<springlist.size();m++){
//                     if((springlist[m].one==backnode1 && springlist[m].two==backnode2) ||
//                         (springlist[m].one==backnode2 && springlist[m].two==backnode1)){
//                         //erase one spring.
//                         wud=springlist[m].wud;
//                         wlr=springlist[m].wlr;
//                         springlist.erase(springlist.begin()+m);   
//                         //replace it by two new ones. One between backnode1 and ELONSTICK[i].nr[j]+background_size
//                         //and ELONSTICK[i].nr[j]+backgroundsize and backgnode2
//                         xbn1=Xb(backnode1);
//                         ybn1=Yb(backnode1);
//                         xbn2=Xb(backnode2);
//                         ybn2=Yb(backnode2);
//                         xm=CURRENTSTICK.x+posonsticki[0]*cos(CURRENTSTICK.th);
//                         ym=CURRENTSTICK.y+posonsticki[0]*sin(CURRENTSTICK.th);
//                         midpoint=ELONSTICK[i].nr[0]+background_size;
//                         //place the midpoint coordinate in the box
//                         if(xm>1){
//                           xm=xm-1;
//                         }
//                         if(xm<0){
//                           xm=xm+1;
//                         }
//                         if(ym>1){
//                           ym=ym-1;
//                         }
//                         if(ym<0){
//                            ym=ym+1;
//                         }
//                         //Check how the original coordinate should be 'fealt' by the midpoint coordinate
//                         if(wlr==1) {
//                             if (fabs(xbn1-xm)<fabs(xbn2-xm)){
//                                 xbn2+=1;
//                             } else {
//                                 xbn1-=1;
//                             }
//                         }
//                         
//                         if(wlr==-1){
//                             if(fabs(xbn1-xm)<fabs(xbn2-xm)){
//                                 xbn2-=1.0;
//                             } else{
//                                 xbn1+=1.0; 
//                             }
//                         }
//                         
//                         if(wud==1) {
//                             if (fabs(ybn1-ym)<fabs(ybn2-ym)){
//                                 ybn2+=1;
//                             } else {
//                                 ybn1-=1;
//                             }
//                         }
//                         if(wud==-1){
//                             if(fabs(ybn1-ym)<fabs(ybn2-ym)){
//                                 ybn2-=1.0;
//                             } else{
//                                 ybn1+=1.0; 
//                             }
//                         }
//                         newspring=makespring(backnode1,midpoint,xbn1,xm,ybn1,ym,-1,1.0,1.0);
//                         springlist.push_back(newspring);
//                         newspring=makespring(midpoint,backnode2,xm,xbn2,ym,ybn2,-1,1.0,1.0);
//                         springlist.push_back(newspring);
//                         break;
//                     }
//                 }
//             }
// 
//             //check for the other elements
//             if(ELONSTICK[i].type[j+1]==1){
//                 makeanddeletebondsonbackground(springlist,ELONSTICK,CURRENTSTICK,posonsticki,Xb,Yb,background_size,i,j+1);
//             }
//             
//             //THIS SHOULD BE ONE FUNCTION
//             if(ELONSTICK[i].type[j+1]==1){
//                 //spring to be removed has nodes
//                 backnode1=ELONSTICK[i].backgroundspringn[j+1][0];
//                 backnode2=ELONSTICK[i].backgroundspringn[j+1][1];
//                 for(int m=0;m<springlist.size();m++){
//                     if((springlist[m].one==backnode1 && springlist[m].two==backnode2) ||
//                         (springlist[m].one==backnode2 && springlist[m].two==backnode1)){
//                         //erase one spring.
//                         wud=springlist[m].wud;
//                         wlr=springlist[m].wlr;
//                         springlist.erase(springlist.begin()+m);   
//                         //replace it by two new ones. One between backnode1 and ELONSTICK[i].nr[j]+background_size
//                         //and ELONSTICK[i].nr[j]+backgroundsize and backgnode2
//                         xbn1=Xb(backnode1);
//                         ybn1=Yb(backnode1);
//                         xbn2=Xb(backnode2);
//                         ybn2=Yb(backnode2);
//                         xm=CURRENTSTICK.x+posonsticki[j+1]*cos(CURRENTSTICK.th);
//                         ym=CURRENTSTICK.y+posonsticki[j+1]*sin(CURRENTSTICK.th);
//                         midpoint=ELONSTICK[i].nr[j+1]+background_size;
//                         
//                         //place the midpoint coordinate in the box
//                         if(xm>1){
//                           xm=xm-1;
//                         }
//                         if(xm<0){
//                           xm=xm+1;
//                         }
//                         if(ym>1){
//                           ym=ym-1;
//                         }
//                         if(ym<0){
//                            ym=ym+1;
//                         }
// 
//                         //Check how the original coordinate should be 'fealt' by the midpoint coordinate
//                         if(wlr==1) {
//                             if (fabs(xbn1-xm)<fabs(xbn2-xm)){
//                                 xbn2+=1;
//                             } else {
//                                 xbn1-=1;
//                             }
//                         }
//                         
//                         if(wlr==-1){
//                             if(fabs(xbn1-xm)<fabs(xbn2-xm)){
//                                 xbn2-=1.0;
//                             } else{
//                                 xbn1+=1.0; 
//                             }
//                         }
//                         
//                         if(wud==1) {
//                             if (fabs(ybn1-ym)<fabs(ybn2-ym)){
//                                 ybn2+=1;
//                             } else {
//                                 ybn1-=1;
//                             }
//                         }
//                         if(wud==-1){
//                             if(fabs(ybn1-ym)<fabs(ybn2-ym)){
//                                 ybn2-=1.0;
//                             } else{
//                                 ybn1+=1.0; 
//                             }
//                         }
//                         cout<<ybn1<<"   "<<ym<<"        "<<ybn2<<endl;
//                         newspring=makespring(backnode1,midpoint,xbn1,xm,ybn1,ym,-1,1.0,1.0);
//                         springlist.push_back(newspring);
//                         newspring=makespring(midpoint,backnode2,xm,xbn2,ym,ybn2,-1,1.0,1.0);
//                         springlist.push_back(newspring);
//                         break;
//                     }
//                 }
//             }
//         }
//     }
//     std::sort(nodes.begin(),nodes.end());
// }
    
    
    

double inbox(double x,double boxsize){
  if(x<0){
    x=x+boxsize;
  }
  if(x>boxsize){
  x=x-boxsize;
  }
return x;
}

void makeSpringpairs(vector<vector<int>> &springpairs,const vector<spring> &springlist)
{
//springpairs is an vector, coding a triplet in each row: 
//springpairs[i][1]=ith triplet first spring
//springpairs[i][2] ith triplet second spring. 
//springpairs[i][3] ith triplet mikado nr.
//the labels for the springs are the indexnumbers of coding the
//springs in the springslist.    
    
    //first sive out the backgroundsprings
    vector<spring> newspringlist;
    for(int i=0;i<springlist.size();i++){
        if(springlist[i].sticki!=-1){
            newspringlist.push_back(springlist[i]);
        }
        
    }

    for(std::size_t i=0;i<springlist.size()-1;i++){
        if(springlist.size()>1){
            for(int j=i+1;j<springlist.size();j++)
            {
            //if(springlist[i].sticki==springlist[i+1].sticki && springlist[i].sticki!=-1){
                if(springlist[i].sticki==springlist[j].sticki && springlist[i].sticki!=-1){
                    vector<int> pair(3);
                    pair[0]=i;
                    pair[1]=j;
                    pair[2]=springlist[i].sticki;
                    springpairs.push_back(pair);
                    break;
                }
            }
        }
    }

}

void makeSticks(vector<stick> &mikado,vector<stick> &mikorig,const int NumberMikado,const double LStick)
{
    mikado=make_sticks(NumberMikado,LStick);
    mikorig=mikado; //The original set of sticks
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
}
    



spring makespring(int node1,int node2,double x1,double x2, double y1, double y2,int stick,double k,double stretchf)
{
    spring newspring;
    newspring.one=node1;
    newspring.two=node2;
    newspring.rlen=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))/stretchf;
    newspring.k=k;
    newspring.sticki=stick;
    
    if((x1<1 && x1>0)&&(x2>0&&x2<1)){ //Check if crossed wlr or wud wall.
        newspring.wlr=0;
    }
    else if((x1>0&&x1<1)&&x2>1){
        newspring.wlr=1;       
    }
    else if((x1>0&&x1<1)&&x2<0){
        newspring.wlr=-1;
    }
    else if(x1<0&&x2<0){
        newspring.wlr=0;
    }
    else if(x2>0 && x2<1.0 && x1<0.0){
        newspring.wlr=1;
    }
    else if(x2>0 && x2<1.0 && x1>1.0){
          newspring.wlr=-1;
    }
    else if(x1>1&&x2>1){
        newspring.wlr=0;
    }
    if((y1<1 && y1>0)&&(y2>0&&y2<1)){
        newspring.wud=0;
    }
    else if((y1>0&&y1<1)&&y2>1){
        newspring.wud=1;
    }
    else if((y1>0&&y1<1)&&y2<0){
        newspring.wud=-1;
    }
    else if(y1<0&&y2<0){
        newspring.wud=0;
    }
    else if(y1>1&&y2>1){
        newspring.wud=0;
    }
    else if(y2>0 && y2<1.0 && y1<0.0){
        newspring.wud=1;
    }
    else if(y2>0 && y2<1.0 && y1>1.0){
        newspring.wud=-1;
    }

    return newspring;
}

void makeanddeletebondsonbackground(vector<spring> &springlist,const vector<elonstick> &ELONSTICK,
                        stick &CURRENTSTICK,vector<double> &posonsticki,VectorXd &Xb,VectorXd &Yb,
                        int background_size,int sticknr,int elnr)
{
    int wud, wlr;
    int backnode1,backnode2,midpoint;
    double xbn1,xbn2,xm,ybn1,ybn2,ym;
    spring newspring;
    //spring to be removed has nodes
    backnode1=ELONSTICK[sticknr].backgroundspringn[elnr][0];
    backnode2=ELONSTICK[sticknr].backgroundspringn[elnr][1];
        for(int m=0;m<springlist.size();m++){
            if((springlist[m].one==backnode1 && springlist[m].two==backnode2) ||
                (springlist[m].one==backnode2 && springlist[m].two==backnode1)){
                //erase one spring.
                wud=springlist[m].wud;
                wlr=springlist[m].wlr;
                //parameters of the to be deleted spring 
                double springk=springlist[m].k;
                double springrlen=springlist[m].rlen;
                springlist.erase(springlist.begin()+m);   
                //replace it by two new ones. One between backnode1 and ELONSTICK[i].nr[j]+background_size
                //and ELONSTICK[i].nr[j]+backgroundsize and backgnode2
                xbn1=Xb(backnode1);
                ybn1=Yb(backnode1);
                xbn2=Xb(backnode2);
                ybn2=Yb(backnode2);
                xm=CURRENTSTICK.x+posonsticki[elnr]*cos(CURRENTSTICK.th);
                ym=CURRENTSTICK.y+posonsticki[elnr]*sin(CURRENTSTICK.th);
                midpoint=ELONSTICK[sticknr].nr[elnr]+background_size;
                //place the midpoint coordinate in the box
                if(xm>1){
                    xm=xm-1;
                }
                if(xm<0){
                    xm=xm+1;
                }
                if(ym>1){
                    ym=ym-1;
                }
                if(ym<0){
                    ym=ym+1;
                }
                //Check how the original coordinate should be 'fealt' by the midpoint coordinate
                if(wlr==1) {
                    if (fabs(xbn1-xm)<fabs(xbn2-xm)){
                        xbn2+=1;
                        } else {
                            xbn1-=1;
                    }
                }
                if(wlr==-1){
                    if(fabs(xbn1-xm)<fabs(xbn2-xm)){
                        xbn2-=1.0;
                    } else{
                        xbn1+=1.0; 
                    }
                }
                        
                        if(wud==1) {
                            if (fabs(ybn1-ym)<fabs(ybn2-ym)){
                                ybn2+=1;
                            } else {
                                ybn1-=1;
                            }
                        }
                        if(wud==-1){
                            if(fabs(ybn1-ym)<fabs(ybn2-ym)){
                                ybn2-=1.0;
                            } else{
                                ybn1+=1.0; 
                            }
                        }
                        double k1,k2;
                        double l1,l2;
                        l1=sqrt(pow((xbn1-xm),2)+pow((ybn1-ym),2));
                        l2=sqrt(pow((xbn2-xm),2)+pow((ybn2-ym),2));
                        k1=springrlen*springk/l1;
                        k2=springrlen*springk/l2;
                        newspring=makespring(backnode1,midpoint,xbn1,xm,ybn1,ym,-1,k1,1.0);
                        springlist.push_back(newspring);
                        newspring=makespring(midpoint,backnode2,xm,xbn2,ym,ybn2,-1,k2,1.0);
                        springlist.push_back(newspring);
                        break;
        }
    }
}

void ordernodes(vector<node> &nodes,vector<spring> &springlist){
    
    std::sort(nodes.begin(),nodes.end());

    vector<int> remove(0);
    
    for(int i=0;        i<nodes.size()-1;       i++){
        if(nodes[i].number==nodes[i+1].number){
            remove.push_back(i+1);
        }
    }
    
    //now erase the elements from the remove vector
    int del;
    while(remove.size()>0){
        del=remove[remove.size()-1];
        remove.pop_back();
        nodes.erase(nodes.begin()+del);
    }
    
    // next step: renumber and do the same to the springs
    int num;
    for(int i=0;        i<nodes.size(); i++){
        if(nodes[i].number!=i){
           num=nodes[i].number; 
           //correct in the springlist
           for(int j=0; j<springlist.size();    j++){
               if(springlist[j].one==num){
                   springlist[j].one=i;
               } 
               if(springlist[j].two==num){
                   springlist[j].two=i;
               }
            }
           nodes[i].number=i;
        }
    }
    

}




