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
std::vector<stick> make_sticks(int N)
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
  for (int i=0;i<m.size()-1;i++){
    for (int j=i+1;j<m.size();j++){
      if(m[i].th!=m[j].th || m[i].nr!=m[j].nr){
            A<<-cos(m[i].th),cos(m[j].th),-sin(m[i].th),sin(m[j].th);
            b<<m[i].x-m[j].x,m[i].y-m[j].y;
            st=A.lu().solve(b);
            if ((st(0)>0.0 && st(0)<LStick)&&(st(1)>0.0 && st(1)<LStick)){
                xtrarow.first=m[i].nr;	//[stick i stick j sij sji] 
                xtrarow.second=m[j].nr;
		cout<<st(0)<<"\t"<<st(1)<<endl;
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
    cout<<"gethere1"<<endl;
    for(int i=0; i<m.size();i++){
        xmik=m[i].x;//+m[i].wlr;
        ymik=m[i].y;//+m[i].wud;
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
  
  cout<<"WE FOUND THE FOLLOWING CONNECTIONS BETWEEN SPRINGS AND MIKADO"<<endl;
  cout<<"Miknr1  Miknr2    NRCon        type    backgroundnode1 node2   s1      s2"<<endl;
for(int i=0;i<Connection.size();i++){
    cout<<Connection[i].first<<"\t"<<Connection[i].second<<"\t"
    <<Connection[i].nrCon<<"\t"<<Connection[i].type<<"\t"<<Connection[i].backgroundspring[0]<<"\t"
    <<Connection[i].backgroundspring[1]<<"\t"<<Connection[i].s1<<"\t"<<Connection[i].s2<<endl;
}
  cout<<"END"<<endl;
}


//Tomorrow make sure that ALL the elements are include 
//in ELONSTICK 
void sortELEMENTSperMIKADO(vector<elonstick> &ELONSTICK,vector<connected> &Connection)
{
    
   // cout<<Connection.size()<<"SSSS"<<endl;

  combineElementsOnMikado(Connection,ELONSTICK);

  for(int i=0;i<ELONSTICK.size();i++){
   cout<<"MIK NR=       "<<ELONSTICK[i].sticki<<endl;
   
    for(int j=0;j<ELONSTICK[i].S.size();j++){
      cout<<ELONSTICK[i].nr[j]<<"	";
    };
    cout<<endl<<"**"<<endl;
    for(int j=0;j<ELONSTICK[i].S.size();j++){
     cout<<ELONSTICK[i].S[j]<<"	";
    };
    cout<<endl<<"**"<<endl;
    for(int j=0;j<ELONSTICK[i].S.size();j++){
     cout<<ELONSTICK[i].type[j]<<"	";
    };
    cout<<endl;
   }
  
  
// elonstick extrarow;
// int sticki;
// 
// for(int i=0;i<Connection.size()-1;i++){
//   vector<double> extrarowpos(1);
//   vector<int> extrarownumber(1);
//   vector<int> elementtype(1);
//   
//   if(Connection[i].type==0){
//       if(Connection[i].recur==1){
//             
//             extrarowpos[0]=Connection[i].s1;
//             extrarownumber[0]=Connection[i].nrCon;
//             elementtype[0]=0;
//             sticki=Connection[i].first;
//             
//             for(int j=i+1;j<Connection.size();j++){
//                 if(Connection[j].recur==1 && Connection[j].type==0){
//                     if(Connection[i].first==Connection[j].first){
//                         Connection[j].recur=0; 
//                         extrarownumber.push_back(Connection[j].nrCon);
//                         extrarowpos.push_back(Connection[j].s1);
//                         elementtype.push_back(0);
//                     }
//                 }
//             }
//        }
//   }
// 
//     if(Connection[i].type==1){
//         extrarownumber.push_back(Connection[i].nrCon);
//         extrarowpos.push_back(Connection[i].s1);
//         elementtype.push_back(1);
//         sticki=Connection[i].first;
//     }
//        extrarow.sticki=sticki;
//        extrarow.nr=extrarownumber;
//        extrarow.S=extrarowpos;
// /*       cout<<"stick     "<<sticki<<endl;
//        for(int m=0; m<extrarowpos.size();m++)
//        {cout<<extrarowpos[m]<<"\t";
//        }
//        cout<<endl*/;
//        extrarow.type=elementtype;   
//        ELONSTICK.push_back(extrarow);
// 
// }



// now sort extrarow on descending order per stick;
for(int j=0; j<ELONSTICK.size();j++){
    vector<double> distances=ELONSTICK[j].S;
    vector<int> numbers=ELONSTICK[j].nr;
    vector<int> type=ELONSTICK[j].type;
// now sort extrarow on descending order per stick;
    int swapped=0;
      do{
        int k=0;
        swapped=0; //this is the control parameter, checks 1 if elements are swapped
          for(int i=0;i<distances.size()-1;i++){ //loop through the list thill the end-k-1 th element;
            if(distances[i]>distances[i+1]){ //checks if neighbours are in right order, if not then swap en change swap parameter
              double aa=distances[i];
              double bb=distances[i+1];
              distances[i]=bb;
              distances[i+1]=aa;
              int a=numbers[i]; 
              int b=numbers[i+1];
              numbers[i]=b; 
              numbers[i+1]=a;
              int aaa=type[i];
	      int bbb=type[i+1];
	      type[i]=bbb;
	      type[i]=aaa;
	            
	      swapped=1;
              k++;
            }
          }
    } while(swapped==1);
        ELONSTICK[j].S=distances; //Put the new data back into the original vectors
        ELONSTICK[j].nr=numbers;
	ELONSTICK[j].type=type;
    }
std::sort(ELONSTICK.begin(),ELONSTICK.end());


//   for(int i=0;i<ELONSTICK.size();i++){
//    cout<<"MIK NR=       "<<ELONSTICK[i].sticki<<endl;
//    
//     for(int j=0;j<ELONSTICK[i].S.size();j++){
//       cout<<ELONSTICK[i].nr[j]<<"    ";
//     };
//     cout<<endl<<"**"<<endl;
//     for(int j=0;j<ELONSTICK[i].S.size();j++){
//      cout<<ELONSTICK[i].S[j]<<"      ";
//     };
//     cout<<endl<<"**"<<endl;
//     for(int j=0;j<ELONSTICK[i].S.size();j++){
//      cout<<ELONSTICK[i].type[j]<<"   ";
//     };
//     cout<<endl;
//    }



// for(int j=0;j<ELONSTICK.size();j++){
//         cout<<"STICKNR  S       type    nrcon"<<endl;
//         cout<<"On Stick "<<ELONSTICK[j].sticki<<endl;
//         for(int m=0;m<ELONSTICK[j].S.size();m++){
//              cout<<"\t"<<ELONSTICK[j].S[m]<<"\t"<<ELONSTICK[j].type[m]<<"\t"<<ELONSTICK[j].nr[m]<<endl;
//         }
//     }
    
//somehow on some elements the sequence {0,0,0} occurs. 
//I cant find an elegant way to remove it, so let's bruteforce it:

// vector<elonstick> ELONSTICK2(0);
// for(int i=0;i<ELONSTICK.size();i++){
//     if(ELONSTICK[i].S[0]!=0){
//         ELONSTICK2.push_back(ELONSTICK[i]);
//     }
// }
// cout<<"****"<<endl;

// for(int j=0;j<ELONSTICK2.size();j++){
//         cout<<"STICKNR  S       type    nrcon"<<endl;
//         cout<<"On Stick "<<ELONSTICK2[j].sticki<<endl;
//         for(int m=0;m<ELONSTICK2[j].S.size();m++){
//              cout<<"\t"<<ELONSTICK2[j].S[m]<<"\t"<<ELONSTICK2[j].type[m]<<"\t"<<ELONSTICK2[j].nr[m]<<endl;
//         }
//     }

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
double rlenshort, double rlenlong,double k1,double k2,double stretchf)
{
for(int i=0;i<ELONSTICK.size();i++){
    int sticknr=ELONSTICK[i].sticki;
    vector<int> nodesonsticki=ELONSTICK[i].nr;
    vector<double> posonsticki=ELONSTICK[i].S;
    stick CURRENTSTICK=mikorig[sticknr];
        
        for(int j=0;j<nodesonsticki.size()-1;j++){
            spring newspring;
            newspring.one=nodesonsticki[j];
            newspring.two=nodesonsticki[j+1];
            
            double x1, x2, y1, y2;
            x1=CURRENTSTICK.x+posonsticki[j]*cos(CURRENTSTICK.th); //calculate the position of the node
            x2=CURRENTSTICK.x+posonsticki[j+1]*cos(CURRENTSTICK.th);//and the position of the adjacent one
            y1=CURRENTSTICK.y+posonsticki[j]*sin(CURRENTSTICK.th);
            y2=CURRENTSTICK.y+posonsticki[j+1]*sin(CURRENTSTICK.th);
            
//             if(sticknr%2==0){
//                 newspring.rlen=rlenshort;
//                 newspring.k=k1;
//             }
//             else{
//                 newspring.rlen=rlenlong;
//                 newspring.k=k2;
//             }
            newspring.rlen=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))/stretchf;
            newspring.k=k1;
//             if(sticknr%2==0){
//             newspring.k=10;}
//             else{
//                 newspring.k=1000;
//             }
//            newspring.rlen=.1;
//            newspring.k=1;

            newspring.sticki=sticknr;
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
            node nodetemp1, nodetemp2;
            nodetemp1.number=newspring.one;
            nodetemp1.x=x1;
            nodetemp1.y=y1;
            nodetemp2.number=newspring.two;
            nodetemp2.x=x2;
            nodetemp2.y=y2;
            nodes.push_back(nodetemp1);
            nodes.push_back(nodetemp2);
            springlist.push_back(newspring);
        }
    }
std::sort(nodes.begin(),nodes.end());

}

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
for(std::size_t i=0;i<springlist.size()-1;i++){
    if(springlist[i].sticki==springlist[i+1].sticki){
        vector<int> pair(3);
        pair[0]=i;
        pair[1]=i+1;
        pair[2]=springlist[i].sticki;
        springpairs.push_back(pair);
    }   
}
}

void makeSticks(vector<stick> &mikado,vector<stick> &mikorig,const int NumberMikado,const double LStick)
{
    mikado=make_sticks(NumberMikado);
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
    

void makeConnections(vector<connected> &Connection,
                     const vector<stick> &mikado,
                     const double LStick,
                     const vector<spring> &background,
                     const Eigen::VectorXd &XYb)
{
    make_connections(Connection,
                     mikado,
                     LStick,
                     background,
                     XYb
                    ); //Make Connections
//     vector<connected> Connection2(1); //sives the double elements from the connections
//     Connection2[0]=Connection[0];
//     for(std::size_t i=0;i<Connection.size();i++){
//         int flag=1;
//         for(std::size_t j=0;j<Connection2.size();j++){
//             if(Connection[i].first==Connection2[j].first && Connection[i].second==Connection2[j].second){
//                 flag=0; 
//                 break;
//             }
//         }
//         if(flag==1){
//             Connection2.push_back(Connection[i]);
//         } 
//     }
// Connection=Connection2;

//Check for pairs in the connected struct if at least a pair exists (2 points on same mikado)
//then connection[j].recur=1; Else recur =0.
//     for(std::size_t i=0;i<Connection.size()-1;i++){
//         for(std::size_t j=i+1;j<Connection.size();j++){
//             if(Connection[i].first==Connection[j].first){
//                 Connection[i].recur=1;  Connection[j].recur=1;
//             }
//         }
//     }
}


