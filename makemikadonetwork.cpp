#include "random.h"
#include <algorithm>
#include "makemikadonetwork.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <vector>
#include <iostream>
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
vector<connected> make_connections(const vector<stick> &m, double LStick){
  vector<connected> Connection;
  Matrix2d A; //The matrix to solve with
  Vector2d b,st; // A*st=b
  int nr=0;
  //Loop over all sticks to find coordinates
  for (int i=0;i<m.size()-1;i++){
    for (int j=i+1;j<m.size();j++){
      if(m[i].th!=m[j].th || m[i].nr!=m[j].nr){
        A<<-cos(m[i].th),cos(m[j].th),-sin(m[i].th),sin(m[j].th);
        b<<m[i].x-m[j].x,m[i].y-m[j].y;
        st=A.lu().solve(b);
        if ((st(0)>0 && st(0)<LStick)&&(st(1)>0 && st(1)<LStick)){
          connected xtrarow, xtrarow2;
          xtrarow.first=m[i].nr;	//[stick i stick j sij sji] 
          xtrarow.second=m[j].nr;
          xtrarow.s1=st(0);
          xtrarow.s2=st(1);
	  xtrarow.nrCon=nr;
	  xtrarow.recur=0;
	  xtrarow2.first=m[j].nr;	//[stick j sticki sji sij]
	  xtrarow2.second=m[i].nr;
	  xtrarow2.s1=st(1);
	  xtrarow2.s2=st(0);
	  xtrarow2.nrCon=nr;
	  xtrarow2.recur=0;
          Connection.push_back(xtrarow);
	  Connection.push_back(xtrarow2);
	  nr++;
        }
      }
    }
  }
  return Connection;
}

//Calculate per mikado how many and where other mikado's cross
vector<elonstick> sortELEMENTSperMIKADO(vector<connected> &Connection)
{
vector<elonstick> ELONSTICK;
for(int i=0;i<Connection.size()-1;i++){
    if(Connection[i].recur==1){
        elonstick extrarow;
        vector<double> extrarowpos(1);
        vector<int> extrarownumber(1);
        extrarowpos[0]=Connection[i].s1;
        extrarownumber[0]=Connection[i].nrCon;
    
        for(int j=i+1;j<Connection.size();j++){
            if(Connection[j].recur==1){
                if(Connection[i].first==Connection[j].first){
                    Connection[j].recur=0; 
                        //Connection[i].recur=0;  
                    extrarownumber.push_back(Connection[j].nrCon);
                    extrarowpos.push_back(Connection[j].s1);
                    extrarow.sticki=Connection[j].first;
                    extrarow.nr=extrarownumber;
                    extrarow.S=extrarowpos;
                }
            }
        }
ELONSTICK.push_back(extrarow);
   }
 }
// now sort extrarow on descending order per stick;
for(int j=0; j<ELONSTICK.size();j++){
    vector<double> distances=ELONSTICK[j].S;
    vector<int> numbers=ELONSTICK[j].nr;
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
	      swapped=1;
	      k++;
	    }
	  }
    } while(swapped==1);
	ELONSTICK[j].S=distances; //Put the new data back into the original vectors
	ELONSTICK[j].nr=numbers;
    }
std::sort(ELONSTICK.begin(),ELONSTICK.end());
return ELONSTICK;
}


//With the points per stick, we can now make nodes and springs.
void SpringsAndNodes(const vector<elonstick> &ELONSTICK,const vector<stick> &mikorig, vector<spring> &springlist,
vector<node> &nodes,
double rlenshort, double rlenlong,double k1,double k2)
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
            newspring.rlen=.1;
            newspring.k=1;

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

