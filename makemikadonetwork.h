#ifndef MAKEMIKADONETWORK_H
#define MAKEMIKADONETWORK_H

#include<vector>

struct stick{
  double x,y,th;
  int nr, wlr, wud;
};

struct connected{
  int first,second,nrCon,recur;
  double s1,s2;
};


struct elonstick{
  int sticki;
  //vector<numberpos> NumberAndPos;
  std::vector<int> nr;
  std::vector<double> S;
};
struct spring{
  int one, two, wlr, wud, sticki;
  double rlen,k;
};
struct node{
  //int number, sticki;
  int number;
  double x;
  double y;
};




std::vector<stick> make_sticks(int N);
std::vector<stick> make_ghost_lr( const std::vector<stick> &m, double LStick, int NumberMikado);
std::vector<stick> make_ghost_ud(const std::vector<stick> &m, double LStick, int NumberMikado);

void make_connections(std::vector<connected> &Connection,const std::vector<stick> &m, double LStick);

void sortELEMENTSperMIKADO(std::vector<elonstick> &ELONSTICK,std::vector<connected> &Connection);
void orderElonstick(std::vector<int> &order,std::vector<elonstick> &ELONSTICK);

void makeSpringsAndNodes(const std::vector<elonstick> &ELONSTICK,
		     const std::vector<stick> &mikorig, std::vector<spring> &springlist, std::vector<node> &nodes,
		     double rlenshort,double rlenlong,double k1, double k2,double stretchf);
bool operator<(const node& first,const node& second);
bool operator<(const elonstick &first, const elonstick &second);
bool operator<(const stick &first, const stick &second);

void makeSticks(std::vector<stick> &mikado,std::vector<stick> &mikorig,const int NumberMikado,const double LStick);

void makeConnections(std::vector<connected> &Connection,
                     const std::vector<stick> &mikado,
                     const double LStick);



double inbox(double x,double boxsize);
void makeSpringpairs(std::vector<std::vector<int>> &springpairs,const std::vector<spring> &springlist);


#endif // MAKEMIKADONETWORK_H