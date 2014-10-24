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

struct numberpos{
  int i; 
  double s;
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
  int number;
  double x;
  double y;
};

struct triplet{
  int one, two, three, sticki;
};



std::vector<stick> make_sticks(int N);
std::vector<stick> make_ghost_lr( const std::vector<stick> &m, double LStick, int NumberMikado);
std::vector<stick> make_ghost_ud(const std::vector<stick> &m, double LStick, int NumberMikado);
std::vector<connected> make_connections(const std::vector<stick> &m, double LStick);
std::vector<elonstick> sortELEMENTSperMIKADO(std::vector<connected> &Connection);
void SpringsAndNodes(const std::vector<elonstick> &ELONSTICK,
		     const std::vector<stick> &mikorig, std::vector<spring> &springlist, std::vector<node> &nodes,
		     double rlenshort,double rlenlong,double k1, double k2);
bool operator<(const node& first,const node& second);
bool operator<(const elonstick& first, const elonstick& second);




#endif // MAKEMIKADONETWORK_H