#ifndef STRUCTS_H
#define STRUCTS_H
#include <vector>

struct stick{
  double x,y,th;
  int nr, wlr, wud;
};

struct connected{
  int first,second,nrCon,recur;
  double s1,s2;
  int type; //0 for mik mik connections 1 for mik-backgroundconnections
  int backgroundspring[2];
  
};


struct elonstick{
  int sticki;
  //vector<numberpos> NumberAndPos;
  std::vector<int> nr;
  std::vector<double> S;
  std::vector<int> type;
};
struct spring{
  int one, two, wlr, wud, sticki;
  double rlen,k;
};
struct node{
  //int number, sticki;
  int number;
  //int stick1;
  //int stick2;
  double x;
  double y;
};




#endif