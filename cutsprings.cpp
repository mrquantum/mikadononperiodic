#include <iostream>
#include <ctime>
#include "structs.h"
#include <vector>
#include "fstream"
#include "math.h"
#include "random.h"
#include <algorithm>


using namespace std;
void cutspring(vector<spring> &springlist,int number){
    
    vector<int> tobedeleted;
    for(int i=0;i<number;i++){
        tobedeleted.push_back(randi<std::size_t>(0,springlist.size()));
    }

    sort(tobedeleted.begin(),tobedeleted.end());
 

    int eltobedeleted;
    while(tobedeleted.size()>1){
        eltobedeleted=tobedeleted[tobedeleted.size()-1];
        cout<<eltobedeleted<<endl;
        springlist.erase(springlist.begin()+eltobedeleted);
        tobedeleted.pop_back();
    }
 
}
