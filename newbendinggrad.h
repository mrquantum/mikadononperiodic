#ifndef NEWBENDINGGRAD_h
#define NEWBENDINGGRAD_h

#include <iostream>
#include "structs.h"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/Sparse"

void bendinggradnew(double* xy,double *grad,int size,std::vector<spring> &springlist,std::vector<std::vector<int> > &springpairs,double b_rigid,double gamma);



#endif