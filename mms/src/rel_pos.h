#ifndef REL_POS_H
#define REL_POS_H

#include <cmath>
#include <cstdlib>
#include <vector>

#include "typedefs.h"


void rel_pos(std::vector<double>& xyz1, std::vector<double>& xyz2, std::vector<double>& rel_pos12)
{
    double x = xyz1[0] - xyz2[0];
    double y = xyz1[1] - xyz2[1];
    double z = xyz1[2] - xyz2[2];

    double r = std::sqrt(x*x + y*y + z*z);
    //double t = std::acos(z/r);
    //double p = std::atan2(y,x);

    rel_pos12[0] = x/r; // std::sin(t)*std::cos(p);
    rel_pos12[1] = y/r; // std::sin(t)*std::sin(p);
    rel_pos12[2] = z/r; // std::cos(t);
    rel_pos12[3] = r;
}



#endif
