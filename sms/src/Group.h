#ifndef GROUP_H
#define GROUP_H


#include <vector>


#include "core/typedefs.h"



class Group
{
    public:
        uint id;
        std::vector<uint> atomic_members;
};


#endif
