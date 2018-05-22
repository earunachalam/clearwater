#ifndef VERTEX_H
#define VERTEX_H


#include <vector>


#include "core/typedefs.h"



class Vertex
{
    public:
        uint idx;                       // index of current vertex in vector of vertices
        uint type;                      // atom type
        bool frozen = false;            // if true, do not integrate motion
        std::vector<uint> neigh_idx;    // indices of neighboring points in vector of vertices
        std::vector<double> xyz,        // xyz position, global coordinates
                            vxyz,       // xyz velocity, global coordinates
                            neigh_R0,   // initial distances between self and each neighboring point, in same order as idx
                            neigh_Rt;   // current distances between self and each neighboring point, in same order as idx
        Vertex(uint idx, uint type=1) :
            idx(idx),    // set id
            type(type),  // set type
            xyz(3, 0.),  // initialize positions to zero
            vxyz(3, 0.)  // initialize velocities to zero
            {};
        Vertex(uint idx, std::vector<double> xyz, uint type=1) :
            idx(idx),    // set id
            type(type),  // set type
            xyz(xyz),    // initialize positions to zero
            vxyz(3, 0.)  // initialize velocities to zero
            {};
};


#endif
