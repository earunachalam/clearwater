#ifndef EXPORT_LAMMPSTRJ_H
#define EXPORT_LAMMPSTRJ_H

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "typedefs.h"
#include "Vertex.h"


void export_lammpstrj(FILE* f, uint tstep, uint natom, std::vector<Vertex*>& vs, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
    fprintf(f, "ITEM: TIMESTEP\n%u\n", tstep);
    fprintf(f, "ITEM: NUMBER OF ATOMS\n%u\n", natom);
    fprintf(f, "ITEM: BOX BOUNDS pp pp pp\n%f %f\n%f %f\n%f %f\n", xmin, xmax, ymin, ymax, zmin, zmax);
    fprintf(f, "ITEM: ATOMS id type x y z\n");

    for (uint iatom = 0; iatom < natom; ++iatom)
        fprintf(f, "%u %u %f %f %f\n", iatom+1, vs.at(iatom)->type, vs.at(iatom)->xyz[0], vs.at(iatom)->xyz[1], vs.at(iatom)->xyz[2]);
}



#endif
