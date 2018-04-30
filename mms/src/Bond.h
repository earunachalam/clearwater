#ifndef BOND_H
#define BOND_H

#include "typedefs.h"

class Bond
{
    public:
        uint id;
        double r  = -1.,
               addforce = 0.;
        uint iatom,
             jatom;
};

class BondHarmonic : public Bond
{
    public:
        double k  = -1.;
        double k0 = -1.;
        double r0 = -1.;
};


#endif
