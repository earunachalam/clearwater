#ifndef EXPORT_LENGTHS_H
#define EXPORT_LENGTHS_H

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "typedefs.h"
#include "Bond.h"


void export_lengths(FILE* f, uint tstep, uint natom, std::vector<Bond*>& bs)
{
    fprintf(f, "%u ", tstep);
    for (uint ibond = 0; ibond < bs.size(); ++ibond)
        fprintf(f, "%f ", bs.at(ibond)->r);
    fprintf(f, "\n");
}


#endif
