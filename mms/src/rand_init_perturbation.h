#ifndef RAND_INIT_PERTURBATION_H
#define RAND_INIT_PERTURBATION_H

#include <random>

#include "typedefs.h"
#include "Vertex.h"


void rand_init_perturbation(std::vector<Vertex*>& vs, double arid)
{
    std::mt19937 rng(time(NULL));
    std::uniform_real_distribution<double> unif(-arid, arid);

    for (auto& v: vs)
    {
        for (uint i = 0; i < 3; ++i)
        {
            v->xyz[0] += unif(rng);
            v->xyz[1] += unif(rng);
            v->xyz[2] += unif(rng);
        }
    }
}

#endif
