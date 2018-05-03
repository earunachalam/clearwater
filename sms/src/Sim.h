#ifndef SIMULATION_H
#define SIMULATION_H


class Vertex;
class Bond;
class Integrator;

class Simulation
{
    size_t ntsteps,
           printvl;
    double dt,
           L = 10.;    
    double xmin = -L, xmax = L;
    double ymin = -L, ymax = L;
    double zmin = -L, zmax = L;

    Integrator* integ;


    std::vector<Vertex*> vs;
    std::vector<Bond*> bs;
    std::vector<Group*> gs;
    std::vector<ExternalForce*> efs;
    size_t natoms,
           nbonds,
           nextfs;
}


#endif
