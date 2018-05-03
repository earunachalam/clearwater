#ifndef BOND_H
#define BOND_H


#include <cmath>


#include "sgn.h"
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

class HarmonicBond : public Bond
{
    public:
        double k  = -1.;
        double k0 = -1.;
        double r0 = -1.;
};

class EvolvingHarmonicBond : public HarmonicBond
{
    public:
        void evolve_r0(double max_ampl, double lambda_star, double dt)
        {
            double deviation = this->r - this->r0;
            double s = sgn(deviation);
            //deviation = std::fabs(deviation);

            double g = s*max_ampl*deviation;
            //double g = (deviation < lambda_star) ? 
                //s*max_ampl : 
                //s*max_ampl*std::exp(-(deviation-lambda_star));

            this->r0 += g*deviation*dt;
        }

};


#endif
