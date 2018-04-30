#ifndef EXTERNALFORCE_H
#define EXTERNALFORCE_H


#include <cmath>
#include <cstdlib>
#include <vector>


#include "Group.h"


#include "typedefs.h"



class ExternalForce
{
    public:

        uint id       = -1;
        uint group_id = -1;
        Group* group;
        
        ExternalForce(uint id, uint group_id, Group* group): id(id), group_id(group_id), group(group) {};
 
        // number of x, y, z mesh points       
        uint nx = -1,
             ny = -1,
             nz = -1;

        // x, y, z limits
        double xmin, ymin, zmin,
               xstep, ystep, zstep;

        std::vector<double> x,
                            y,
                            z;
 
        std::vector<double> fx,
                            fy,
                            fz;

        void find_idx(std::vector<double> xyz, std::vector<double>& rem, uint& ceil, uint& floor)
        {
            rem[0] = (this->xstep == 0.) ? 0. : fmod(xyz[0]-this->xmin, this->xstep);
            rem[1] = (this->ystep == 0.) ? 0. : fmod(xyz[1]-this->ymin, this->ystep);
            rem[2] = (this->zstep == 0.) ? 0. : fmod(xyz[2]-this->zmin, this->zstep);

            // printf("%f %f %f\n", xyz[0]-this->xmin, xyz[1]-this->ymin, xyz[2]-this->zmin);
            // printf("%f %f %f\n", xstep, ystep, zstep);
            // printf("%f %f %f\n", rem[0], rem[1], rem[2]);
            // abort();

            uint xidx, yidx, zidx;

            xidx = (uint) std::ceil((xyz[0]-this->xmin)/this->xstep);
            yidx = (uint) std::ceil((xyz[1]-this->ymin)/this->ystep);
            zidx = (uint) std::ceil((xyz[2]-this->zmin)/this->zstep);
            ceil  = zidx*this->nx*this->ny + yidx*this->nx + xidx;
// std::cout << "zidx = " << zidx << std::endl;
            xidx = (uint) std::floor((xyz[0]-this->xmin)/this->xstep);
            yidx = (uint) std::floor((xyz[1]-this->ymin)/this->ystep);
            zidx = (uint) std::floor((xyz[2]-this->zmin)/this->zstep);
            floor = zidx*this->nx*this->ny + yidx*this->nx + xidx;
// std::cout << "zidx = " << zidx << std::endl;
// std::cout << "rem[2] = " << rem[2] << std::endl;
//             std::cout << "z = " << xyz[2] << std::endl;

//             abort();
        }

        void interpolate_force(std::vector<double> xyz, std::vector<double>& fxyz)
        {
            std::vector<double> rem(3);
            uint ceil, floor;

            find_idx(xyz, rem, ceil, floor);
            printf("cf = %u %u\n", ceil, floor);
            
            fxyz[0] = fx[floor] + rem[0]*(fx[ceil]-fx[floor]);
            fxyz[1] = fy[floor] + rem[1]*(fy[ceil]-fy[floor]);
            fxyz[2] = fz[floor] + rem[2]*(fz[ceil]-fz[floor]);
        }


};

#endif
