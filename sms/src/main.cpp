#include <algorithm>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


#include "Bond.h"
#include "ExternalForce.h"
#include "Group.h"
#include "Vertex.h"


#include "export_lammpstrj.h"
#include "export_lengths.h"
#include "load_table.h"
#include "rand_init_perturbation.h"
#include "rel_pos.h"
#include "to_numeric.h"
#include "typedefs.h"



enum ENV_BOX_BEHAVIOR
{
    fit
} env_box_behavior;

void box_fit(
    double& xmin, double& xmax,
    double& ymin, double& ymax,
    double& zmin, double& zmax,
    const std::vector<Vertex*> vs
)
{
    if (env_box_behavior == fit)
    {
        xmin = vs.front()->xyz[0];
        xmax = vs.front()->xyz[0];
        ymin = vs.front()->xyz[1];
        ymax = vs.front()->xyz[1];
        zmin = vs.front()->xyz[2];
        zmax = vs.front()->xyz[2];

        for (size_t iv = 1; iv < vs.size(); ++ iv)
        {
            xmin = std::min(xmin, vs.at(iv)->xyz[0]);
            xmax = std::max(xmax, vs.at(iv)->xyz[0]);
            ymin = std::min(ymin, vs.at(iv)->xyz[1]);
            ymax = std::max(ymax, vs.at(iv)->xyz[1]);
            zmin = std::min(zmin, vs.at(iv)->xyz[2]);
            zmax = std::max(zmax, vs.at(iv)->xyz[2]);
        }
    }
}

class Integrator {};
class Overdamped_Integrator : public Integrator
{
    public:
        double visc = 0.;
};

int main(int argc, char* argv[])
{
    double max_ampl = std::stod(argv[1]);
    double lambda_star = std::stod(argv[2]);

    std::string inputfilename = "input.mms";
    size_t ntsteps = 1000,
           printvl = 100;
    double dt = 5e-3,
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


    // get neighbor lists
    std::ifstream inputfile(inputfilename);

    std::string line, token;

    while(std::getline(inputfile, line))
    {
        std::cout << line << "\n";
        std::istringstream iss(line);

        // a series of string tokens comprises a single instruction
        std::vector<std::string> instruction;
        while(std::getline(iss, token, ' '))
        {
            instruction.push_back(token);
        }

        if (instruction.empty())      continue;
        if (instruction.at(0) == "#") continue;

        if (instruction.at(0) == "box_behavior")
        {
            if (instruction.at(1) == "fill")
            {
                env_box_behavior = fit;
            }
        }
        else if (instruction.at(0) == "atom")
        {
            uint iatom = (uint) std::stoi(instruction.at(1));
            uint type  = (uint) std::stoi(instruction.at(2));
            double x = std::stod(instruction.at(3));
            double y = std::stod(instruction.at(4));
            double z = std::stod(instruction.at(5));
            vs.push_back(
                            new Vertex(iatom, {{x,y,z}}, type)
                        );
            
            
            for (size_t iopt = 6; iopt < instruction.size(); )
            {
                if (instruction.at(iopt) == "frozen")
                {
                    iopt++;

                    vs.back()->frozen = true;
                }
                else
                {
                    printf("Unrecognized option: %s.\n", instruction.at(iopt).c_str());
                    return 1;
                }
            }
        }
        else if (instruction.at(0) == "group")
        {
            uint id   = (uint) std::stoi(instruction.at(1));
            gs.push_back(new Group);
            auto g_ptr = gs.back();
            g_ptr->id = id;

            for (size_t iopt = 2; iopt < instruction.size(); )
            {
                if (instruction.at(iopt) == "range")
                {
                    iopt++;

                    uint start = (uint) std::stoi(instruction.at(iopt++));
                    uint end   = (uint) std::stoi(instruction.at(iopt++));

                    for (uint jatom = start; jatom <= end; ++jatom)
                    {
                        for (size_t kidx = 0; kidx < vs.size(); ++kidx)
                        {
                            if (vs.at(kidx)->idx == jatom)
                            {
                                g_ptr->atomic_members.push_back(kidx);
                                break;
                            }
                        }
                    }
                }
                else
                {
                    printf("Unrecognized option: %s.\n", instruction.at(iopt).c_str());
                    return 1;
                }
            }
        }
        else if (instruction.at(0) == "bond")
        {
            uint id   = (uint) std::stoi(instruction.at(1));
            uint iatom = (uint) std::stoi(instruction.at(2));
            uint jatom = (uint) std::stoi(instruction.at(3));

            size_t iopt = 4;
            if (instruction.at(iopt) == "harmonic")
            {
                iopt++;

                bs.push_back(new HarmonicBond);
                auto ibond_ptr = static_cast<HarmonicBond*>(bs.back());
                ibond_ptr->id   = id;

                bool ifound = false,
                     jfound = false;

                for (size_t i = 0; i < vs.size(); ++i)
                {
                    if (vs.at(i)->idx == iatom)
                    {
                        ibond_ptr->iatom = i;
                        ifound = true;
                        break;
                    }
                }
   
                if (!ifound)
                {
                    fprintf(stderr, "Error: Vertex with id %u not found.\n", iatom);
                    return 1;
                }

                for (size_t j = 0; j < vs.size(); ++j)
                {
                    if (vs.at(j)->idx == jatom)
                    {
                        ibond_ptr->jatom = j;
                        jfound = true;
                        break;
                    }
                }

                if (!jfound)
                {
                    fprintf(stderr, "Error: Vertex with id %u not found.\n", jatom);
                    return 1;
                }


                for ( ; iopt < instruction.size(); )
                {
                    if (instruction.at(iopt) == "L0")
                    {
                        iopt++;

                        if (instruction.at(iopt) == "initial_scale")
                        {
                            iopt++;
                            
                            std::vector<double> r1 = vs.at(ibond_ptr->iatom)->xyz;
                            std::vector<double> r2 = vs.at(ibond_ptr->jatom)->xyz;
                            std::vector<double> rp(4,0.);
                            rel_pos(r1, r2, rp);
                            
                            double scale   = std::stod(instruction.at(iopt++));
                            ibond_ptr->r0  = rp.at(3)*scale;
                            ibond_ptr->r   = rp.at(3);
                        }
                        if (instruction.at(iopt) == "initial")
                        {
                            iopt++;

                            std::vector<double> r1 = vs.at(ibond_ptr->iatom)->xyz;
                            std::vector<double> r2 = vs.at(ibond_ptr->jatom)->xyz;
                            std::vector<double> rp(4,0.);
                            rel_pos(r1, r2, rp);
                            
                            ibond_ptr->r0  = rp.at(3);
                            ibond_ptr->r   = rp.at(3);
                        }
                        if (instruction.at(iopt) == "k")
                        {
                            iopt++;

                            ibond_ptr->k  = std::stod(instruction.at(iopt++));
                            ibond_ptr->k0 = ibond_ptr->k;
                        }
                        else
                        {
                            printf("Unrecognized option: %s.\n", instruction.at(iopt).c_str());
                            return 1;
                        }
                    }
                    else if (instruction.at(iopt) == "addforce")
                    {
                        iopt++;

                        if (instruction.at(iopt) == "constant")
                        {
                            iopt++;

                            ibond_ptr->addforce = std::stod(instruction.at(iopt++));
                        }
                        else
                        {
                            printf("Unrecognized option: %s.\n", instruction.at(iopt).c_str());
                            return 1;
                        }
                    }
                    else
                    {
                        printf("Unrecognized option: %s.\n", instruction.at(iopt).c_str());
                        return 1;
                    }
                }
            }
        }
        else if (instruction.at(0) == "externalforce")
        {
            uint id       = (uint) std::stoi(instruction.at(1));
            uint group_id = (uint) std::stoi(instruction.at(2));
            Group* group;

            for (Group* igroup: gs)
            {
                if (igroup->id == group_id)
                {
                    group = igroup;
                    break;
                }
            }

            efs.push_back(new ExternalForce(id, group_id, group));
            ExternalForce* e = efs.back();

            size_t iopt = 3;
            if (instruction.at(iopt) == "type")
            {
                iopt++;
                
                if (instruction.at(iopt) == "table")
                {
                    iopt++;

                    std::string dirname = instruction.at(iopt++);
 
                    std::string cname   = dirname + "/config.dat";
 
                    std::string xname   = dirname + "/x.dat";
                    std::string yname   = dirname + "/y.dat";
                    //std::string zname   = dirname + "/z.dat";
 
                    std::string fxname  = dirname + "/fx.dat";
                    std::string fyname  = dirname + "/fy.dat";
                    //std::string fzname  = dirname + "/fz.dat";
                    
                    uint nrows, ncols;
                    std::vector<std::vector<double> > config;

                    load_table(cname, config, nrows, ncols);
                    
                    e->nx    = config[0][0] - 1;
                    e->ny    = config[0][1] - 1;
                    //e->nz    = config[0][2] - 1;
                    
                    e->xstep = config[2][0];
                    e->ystep = config[2][1];
                    //e->zstep = config[2][2];

                    e->xmin  = config[1][0] + 0.5*e->xstep;
                    e->ymin  = config[1][1] + 0.5*e->ystep;
                    //e->zmin  = config[1][2] + 0.5*e->zstep;

                    load_flattened_table(xname, e->x, nrows, ncols);
                    load_flattened_table(yname, e->y, nrows, ncols);
                    //load_flattened_table(zname, e->z, nrows, ncols);

                    load_flattened_table(fxname, e->fx, nrows, ncols);
                    load_flattened_table(fyname, e->fy, nrows, ncols);
                    //load_flattened_table(fzname, e->fz, nrows, ncols);
                }
            }
        }
        else if (instruction.at(0) == "integrator")
        {
            size_t iopt = 1;
            if (instruction.at(iopt) == "overdamped")
            {
                iopt++;

                integ = new Overdamped_Integrator();

                for ( ; iopt < instruction.size(); )
                {
                    if (instruction.at(iopt) == "visc")
                    {
                        iopt++;
                        static_cast<Overdamped_Integrator*>(integ)->visc = std::stod(instruction.at(iopt++));
                    }
                }
            }
        }
        else if (instruction.at(0) == "timestep")
        {
            dt = std::stod(instruction.at(1));
        }
        else if (instruction.at(0) == "printvl")
        {
            printvl = std::stod(instruction.at(1));
        }
        else if (instruction.at(0) == "run")
        {
            ntsteps = std::stod(instruction.at(1));
        }
        else
        {
            printf("Unrecognized command: %s.\n", instruction.at(0).c_str());
            return 1;
        }
                    
        //else if (instruction.at(0) == "integration")
        //{
            //if (instruction.at(1) == "overdamped")
            //{
                //// viscosity
                
                //double g = std::stod(instruction.at(2));
            //}
        //}
    }
    
    inputfile.close();

    natoms = vs.size();
    nbonds = bs.size();
    nextfs = efs.size();

    //abort();


 

    FILE* ofp_traj   = fopen("traj.lammpstrj", "w");
    FILE* ofp_length = fopen("length.dat", "w");
    FILE* ofp_torque = fopen("torque.dat", "w");
    
    // export initial coordinates
    if (env_box_behavior == fit) (box_fit(xmin, xmax, ymin, ymax, zmin, zmax, vs));
    export_lammpstrj(ofp_traj, 0, natoms, vs, xmin, xmax, ymin, ymax, zmin, zmax);
    export_lengths(ofp_length, 0, natoms, bs);

    // integrate
    for (uint tstep = 1; tstep <= ntsteps; ++tstep)
    {
        // calculate center of mass
        double x_cm = 0., y_cm = 0., z_cm = 0.;
        for (size_t iatom = 0; iatom < natoms; ++iatom)
        {
            x_cm += vs.at(iatom)->xyz[0];
            y_cm += vs.at(iatom)->xyz[1];
            z_cm += vs.at(iatom)->xyz[2];
        }
        x_cm /= natoms; y_cm /= natoms; z_cm /= natoms;
        double vx_cm = 0., vy_cm = 0., vz_cm = 0.;
        for (size_t iatom = 0; iatom < natoms; ++iatom)
        {
            vx_cm += vs.at(iatom)->xyz[0];
            vy_cm += vs.at(iatom)->xyz[1];
            vz_cm += vs.at(iatom)->xyz[2];
        }
        vx_cm /= natoms; vy_cm /= natoms; vz_cm /= natoms;
 
        // variable to hold torque x, y, z components       
        double torque[4];
        for (size_t icoord = 0; icoord < 4; ++icoord) torque[icoord] = 0.;

        // + initialize velocities to zero at beginning of each timestep
        // before summing all contributions to force from neighbors
        // and background force
        // + velocity is therefore memoryless (noninertial dynamics)
        for (size_t iatom = 0; iatom < natoms; ++iatom)
            std::fill(vs.at(iatom)->vxyz.begin(), vs.at(iatom)->vxyz.end(), 0.);

//#pragma omp parallel for
        for (size_t ibond = 0; ibond < nbonds; ++ibond)
        {
            auto iBond = static_cast<EvolvingHarmonicBond*>(bs.at(ibond));
            uint iatom = iBond->iatom;
            uint jatom = iBond->jatom;
            
            if (vs.at(iatom)->frozen && vs.at(jatom)->frozen)
            {
                continue;
            }

            std::vector<double> r1 = vs.at(iatom)->xyz;
            std::vector<double> r2 = vs.at(jatom)->xyz;
            std::vector<double> rp(4,0.);
            
            rel_pos(r1, r2, rp);
            iBond->r = rp[3];
            iBond->evolve_r0(max_ampl, lambda_star, dt);

            // add harmonic restoring force
            double deltaR = iBond->r - iBond->r0;
            for (uint icoord = 0; icoord < 3; ++icoord)
            {
                if (!vs.at(iatom)->frozen)
                {
                    vs.at(iatom)->vxyz.at(icoord) -= (iBond->k*deltaR + iBond->addforce*rp[3])*rp.at(icoord);
                }
                if (!vs.at(jatom)->frozen)
                {
                    vs.at(jatom)->vxyz.at(icoord) += (iBond->k*deltaR + iBond->addforce*rp[3])*rp.at(icoord);
                }
            }
        }

        //for (size_t iextf = 0; iextf < nextfs; ++iextf)
        //{
            //ExternalForce* iExtF        = efs.at(iextf);
            //std::vector<uint> atom_idxs = iExtF->group->atomic_members;

            //for (size_t jidx = 0; jidx < atom_idxs.size(); ++jidx)
            //{
                //uint jatom = atom_idxs.at(jidx);
                
                //std::vector<double> fxyz(3);
                //iExtF->interpolate_force(vs.at(jatom)->xyz, fxyz);

                //std::cout << vs.at(jatom)->vxyz[0] << " " << vs.at(jatom)->vxyz[1] << " " << vs.at(jatom)->vxyz[2] << "\n";
                //std::cout << fxyz[0] << " " << fxyz[1] << " " << fxyz[2] << "\n";
                //double a = 1e-16;
                //vs.at(jatom)->vxyz[0] += a*fxyz[0];
                //vs.at(jatom)->vxyz[1] += a*fxyz[1];
                //vs.at(jatom)->vxyz[2] += a*fxyz[2];
            //}
        //}

        for (size_t iatom = 0; iatom < natoms; ++iatom)
        {
            auto ri = vs.at(iatom)->xyz;
            auto vi = vs.at(iatom)->vxyz;
            torque[0] += (vi[1]-vy_cm)*(ri[2]-z_cm) - (vi[2]-vz_cm)*(ri[1]-y_cm);
            torque[1] += (vi[2]-vz_cm)*(ri[0]-x_cm) - (vi[0]-vx_cm)*(ri[2]-z_cm);
            torque[2] += (vi[0]-vx_cm)*(ri[1]-y_cm) - (vi[1]-vy_cm)*(ri[0]-x_cm);
        }

        for (size_t iatom = 0; iatom < natoms; ++iatom)
        {
            for (uint icoord = 0; icoord < 3; ++icoord)
            {
                vs.at(iatom)->vxyz.at(icoord) *= static_cast<Overdamped_Integrator*>(integ)->visc;
            }
        }

        //if (tstep == 2)
            //for (int z = 0; z < (int) natoms; ++z)
            //{
                //printf("%f %f %f %f\n",
                        //vs.at(z)->xyz[0],
                        //vs.at(z)->xyz[1],
                        //vs.at(z)->vxyz[0],
                        //vs.at(z)->vxyz[1]);
            //}

        for (size_t iatom = 0; iatom < natoms; ++iatom)
        {
            for (uint icoord = 0; icoord < 3; ++icoord)
            {
                vs.at(iatom)->xyz.at(icoord) += dt*vs.at(iatom)->vxyz.at(icoord);
            }
        }

        if ((tstep%printvl == 0) || (tstep==ntsteps))
        {
            if (env_box_behavior == fit) box_fit(xmin, xmax, ymin, ymax, zmin, zmax, vs);

            export_lammpstrj(ofp_traj, tstep, natoms, vs, xmin, xmax, ymin, ymax, zmin, zmax);
            export_lengths(ofp_length, tstep, natoms, bs);
        }
    
        torque[3] = std::sqrt(torque[0]*torque[0] + torque[1]*torque[1] + torque[2]*torque[2]);
        fprintf(ofp_torque, "%d %f\n", tstep, torque[3]);

        if (tstep%(ntsteps/100) == 0)
            printf("%d\n", tstep);
    }

    fclose(ofp_traj);
    fclose(ofp_length);
    fclose(ofp_torque);

    // free memory
    for (uint iatom = 0; iatom < natoms; ++iatom)
        delete vs.at(iatom);
    
    return 0;
}
