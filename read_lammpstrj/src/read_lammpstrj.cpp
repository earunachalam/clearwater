#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>


#include <core/common.h>
#include <core/numeric.h>


namespace nm = numeric;

int main(int argc, char** argv)
{
	uint natoms;

	if (argc == 2)
	{
		natoms = std::atoi(argv[1]);
		printf("natoms = %d\n",natoms);
	}
	else
	{
		printf("Syntax: \n");
		return 1;
	}


	file::infile traj("traj.lammpstrj");
	FILE* ofp_thick = fopen("thick.dat", "w");

	//std::string str;
	uint timestep, curr_natoms;
	std::vector<std::string> strs;
	std::vector<double> vec;


	while (traj.is_open())
	{
		traj.readline();
		traj.readline(timestep);
		traj.readline();
		traj.readline(curr_natoms);

		if (curr_natoms != natoms)
		{
			printf("Error: atom number mismatch.\n");
			return 2;
		}

		traj.readline();
		traj.readline();
		traj.readline();
		traj.readline();
		traj.readline();

		std::vector<int> n(natoms), t(natoms);
		std::vector<double> x(natoms), y(natoms), z(natoms);

		for (uint iatom = 0; iatom < natoms; ++iatom)
		{
			traj.readline(vec);
			n[iatom] = vec[0];
			t[iatom] = vec[1];
			x[iatom] = vec[2];
			y[iatom] = vec[3];
			z[iatom] = vec[4];
		}

		std::vector<double> x1 = use(x, t==1);
		std::vector<double> x2 = use(x, t==2);

		std::vector<double> y1 = use(y, t==1);
		std::vector<double> y2 = use(y, t==2);

		std::vector<double> z1 = use(z, t==1);
		std::vector<double> z2 = use(z, t==2);

		unsigned int i1, i2;
		double dx, dy, dz, thick_side, thick_pole;

		std::vector<double> x1abs = nm::abs(x1);
		std::vector<double> x2abs = nm::abs(x2);

		std::vector<double> y1abs = nm::abs(y1);
		std::vector<double> y2abs = nm::abs(y2);

		std::vector<double> z1abs = nm::abs(z1);
		std::vector<double> z2abs = nm::abs(z2);

		i1 = nm::argmin(z1abs);
		i2 = nm::argmin(z2abs);
		
		dx = x2abs[i2] - x1abs[i1];
		dy = y2abs[i2] - y1abs[i1];
		//dz = z2abs[i2] - z1abs[i1];

		thick_side = std::sqrt(dx*dx + dy*dy);
		
		std::vector<double> x1gt0 = use(x1, x1>0.);
		std::vector<double> x2gt0 = use(x2, x2>0.);

		std::vector<double> y1gt0 = use(y1, y1>0.);
		std::vector<double> y2gt0 = use(y2, y2>0.);

		std::vector<double> z1gt0 = use(z1, z1>0.);
		std::vector<double> z2gt0 = use(z2, z2>0.);

		i1 = nm::argmax(z1gt0);
		i2 = nm::argmax(z2gt0);
		
		//dx = x2gt0[i2] - x1gt0[i1];
		//dy = y2gt0[i2] - y1gt0[i1];
		dz = z2gt0[i2] - z1gt0[i1];

		thick_pole = std::sqrt(dz*dz);
		
		fprintf(ofp_thick, "%lf %lf\n", thick_side, thick_pole);
	}

	fclose(ofp_thick);

	return 0;
}
