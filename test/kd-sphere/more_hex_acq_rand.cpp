#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
using namespace std ;

char* itoa(int val, int base) ;

int main()
{
	//=============================================================================================================================
	//
	const double KS = 0.10 ;

	/*const*/ double L0 = 1.0 ;

	// const double Toff = 1000;

	//const double phi = 0.10 ;

	const double dt = 0.005;
	//=============================================================================================================================

	//=============================================================================================================================
	//
	double* x;
	double* y;
	double* z;

	double* fx;
	double* fy;
	double* fz;

	int i=0, j=0;

	double x1=0, y1=0, z1=0, x2=0, y2=0, z2=0;
	double r2=0, r=0;
	int idx1=0, idx2=0;

	double rx=0, ry=0, rz=0;

	int NumPoints = 0;

	char str1[100];
	char str2[100];
	char str3[100];

	ifstream* rd1;
	ifstream* rd2;
	ifstream* rd3;

	rd1 = new ifstream("./meshfiles/xs.dat");

	while( rd1->getline(str1,100) )
	{
		NumPoints++;
	}
	rd1->close();
	delete rd1;

	// cout<<NumPoints<<endl;

	x = new double[NumPoints];
	y = new double[NumPoints];
	z = new double[NumPoints];

	fx = new double[NumPoints];
	fy = new double[NumPoints];
	fz = new double[NumPoints];

	rd1 = new ifstream("./meshfiles/xs.dat");
	rd2 = new ifstream("./meshfiles/ys.dat");
	rd3 = new ifstream("./meshfiles/zs.dat");

	for(i=0; i<NumPoints; i++)
	{
		rd1->getline(str1,100);
		rd2->getline(str2,100);
		rd3->getline(str3,100);

		x[i] = atof(str1);
		y[i] = atof(str2);
		z[i] = atof(str3);
	}
	rd1->close();
	delete rd1;

	rd2->close();
	delete rd2;

	int NumEdges = 0;
	rd1 = new ifstream("./meshfiles/edg1.dat");

	while( rd1->getline(str1,100) )
	{
		NumEdges++;
	}
	rd1->close();
	delete rd1;

	int* edg1 = new int[NumEdges];
	int* edg2 = new int[NumEdges];

	rd1 = new ifstream("./meshfiles/edg1.dat");
	rd2 = new ifstream("./meshfiles/edg2.dat");

	for(i=0; i<NumEdges; i++)
	{
		rd1->getline(str1,100);
		rd2->getline(str2,100);

		edg1[i] = int( atof(str1) );
		edg2[i] = int( atof(str2) );
	}
	rd1->close();
	delete rd1;

	rd2->close();
	delete rd2;

	rd3->close();
	delete rd3;

	double* EdgeLengths0 = new double[NumEdges];

	for(i=0; i<NumEdges; i++)
	{
		idx1 = edg1[i];
		idx2 = edg2[i];

		x1 = x[idx1];
		y1 = y[idx1];
		z1 = z[idx1];

		x2 = x[idx2];
		y2 = y[idx2];
		z2 = z[idx2];

		rx = x2-x1;
		ry = y2-y1;
		rz = z2-z1;

		r = sqrt(rx*rx+ry*ry+rz*rz);

		EdgeLengths0[i] = r;
	}

	for(i=0; i<NumEdges; i++)
	{

		idx1 = edg1[i];
		idx2 = edg2[i];

		x1 = x[idx1];
		y1 = y[idx1];
		z1 = z[idx1];

		x2 = x[idx2];
		y2 = y[idx2];
		z2 = z[idx2];

	}

	//ofstream * file = new ofstream("out.m");
	//(*file)<<"xs=[];"<<endl;
	//file->close();
	//delete file;
	//========================================================================================================

	ofstream * file = new ofstream("traj.lammpstrj");
    (*file) << "ITEM: TIMESTEP\n" << 0 << "\n";
    (*file) << "ITEM: NUMBER OF ATOMS\n" << NumPoints << "\n";
    (*file) << "ITEM: BOX BOUNDS pp pp pp\n-100.0 100.0\n-100.0 100.0\n-100.0 100.0\n";
    (*file) << "ITEM: ATOMS id type x y z\n";
    for(j=0; j<NumPoints; j++)
    {
        (*file) << j
                << " 1 " << 
           x[j] << " " <<
           y[j] << " " <<
           z[j] << "\n";
    }

	for(i=1; ;i++)
	{
        if (i%1000 == 0)
            cout<<i<<"\n";
		//if( !(i/1000.0-i/1000) ) { cout<<"Iteration # "<<i<<endl;
		//}

		for(j=0; j<NumPoints; j++)
		{
			fx[j] = 0;
			fy[j] = 0;
			fz[j] = 0;
		}

		for(j=0; j<NumEdges; j++)
		{

			idx1 = edg1[j];
			idx2 = edg2[j];

			L0 = EdgeLengths0[j];

			x1 = x[idx1];
			y1 = y[idx1];
			z1 = z[idx1];

			x2 = x[idx2];
			y2 = y[idx2];
			z2 = z[idx2];

			r2 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1);

			r = sqrt(r2);
            //printf("%.6f ", r);

			rx = (x2-x1)/r;
			ry = (y2-y1)/r;
			rz = (z2-z1)/r;
            
            //if (idx1==0 || idx2==0)
            //{
                //printf("%d %d %d\n", j, idx1, idx2); 
                //printf("%f %f %f\n", x1, y1, z1);
                //printf("%f %f %f\n", x2, y2, z2);
                //printf("%f %f %f\n", rx, ry, rz);
                //printf("%f %f\n\n", r2, r); 
                ////abort();
            //}

			fx[idx1] +=  KS* rx*(r-L0);
			fy[idx1] +=  KS* ry*(r-L0);
			fz[idx1] +=  KS* rz*(r-L0);

			fx[idx2] += -KS* rx*(r-L0);
			fy[idx2] += -KS* ry*(r-L0);
			fz[idx2] += -KS* rz*(r-L0);

            fx[idx1] +=  0.1*rx ;
            fy[idx1] +=  0.1*ry ;
            fz[idx1] +=  0.1*rz ;

            fx[idx2] += -0.1*rx ;
            fy[idx2] += -0.1*ry ;
            fz[idx2] += -0.1*rz ;
		}
        //abort();
        //cout << "\n";

        //if (i == 2)
            //for (j = 0; j < NumPoints; ++j)
                //printf("%f %f %f %f\n", x[j], y[j], fx[j], fy[j]);

		for(j=0; j<NumPoints; j++)
		{
			x[j] += fx[j] *dt ;
			y[j] += fy[j] *dt ;
			z[j] += fz[j] *dt ;

            //if (j == 0)
            //{
                //printf("%f %f %f\n", fx[j], fy[j], fz[j]);
                //abort();
            //}
		}

        //#include "outpt.h"
        if (i%1 == 0)
        {
            (*file) << "ITEM: TIMESTEP\n" << i << "\n";
            (*file) << "ITEM: NUMBER OF ATOMS\n" << NumPoints << "\n";
            (*file) << "ITEM: BOX BOUNDS pp pp pp\n-100.0 100.0\n-100.0 100.0\n-100.0 100.0\n";
            (*file) << "ITEM: ATOMS id type x y z\n";
            for(j=0; j<NumPoints; j++)
            {
                (*file) << j
                        << " 1 " << 
                   x[j] << " " <<
                   y[j] << " " <<
                   z[j] << "\n";
            }
        }
	}

    file->close();

	return 1;
}

char* itoa(int val, int base) {
     static char buf[32] = {0};
     int i = 30;
     	if(val==0)
     	{
     	buf[i]='0';
	return &buf[i];
	}
     for (; val && i; --i, val /= base)
         buf[i] = "0123456789abcdef"[val % base];
     return &buf[i+1];
}
