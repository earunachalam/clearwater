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

	const double phi = 0.10 ;

	const double dt = 0.10 /2;
	//=============================================================================================================================

	//=============================================================================================================================
	//
	double* x;
	double* y;

	double* fx;
	double* fy;

	int i=0, j=0, k=0;

	double x1=0, y1=0, x2=0, y2=0;
	double r2=0, r=0;
	int idx1=0, idx2=0;

	double rx=0, ry=0;

	int NumPoints = 0;

	char str1[100];
	char str2[100];

	ifstream* rd1;
	ifstream* rd2;

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

	fx = new double[NumPoints];
	fy = new double[NumPoints];

	rd1 = new ifstream("./meshfiles/xs.dat");
	rd2 = new ifstream("./meshfiles/ys.dat");

	for(i=0; i<NumPoints; i++)
	{
		rd1->getline(str1,100);
		rd2->getline(str2,100);

		x[i] = atof(str1);
		y[i] = atof(str2);
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

		edg1[i] = int( atof(str1) )-1;
		edg2[i] = int( atof(str2) )-1;
	}
	rd1->close();
	delete rd1;

	rd2->close();
	delete rd2;

	double* EdgeLengths0 = new double[NumEdges];

	for(i=0; i<NumEdges; i++)
	{
		idx1 = edg1[i];
		idx2 = edg2[i];

		x1 = x[idx1];
		y1 = y[idx1];

		x2 = x[idx2];
		y2 = y[idx2];

		rx = x2-x1;
		ry = y2-y1;

		r = sqrt(rx*rx+ry*ry);

		EdgeLengths0[i] = r;
	}

	int* EdgeUnderTension = new int[NumEdges];

	for(i=0; i<NumEdges; i++)
	{
		EdgeUnderTension[i] = 0;

		idx1 = edg1[i];
		idx2 = edg2[i];

		x1 = x[idx1];
		y1 = y[idx1];

		x2 = x[idx2];
		y2 = y[idx2];

		if( x1>25.0 && x1<55.0 && y1>25.0 && y1<55.0 )
		{
			EdgeUnderTension[i] = 1;
		}
	}

	ofstream * file = new ofstream("out.m");
	(*file)<<"xs=[];"<<endl;
	file->close();
	delete file;
	//=============================================================================================================================

	for(i=0; ;i++)
	{
		if( !(i/1000.0-i/1000) ) { cout<<"Iteration # "<<i<<endl;
		}

		for(j=0; j<NumPoints; j++)
		{
			fx[j] = 0;
			fy[j] = 0;
		}

		for(j=0; j<NumEdges; j++)
		{

			idx1 = edg1[j];
			idx2 = edg2[j];

			L0 = EdgeLengths0[j];

			x1 = x[idx1];
			y1 = y[idx1];

			x2 = x[idx2];
			y2 = y[idx2];

			r2 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) ;

			r = sqrt(r2);

			rx = (x2-x1)/r;
			ry = (y2-y1)/r;

			/*
			double dL = r-L0 ;
			if( dL<0 ) dL *= -1.0;

			double KS_ = 50*KS / ( 1.0+dL*dL*5 ) ;
			*/

			fx[idx1] +=  KS* rx*(r-L0);
			fy[idx1] +=  KS* ry*(r-L0);

			fx[idx2] += -KS* rx*(r-L0);
			fy[idx2] += -KS* ry*(r-L0);

			if( EdgeUnderTension[j]==1 )
			{
				fx[idx1] +=  0.1*rx ;
				fy[idx1] +=  0.1*ry ;

				fx[idx2] += -0.1*rx ;
				fy[idx2] += -0.1*ry ;
			}

			/*
			fx[idx1] +=  mu* rx;
			fy[idx1] +=  mu* ry;

			fx[idx2] += -mu* rx;
			fy[idx2] += -mu* ry;
			*/
		}

		// if( (i<1.0 * 1*100000) || ((i> 1.0 * 20*100000) && (i< 1.0 * 21*100000)) )
		// if( i<1.0 * 1*100000/100 )
		/*
		if( i*dt<Toff )
		{
			// fx[3630-1] += phi ; // 15.0 ;

			fx[4020-1] += phi ;
		}
		*/

		// fx[4020-1] += phi ;

		// fx[16580-1] += phi ;

		// fx[5754-1] += phi ;

		for(j=0; j<NumPoints; j++)
		{
			/*
			if( x[j]<1.0 ) continue;
			if( x[j]>102.0 ) continue;

			if( y[j]<1.0 ) continue;
			if( y[j]>88.0 ) continue;
			*/

			x[j] += fx[j] *dt ;
			y[j] += fy[j] *dt ;
		}

		#include "outpt.h"
	}

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
