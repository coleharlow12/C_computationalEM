// Author Cole Harlow
// Project: HW2 ECE 7011

//Compile with
//Note the boost and gnuplot are dependencies
//g++ HW2_Q2_V2.cpp -lboost_iostreams -lboost_system -lboost_filesystem -lfftw3 -lm

#include <iostream>
#include <iomanip>
#include <cmath>
#include <boost/tuple/tuple.hpp>
#include "gnuplot-iostream.h" //Must be linked with -lboost_iostreams -lboost_system -lboost_filesystem
#include <fftw3.h> //Must be linked with -lfftw3 -lm

//Sets the number of nodes
#define numNodesX 101
#define numNodesY 101
#define numTimeSteps 8192 //Must be an even number
#define Cf 0.8

//These parameters are used in the initial conditions and the exact solution
double Eps0 = 8.8541878128*pow(10,-12);
double Mu0 = 1.25663706*pow(10,-6); 
double C = 1.0/(sqrt(Mu0*Eps0)); // Sets the speed of light

//Grid Spacings
double delX = 50*pow(10,-3)*pow(10,-2);
double delY = 100*pow(10,-3)*pow(10,-2);
double delt = delX*delY*Cf/(delY*C+delX*C);

//Sets the ranges of the X and Y coordinates of the 
double Xstart = 0;
double Xend = 100*pow(10,-3);
double Ystart = 0;
double Yend = 50*pow(10,-3);

//Creates the H and E vectors
//Remember the H fields are on the half grid points so there is one less row/column
//Entries with a common value Y are described as ROWS
//Entries with a common value of X are described as COLUMNS
double Hx[(numNodesX)*(numNodesY)] = {0}; //For Simplicity Later I keep the number of Rows and columns the same in Hx and Hy. Imagine Hx[0] is a half grid point in the y direction above Ez[0]
double Hy[(numNodesX)*(numNodesY)] = {0}; //Imagine Hy[0] is a half grid point in the x direction to the right of Ez[0]
double Ez[numNodesX*numNodesY] = {0};
double E_at_13_29[numTimeSteps][2] = {0};
double MagFFT[numTimeSteps/2+1][2] = {0};

int numNodes = numNodesX*numNodesY;

int main()
{
	fftw_complex *in, *out;
	fftw_plan p;

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numTimeSteps);		//Allocates a block of memory size to hold the FFT input data
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * numTimeSteps);		//Allocates a block of memory size to hold the FFT output data
	p = fftw_plan_dft_1d(numTimeSteps, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	Gnuplot gp; //Used to plot a figure using GNUplot

	//This for loop initializes the first row of the array E_at_13_29
	for(int i = 0; i<numTimeSteps;i++)
	{
		E_at_13_29[i][0] = i;
	}
	for(int i = 0; i<numTimeSteps/2+1;i++)
	{
		MagFFT[i][0] = i*(1/delt) / numTimeSteps;
	}

	/*----------------SETS THE INITIAL CONDITIONS FOR FIELDS--------------*/
	srand(time(NULL)); //sets the seed of the random number generator
	//The field values are random across the entire grid.
	for(int i=0;i<numNodes;i++) {Ez[i] = (rand()%100)+100;}

	/*--------SETS THE BOUNDARY CONDITIONS FO THE CONDUCTING WALLS----------*/
	for(int i=0; i<numNodesY;i++)
	{
		//Sets the electric field at the left conductive boundary to zero
		Ez[numNodesX*i] = 0;
		//Sets the electric field at the right conductive boundary to zero
		Ez[(numNodesX-1)+numNodesX*i] = 0;
	}

	for(int i=0; i<numNodesX;i++)
	{
		//Sets the electric field at the bottom conductive boundary to zero
		Ez[i] = 0;
		//Sets the electric field at the top conductive boundary to zero
		Ez[i+(numNodesX)*(numNodesY-1)] = 0;
	}	

	/*---------------------------------------------------------------------------BEGINS ITERATIONS TO CALCULATE THE FIELD DIST-------------------------------------------*/
	for(int t = 0; t<numTimeSteps; t++)
	{
		//Stores the value of the field at grid point (13,29) or (i = 29*numNodesX + 13) = 2942
		E_at_13_29[t][1] = Ez[2942];
		in[t][0] = Ez[2942]; //Fills FFT array with the field information
		in[t][1] = 0; 

		//First we start by updating the Y component of the magnetic field
		for(int i = 0; i < numNodesX*numNodesY; i++)
		{
			Hy[i] = Hy[i] + delt/(Mu0*delX) * (Ez[i+1]-Ez[i]);
		}
		//This section updates all of the Hx values this is easy since Hx just has one less row then Ex
		for(int i = 0; i < numNodesX*numNodesY; i++)
		{
			Hx[i] = Hx[i] + delt/(Mu0*delY) * (Ez[i]-Ez[i+numNodesX]);
		}

		//Finally we update the Electric Field at all points on the grid. Note we skip the first and last rows of the grid. They are conductive B.Cs 
		for(int i = numNodesX; i< numNodesX*(numNodesY-1); i++)
		{
			Ez[i] = Ez[i] + delt/(Eps0*delX)*(Hy[i]-Hy[i-1]) - delt/(Eps0*delY)*(Hx[i]-Hx[i-numNodesX]);
		}

		//Resets the left and right wall boundary conditions for the next iteration
		for(int i=0; i<numNodesY;i++)
		{
		//Sets the electric field at the left conductive boundary to zero
		Ez[numNodesX*i] = 0;
		//Sets the electric field at the right conductive boundary to zero
		Ez[(numNodesX-1)+numNodesX*i] = 0;
		}

	}

	fftw_execute(p);

	for(int i = 0; i<numTimeSteps/2+1;i++)
	{
		//std::cout << "The time index is " << E_at_13_29[i][0] << " The values are " << E_at_13_29[i][1] << std::endl;
		std::cout << "Real Part: " <<  out[i][0] << " Imag Part: " << out[i][1] << std::endl;
		MagFFT[i][1] = sqrt(pow(out[i][0],2)+pow(out[i][1],2));
	}

	//Turns on the Autoscale
	gp << "set term qt 0 font 'Helvetica,20'\n";
	gp << "set xrange [0:10000000000]\n";
	gp << "set xlabel 'Frequency (Hz)'\n";
	gp << "set ylabel 'Magnitude'\n";
	gp << "set autoscale y\n";
	gp << "plot '-' u 1:2 w l title 'FDTD'\n";
	gp.send1d(MagFFT);
	gp << "\n";
}
