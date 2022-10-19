// Author Cole Harlow
// Project: HW2 ECE 7011

//Compile with
//Note the boost and gnuplot are dependencies
//g++ HW2_Q1.cpp -lboost_iostreams -lboost_system -lboost_filesystem

#include <iostream>
#include <iomanip>
#include <cmath>
#include <boost/tuple/tuple.hpp>
#include "gnuplot-iostream.h"

#define NumEFieldPoints 501
#define NumHFieldPoints 500
#define NumSteps 700 //The number of steps 
#define initial_Cond 1 //Determines which initial conditions we use

//These parameters are used in the initial conditions and the exact solution
double Eps0 = 8.8541878128*pow(10,-12);
double Mu0 = 1.25663706*pow(10,-6); 

double Ex[NumEFieldPoints] = {0}; //The number of electric field points that there are
double Hy[NumHFieldPoints] = {0}; //The number of magnetic HALF field points that there are
double exactEx[NumEFieldPoints] = {0}; //same as previous two lines but this array holds the exact values
double exactHy[NumHFieldPoints] = {0};

double ic = 80;
double iw = 4;
double eta0 = sqrt(Mu0/Eps0);
double c = 1.0/(sqrt(Mu0*Eps0)); // Sets the speed of light
double Cf = 1.0/2.0; //Courant Number

//Used for calculation of exactEx values
double delZ = pow(10,-3); //sets the size of space between each node
double delt = 1.67*pow(10,-12); //Sets the time interval between iterations
double Zstart = 0; //sets the starting point on the Z-axis
double Zend = 0.5; //set the ending point on the Z-axis
double zc = ic*delZ;
double w = iw*delZ;

double totErr[NumSteps][2] = {0};

int main()
{
	Gnuplot gp; //Creates the interface for communicating with gnuplot
	Gnuplot gp2; //Creates a second interface for communicating with gnuplot

	//Creates an array which tells us when to sample the field solutions for the HW
	int take_sample[8] = {0,100,200,300,400,500,600,700};

	//Sets up 2D arrays whihc are used to hold the field values for each of the 8 samples
	double EfieldSamples[NumEFieldPoints+200][9] = {0};
	double exactEfieldSamples[NumEFieldPoints+200][9] = {0};

	//Initializes the first Column of each array so that it holds the grid point or half grid point number
	//This is used later for graphing purposes 
	for(int i = 0; i<NumEFieldPoints+200; i++)
	{
		EfieldSamples[i][0] = i*delZ;
		exactEfieldSamples[i][0] = i*delZ;
	}

	for(int i = 0; i<NumSteps;i++)
	{
		totErr[i][0] = i;
	}

	/*-----------------------------CALCULATES THE INITIAL CONDITIONS---------------*/
	for(int i = 0; i<NumEFieldPoints ; i++)
	{
		if(initial_Cond == 0)
		{
			Ex[i] = exp(-pow((i-ic)/iw,2));
		}
		else if(initial_Cond == 1)
		{
			Ex[i] = Ex[i] = exp(-pow((i-ic)/iw,2));
		}
	}
	//Calculates the initial condition for the magnetic fields
	for(int i = 0; i<NumHFieldPoints; i++)
	{
		if(initial_Cond == 0)
		{
			Hy[i] = 1/eta0 * exp(-pow(((i+1.0/2.0-ic-Cf/2.0)/iw),2));
		}
		else if(initial_Cond == 1)
		{
			Hy[i] = 0;
		}
	}

	/*-----------------------------THIS IS THE SECTION IN WHICH WE ITERATE----------*/
	//Used to keep track of what timestep we are on. 
	int timeStep = 0;
	//Holds the number of saved samples in the samples matrixs
	int numSaved = 0;

	//Starts looping through the desired number of time steps
	while(timeStep <= NumSteps)
	{
		/*-----------------------Subsection to calculate the theoretical value---*/
		double z = 0; //Used in part to calculate the fields
		double t = delt*timeStep; //Calculates what time we are currently at

		for(int i=0; i< NumEFieldPoints; i++)
		{
			z = delZ*i;	//Calculates the Z location based on the index i 
			exactEx[i] = exp(-pow((z-zc-c*t)/w,2)); //Calculates the electric field at z index i
		}

		/*-----------------------Subsection to store the desired iterations-----------*/
		if(timeStep == take_sample[numSaved])
		{
			std::cout << "Sample Saved \n";
			//Saves the FDTD calculated and exact E fields
			for(int i = 0; i< NumEFieldPoints; i++)
			{
				EfieldSamples[i][numSaved+1] = Ex[i];
				exactEfieldSamples[i][numSaved+1] = exactEx[i];
			}
			numSaved++;
		}

		/*-------------------Subsection to calculate the error--------------------*/
		for(int i=0; i < NumEFieldPoints; i++)
		{
			totErr[timeStep][1] = totErr[timeStep][1]+pow((exactEx[i] - Ex[i]),2);
		}

		totErr[timeStep][1] = totErr[timeStep][1]/NumSteps;

		std::cout<< timeStep << std::endl;
		timeStep++;

		/*-----------------------Subsection to calculate the next FDTD value ---------*/
		//Updates the Magnetic Field (Hy) to the next half time step
		//Enforces the dirichlet boundary conditions
		Ex[0] = 0;
		Ex[500] = 0;

		for(int i = 0; i < NumHFieldPoints; i++)
		{
			//Iterates to the next half time step for the magnetic field
			//Hy[0] in reality is Hy[1/2] in the same way Hy[1] is Hy[1+1/2] etc
			Hy[i] = Hy[i]-(delt/(Mu0*delZ))*(Ex[i+1]-Ex[i]);
		}
		//Updates the Electric Field (Ex) to the next half time step
		for(int i = 1; i < NumEFieldPoints-1; i++)
		{
			//Updates the electric field 
			Ex[i]= Ex[i]-(delt/(Eps0*delZ))*(Hy[i]-Hy[i-1]);
		}
	}

	std::cout<<"Done \n";
	std::cout<<"Begin Plotting ... \n";

	//Turns on the Autoscale
	gp << "set term qt 0 font 'Helvetica,20'\n";
	gp << "set xlabel 'location (m)'\n";
	gp << "set ylabel 'E field (V/m)'\n";
	gp << "set autoscale\n";
	gp << "plot '-' u 1:2 w l title 'FDTD 0','-' u 1:2 w l title 'Exact 0',\
		  '-' u 1:3 w l title 'FDTD 100','-' u 1:3 w l title 'Exact 100',\
		  '-' u 1:4 w l title 'FDTD 200','-' u 1:4 w l title 'Exact 200',\
		  '-' u 1:5 w l title 'FDTD 300','-' u 1:5 w l title 'Exact 300'\n";
	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);
	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);
	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);
	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);
	

	gp << "set term qt 1 font 'Helvetica,20'\n";
	gp << "set xlabel 'location (m)'\n";
	gp << "set ylabel 'E field (V/m)'\n";
	gp << "set autoscale\n";
	gp << "plot '-' u 1:6 w l title 'FDTD 400','-' u 1:6 w l title 'Exact 400',\
		  '-' u 1:7 w l title 'FDTD 500','-' u 1:7 w l title 'Exact 500',\
	  	  '-' u 1:8 w l title 'FDTD 600','-' u 1:8 w l title 'Exact 600',\
	      '-' u 1:9 w l title 'FDTD 700','-' u 1:9 w l title 'Exact 700'\n";
 	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);
	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);
	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);
	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);

	gp << "set term qt 2 font 'Helvetica,20'\n";
	gp << "set xlabel 'Iteration'\n";
	gp << "set ylabel 'Error'\n";
	gp << "set autoscale\n";
	gp << "plot '-' u 1:2 w l title 'Error At Iteration'\n";
	gp.send1d(totErr);

/*
	gp << "set term qt 1 font 'Helvetica,20'\n";
	gp << "set autoscale\n";
	gp << "plot '-' u 1:3 w l title 'FDTD 100','-' u 1:3 w l title 'Exact 100'\n";
	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);


	gp << "set term qt 2\n";
	gp << "set autoscale\n";
	gp << "plot '-' u 1:4 w l title 'FDTD 200','-' u 1:4 w l title 'Exact 200'\n";
	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);

	gp << "set term qt 3\n";
	gp << "set autoscale\n";
	gp << "plot '-' u 1:5 w l title 'FDTD 300','-' u 1:5 w l title 'Exact 300'\n";
	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);

	gp << "set term qt 4\n";
	gp << "set autoscale\n";
	gp << "plot '-' u 1:6 w l title 'FDTD 400','-' u 1:6 w l title 'Exact 400'\n";
	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);

	gp << "set term qt 5\n";
	gp << "set autoscale\n";
	gp << "plot '-' u 1:7 w l title 'FDTD 500','-' u 1:7 w l title 'Exact 500'\n";
	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);

	gp << "set term qt 6\n";
	gp << "set autoscale\n";
	gp << "plot '-' u 1:8 w l title 'FDTD 600','-' u 1:8 w l title 'Exact 600'\n";
	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);

	gp << "set term qt 7\n";
	gp << "set autoscale\n";
	gp << "plot '-' u 1:9 w l title 'FDTD 700','-' u 1:9 w l title 'Exact 700'\n";
	gp.send1d(EfieldSamples);
	gp.send1d(exactEfieldSamples);

	gp << "set term qt 8\n";
	gp << "set autoscale\n";
	gp << "plot '-' u 1:2 w l title 'Error At Iteration'\n";
	gp.send1d(totErr);
*/
}