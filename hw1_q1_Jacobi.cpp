// Author Cole Harlow
// Project: HW1 ECE 7011 PartA Jacobi

//Compile with
//Note the boost and gnuplot are dependencies
//g++ hw1_q1_Jacobi.cpp -lboost_iostreams -lboost_system -lboost_filesystem

#include <iostream>
#include <iomanip>
#include <cmath>
#include <boost/tuple/tuple.hpp>
#include "gnuplot-iostream.h"

#define N 30 //Number of points in the grid in the X direction
#define M 30 //Number of points in the grid in the Y direction
#define numBCX 2 //Number of boundary conditions in X direction
#define numBCY 2 //Number of boundary conditions in the Y directions
#define runs 500 //Number of Iterations
#define PI 3.14159265 //Defines Pi for 
#define non_Zero_Guess

float Phi_I[(N+numBCX)*(M+numBCY)] = {0}; //Creates an "array" to hold all of the information including the b.cs
float TempPhi[(N+numBCX)*(M+numBCY)] = {0}; //Creates a temporary "array" to hold the information during calculations
float Theoretical[(N+numBCX)*(M+numBCY)] = {0}; //Creates an "array" to hold the theoretically calculated values

float error_At_Point[M*N] = {0}; //Holds the error at each individual point during an iteration
float totErr[runs] = {0}; //Creates an empty array to hold the error of each iteration

int dimX = N+numBCY; //Total number of grid points on the X axis I assume the boundary values don't count in the grid
int dimY = M+numBCX; //Total number of grid points on the Y axis ''

int main()
{
	Gnuplot gp; //Creates the interface for communicating with gnuplot
	std::vector<std::pair<double, double> > Iteration_Error; //Creates a vector to store the error information in
	double ind_Point_Error[N*M][5]; //Creates a 2d array for storing the error at each point for a select number of iterations.
	int take_sample[4] = {5,50,80,140}; //Sets the iteration numbers at which we want to examine individual errors. 
	int j = 0;

	//Initializes the first row of the array so that it holds the grid point the error is reported on
	for(int i=0; i<M*N; i++)
		{ind_Point_Error[i][0]=i;}

	/*---------------------------------------CALCULATES THE THEORETICAL MATRIX--------------------------------------*/
	for (int i = 0; i<dimY*dimX; i++)
	{
		 //Note the backlashes(\) just allow the 3 following lines to be treated as a single line
		Theoretical[i] = sin(PI/(dimX-1)*(i%dimX))*\
						 sinh(PI/(dimY-1)*int(i/dimX))/\ 
						 sinh(PI); 
	}

	/*--------------------------------------SETS THE INITIAL GUESS OF THE SOLUTION-------------------------------*/
	//This ensures that the first and last rows are not considered since they are B.Cs
	#ifdef non_Zero_Guess
		for(int i = dimX+1; i<dimX*(dimY-1); i++)
		{
			if(i % dimX != 0 && i % dimX !=1) //This if statement ensures that the point is not at the end or beginning of the row (i.e exlcudes the first and last columns since they are b.c)
			{
				//The row y can be determined by dividing i by dimX and type casting to integer (this always rounds down)
				Phi_I[i] = sinh(PI/(dimY-1)*int(i/dimX))/sinh(PI);
			}
		}
	#endif

	/*---------------------------------------SETS UP THE BOUNDARY CONDITIONS----------------------------------------*/
	//Note that since at the x boundaries the values are just zero we do not need to calculate these B.Cs since this is how the vector was intiialized.

	//This section calculates the values at the Y boundarys 
	for(int i = dimY*(dimX-1); i<(dimY*dimX); i++)
	{
		Phi_I[i] = sin(PI*(i-dimY*(dimX-1))/(dimX-1)); //Note that we divide by (dimX-1) and not just dimX this is because there are 11 pointsso to have the 11th point occur at 0 we need 10 increments
	}

	/*---------------------------------------BEGINS ITERATING TO FIND AN ANSWER-----------------------------------*/
	for(int k = 0; k < runs; k++)
	{
		int l = 0; //Used to store information in the ind_Point_Error array
		
		//This ensures that the first and last column are not considered since they are B.Cs
		for(int i = dimX+1; i<dimX*(dimY-1); i++)
		{
			if(i % dimX != 0 && i % dimX !=1)
			{
				TempPhi[i] = 1/float(4)*(Phi_I[i-1]+Phi_I[i+1]+Phi_I[i+dimX]+Phi_I[i-dimX]); //Uses the difference equation to calculate the next guess of phi
				error_At_Point[l] = pow((Phi_I[i]-Theoretical[i]),2)/(M*N); //Calculates the error at each point
				l++;
			}
			else
			{
				TempPhi[i] = Phi_I[i];
			}
		}

		//Writes the new values of Phi into the Phi vector
		for (int i = dimX+1; i<dimX*(dimY-1);i++)
		{
			Phi_I[i] = TempPhi[i];
		}

			//Stores the data if this is one of the iterations of interest. 
		if(k == take_sample[j])
		{
			j++; //Incremets J note the first row holds the grid point number already
			for(int is=0;is<M*N;is++)
			{
				ind_Point_Error[is][j] = error_At_Point[is];
			}
		}

		//Calculates the error at points not of interest
		for(int is=0;is<M*N;is++)
			{
				totErr[k] = totErr[k] + error_At_Point[is];
			}
		Iteration_Error.push_back(std::make_pair(k,totErr[k]));
		//std::cout << totErr[k] << std::endl;;
	}

	gp << "set multiplot layout 1, 2\n";
	gp << "set title 'Global Error vs Iteration'\n";
	gp << "set autoscale\n";
	gp << "set grid\n";
	gp << "plot '-' with lines title 'Iteration Error'\n";
	gp.send1d(Iteration_Error);

	gp << "set title 'Local Error'\n";
	gp << "set autoscale\n";
	gp << "set grid\n";
	gp << "plot '-' u 1:2 w l title 'i = 5', '-' u 1:3 w l title 'i = 50', '-' u 1:4 w l title 'i = 80', '-' u 1:5 w l title 'i = 150'\n";
	gp.send1d(ind_Point_Error);
	gp.send1d(ind_Point_Error);
	gp.send1d(ind_Point_Error);
	gp.send1d(ind_Point_Error);
	gp << "unset multiplot";

return 0;
}