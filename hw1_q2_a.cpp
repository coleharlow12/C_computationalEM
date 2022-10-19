// Compile it with:
// g++ hw1_q2_a.cpp -lboost_iostreams -lboost_system -lboost_filesystem

#include<iostream>
#include <iomanip>
#include <cmath>
#include "gnuplot-iostream.h"

#define runs 200
float r = 0.1;

float Phi_n[11] = {0};
float Phi_n_p_1[11] = {0};
double Data[11][5];

int main()
{
	Gnuplot gp;
	//Creates a vector for storing the iteration number and the iteration error. Used to plot later
	std::vector<boost::tuple<int,double,double,double,double> > location;

	int take_sample[4] = {5,50,80,140};
	int j = 0;

	/*------------------------------SETS UP INITIAL CONDITIONS-------------------*/
	Phi_n[5] = 10;

	//Fills up the first column of the matrix with the desired data
	for(int i=0;i<=11;i++)
		{Data[i][0]=i;}


	/*------------------------------BEGINS CALCULATIONS--------------------------*/
	for(int k = 0; k < runs; k++)
	{
		for(int i = 1; i < 11; i++)
		{
			//Calculates the next value in time of the simulation. 
			Phi_n_p_1[i] = Phi_n[i]-r*(Phi_n[i+1]-Phi_n[i-1]);
		}

		for(int i = 1; i < 11; i++)
		{
			Phi_n[i] = Phi_n_p_1[i];
		}

		//Determines if this is one of the iterations to save and grap based on take_sample
		if(k == take_sample[j])
		{
			j++; //Increments j 
			for(int l=0;l<11;l++)
			{
				Data[l][j]=Phi_n[l];
			}
		}
	}

	for(int i = 0; i<=11;i++)
	{
		location.push_back(boost::make_tuple(Data[i][0],Data[i][1],Data[i][2],Data[i][3],Data[i][4]));
	}


	gp << "set autoscale\n";
	gp << "plot '-' u 1:2 w l title 'First','-' u 1:3 w l title 'Second','-' u 1:4 w l title 'Third','-' u 1:5 w l title 'Fourth'\n";
	gp.send1d(location);
	gp.send1d(location);
	gp.send1d(location);
	gp.send1d(location);
}