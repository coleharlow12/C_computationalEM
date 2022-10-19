// Compile it with:
// g++ hw1_q2_b.cpp -lboost_iostreams -lboost_system -lboost_filesystem

#include <iostream>
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
			Phi_n_p_1[i] = 0.5*(Phi_n[i+1]+Phi_n[i-1])-r*(Phi_n[i+1]-Phi_n[i-1]);
		}

		for(int i = 1; i < 11; i++)
		{
			Phi_n[i] = Phi_n_p_1[i];
		}

		if(k == take_sample[j])
		{
			j++;
			for(int l=0;l<11;l++)
			{
				Data[l][j]=Phi_n[l];
			}
		}
	}

	gp << "set autoscale\n";
	gp << "plot '-' u 1:2 w l title 'First','-' u 1:3 w l title 'Second','-' u 1:4 w l title 'Third','-' u 1:5 w l title 'Fourth'\n";
	gp.send1d(Data);
	gp.send1d(Data);
	gp.send1d(Data);
	gp.send1d(Data);
}
