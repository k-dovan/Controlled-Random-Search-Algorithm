#include "crs.h";
#include "crs_more_computation.h";
#include <random>
#include <iostream>
#include <omp.h>;
using namespace std;

int main(int argc, char** argv)
{
	if (argc == 2)
	{
		// get parameter n
		int n = atoi(argv[1]);

		// init parameters
		//crs c = crs();
		crs_more_computation c = crs_more_computation();

		double wtime;
		double wtime1;
		double wtime2;
		wtime1 = omp_get_wtime();

		c.init(n);

		c.crs_exec();

		wtime2 = omp_get_wtime();
		wtime = wtime2 - wtime1;

		std::cout << "************** SOLUTION *******************\n";
		cout << "Execution time = " << wtime << "\n";
		c.print_solution();

		return 0;
	}
	else return -1;	
};