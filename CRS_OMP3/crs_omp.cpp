#include "crs_omp.h"
#include <iostream>
#include <random>
#include <vector>
#include <math.h>
#include <omp.h>

double crs_omp::rand_val()
{
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 eng(rd()); // seed the generator
	std::uniform_real_distribution<> distr(-D_max, D_max); // define the range

	return distr(eng);
}

double crs_omp::evaluate_F(std::vector<double> p)
{
	double val_sum = 0, val_mul = 1;
	int size = p.size();
	int i;
//#pragma omp parallel for reduction (+:val_sum) reduction (*:val_mul)
	for (i = 0; i < size; i++)
	{
		val_sum += p[i] * p[i];
		val_mul *= cos(p[i] / ((double)i + 1));
	}

	return (val_sum / 40 + 1 - val_mul);
};

bool crs_omp::is_in_domain(std::vector<double> r)
{
	bool val = true;
	for (int i = 0; i < n; i++)
	{
		if (r[i] < -D_max || r[i] > D_max)
		{
			val = false;
			break;
		}
	}

	return val;
};

void crs_omp::init(int dim)
{
	double wtime;
	double wtime1;
	double wtime2;
	wtime1 = omp_get_wtime();
	// assign n
	this->n = dim;
	// take N sampled points in the domain
	N = multiplier * (n + 1);

	// initialize set sampled points P
	int i, j;
	P = std::vector<std::vector<double>>(N);
	//int no_threads = omp_get_num_threads() * 2;
#pragma omp parallel private(i,j)
	{
		unsigned int myseed = omp_get_thread_num();

	#pragma omp for
		for (i = 0; i < N; i++) {
			std::vector<double> p = std::vector<double>(n);
			for (j = 0; j < n; j++) {
				p[j] = rand_val();
			}
			P[i] = p;
		}
	}

	// compute evaluation values of the sample points and find the best, worst values
	// add the function value of the first point to F
	F = std::vector<double>(N);
	F[0] = evaluate_F(P[0]);
	// assign initialized value to F_min, F_max
	F_min = F_max = F[0];
	idx_min = idx_max = 0;

#pragma	omp parallel
	{
		double F_min1 = F[0], F_max1 = F[0];
		int idx_min1 = 0, idx_max1 = 0;

#pragma omp for private(i) nowait 
		for (int i = 1; i < N; i++)
		{
			double val = evaluate_F(P[i]);
			F[i] = val;

			if (val > F_max1) {
				F_max1 = val;
				idx_max1 = i;
			}
			if (val < F_min1) {
				F_min1 = val;
				idx_min1 = i;
			}
		}

	#pragma omp critical 
		if (F_max1 > F_max) {
			F_max = F_max1;
			idx_max = idx_max1;

		}
	#pragma omp critical 
		if (F_min1 < F_min) {
			F_min = F_min1;
			idx_min = idx_min1;
		}
	}
	
	//numb_of_evals = N;

	wtime2 = omp_get_wtime();
	wtime = wtime2 - wtime1;
	std::cout << "\nInitialization time = " << wtime << "\n\n";
};


void crs_omp::crs_exec()
{
	// check stopping condition
	if (F_max - F_min < eps)
	{
		print_solution();
		return;
	}

	int buffer_size = 4;
	int LOOP_MAX = 100000/ buffer_size;
	int loop_count = 0;
	do
	{
		/*double wtime;
		double wtime1;
		double wtime2;
		wtime1 = omp_get_wtime();*/

		loop_count++;

		// init [n*no_threads] trial points set and new solution points
		std::vector<std::vector<std::vector<double>>> T = std::vector<std::vector<std::vector<double>>>(buffer_size);
		std::vector<std::vector<double>> r = std::vector<std::vector<double>>(buffer_size);
		std::vector<double> r_val = std::vector<double>(buffer_size);
		// set default values
	#pragma omp parallel for
		for (int t = 0;t < buffer_size;t++)
		{
			T[t] = std::vector<std::vector<double>>(n);
			r[t] = std::vector<double>(n);
			r_val[t] = -1;
		}

		// choose randomly from P [n*no_threads] trial points
	#pragma omp parallel for
		for (int t = 0; t < buffer_size; t++) {
			for (int idx = 0; idx < n; idx++) {
				T[t][idx] = P[rand() % N];
			}

			// add the best point to trial points (n -> n+1)
			T[t].insert(T[t].begin(), P[idx_min]);
		}		

		// compute its centroid, and new solution points
		int t, eval_count = 0;
	#pragma omp parallel for private (t) reduction(+: eval_count)
		for (t = 0; t < buffer_size; t++)
		{
			std::vector<double> c = std::vector<double>(n);

			for (int dim = 0; dim < n; dim++) // dimensions
			{
				double sum = 0;
				for (int idx = 0; idx < n; idx++) // points of T
				{
					sum += T[t][idx][dim];
				}
				c[dim] = sum / n;
			}

			// compute new solution point -> r = 2*c - T(n)
			for (int dim = 0; dim < n; dim++)
			{
				r[t][dim] = 2 * c[dim] - T[t][n][dim];
			}

			if (is_in_domain(r[t])) {
				// compute evaluation value at r
				r_val[t] = evaluate_F(r[t]);
				eval_count++;
			}
			else
			{
				// another way to compute new solution point -> r = (c + T(n))/2
				for (int dim = 0; dim < n; dim++)
				{
					r[t][dim] = (c[dim] + T[t][n][dim]) / 2;
				}

				if (is_in_domain(r[t]))
				{
					// compute evaluation value at r
					r_val[t] = evaluate_F(r[t]);
					eval_count++;
				}
			}
		}
		numb_of_evals += eval_count;

		for (int t = 0; t < buffer_size; t++)
		{
			if (r_val[t] == -1) continue;

			if (r_val[t] < F_max)
			{
				// update P, best and worst points
				// update P
				// replace the worst point and it's func value by 
				// new solution point [p] and it's func val [f_val]
				P.erase(P.begin() + idx_max);
				F.erase(F.begin() + idx_max);
				P.insert(P.begin() + idx_max, r[t]);
				F.insert(F.begin() + idx_max, r_val[t]);

				// update best point on P
				if (r_val[t] < F_min)
				{
					idx_min = idx_max;
					F_min = r_val[t];
				}

				// refind worst point on P
				F_max = r_val[t];
				for (int i = 0; i < N; i++)
				{
					if (F[i] > F_max)
					{
						F_max = F[i];
						idx_max = i;
					}
				}
			}
		}

	//} while (F_max - F_min >= eps);
	} while (F_max - F_min >= eps && loop_count < LOOP_MAX);
};

void crs_omp::print_solution()
{
	std::cout << "x=[";
	for (size_t i = 0; i < n-1; i++)
	{
		printf("%f,", P[idx_min][i]);
	}
	printf("%f]\n", P[idx_min][n-1]);
	std::cout << "Number of evaluations: " << numb_of_evals << "\n";
	std::cout << "Min value = " << F_min << "\n";
	std::cout << "Max value = " << F_max << "\n";
};
