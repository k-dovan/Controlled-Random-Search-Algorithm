#include "crs.h";
#include <iostream>;
#include <random>;
#include <vector>;
#include <math.h>;
#include <omp.h>;
using namespace std;

double crs::rand_val()
{
	random_device rd; // obtain a random number from hardware
	mt19937 eng(rd()); // seed the generator
	uniform_real_distribution<> distr(-D_max, D_max); // define the range

	return distr(eng);
};

double crs::evaluate_F(vector<double> p)
{
	double val_sum = 0, val_mul = 1;
	int size = p.size();
	for (int i = 0; i < size; i++)
	{
		val_sum += p.at(i) * p.at(i);
		val_mul *= cos(p.at(i) / ((double)i + 1));
	}

	return (val_sum / 40 + 1 - val_mul);
};

bool crs::is_in_domain(std::vector<double> r)
{
	bool val = true;
	for (int i = 0; i < n; i++)
	{
		if (r.at(i) < -D_max || r.at(i) > D_max)
		{
			val = false;
			break;
		}
	}

	return val;
};

void crs::init(int n)
{
	double wtime;
	double wtime1;
	double wtime2;
	wtime1 = omp_get_wtime();
	// assign n
	this->n = n;
	// take N sampled points in the domain
	N = multiplier * (n + 1);

	// initialize set sampled points P
	for (int i = 0; i < N; i++) {
		vector<double> p;
		for (int j = 0; j < n; j++) {
			p.push_back(rand_val());
		}
		P.push_back(p);
	}

	// compute evaluation values of the sample points and find the best, worst values
	// add the function value of the first point to F
	F.push_back(evaluate_F(P.at(0)));
	// assign initialized value to F_min, F_max
	F_min = F_max = F.at(0);
	idx_min = idx_max = 0;

	for (int i = 1; i < N; i++)
	{
		double val = evaluate_F(P.at(i));
		F.push_back(val);
		
		if (val > F_max) { 
			F_max = val; 
			idx_max = i;
		}
		if (val < F_min) {
			F_min = val;
			idx_min = i;
		}
	}

	//numb_of_evals = N;

	wtime2 = omp_get_wtime();
	wtime = wtime2 - wtime1;
	std::cout << "\nInitialization time = " << wtime << "\n\n";
};

void crs::update_P_best_worst_point()
{
	// update P
	// replace the worst point and it's func value by 
	// new solution point [r] and it's func val [r_val]
	P.erase(P.begin() + idx_max);
	F.erase(F.begin() + idx_max);
	P.insert(P.begin() + idx_max, r);
	F.insert(F.begin() + idx_max, r_val);

	// update best point on P
	if (r_val < F_min)
	{
		idx_min = idx_max;
		F_min = r_val;
	}

	// refind worst point on P
	F_max = r_val;
	for (int i = 0; i < N; i++)
	{
		if (F.at(i) > F_max)
		{
			F_max = F.at(i);
			idx_max = i;
		}
	}	
};

void crs::crs_exec()
{
	// check stopping condition
	if (F_max - F_min < eps)
	{
		print_solution();
		return;
	}

	int LOOP_MAX = 100000;
	int loop_count = 0;
	do
	{
		/*double wtime;
		double wtime1;
		double wtime2;
		wtime1 = omp_get_wtime();*/

		loop_count++;

		// choose randomly from P n trial points (set T-indexes referenced to P)
		T.clear();
		for (int i = 0; i < n; i++)
			T.push_back(rand() % N);
		
		// add the best point to trial points (n -> n+1)
		T.insert(T.begin(), idx_min);

		// compute its centroid except the last point of T (T(n))
		vector<double> c;
		for (int dim = 0; dim < n; dim++) // dimensions
		{
			double sum = 0;
			for (int idx = 0; idx < n; idx++) // points of T
			{
				sum += P.at(T.at(idx)).at(dim);
			}
			c.push_back(sum / n);
		}

		// compute new solution point -> r = 2*c - T(n)
		r.clear();
		for (int i = 0; i < n; i++)
		{
			r.push_back(2 * c.at(i) - P.at(T.at(n)).at(i));
		}

		if (is_in_domain(r))
		{
			// compute evaluation value at r
			r_val = evaluate_F(r);

			numb_of_evals++;

			if (r_val < F_max)
			{
				// update P, best and worst points
				update_P_best_worst_point();
			}
		}
		else
		{
			// compute new solution point -> r = (c + T(n))/2
			r.clear();
			for (int i = 0; i < n; i++)
			{
				r.push_back((c.at(i) + P.at(T.at(n)).at(i))/2);
			}

			if (is_in_domain(r))
			{
				// compute evaluation value at r
				r_val = evaluate_F(r);

				numb_of_evals++;

				if (r_val < F_max)
				{
					// update P, best and worst points
					update_P_best_worst_point();
				}
			}
		}

		/*wtime2 = omp_get_wtime();
		wtime = wtime2 - wtime1;
		cout << "Each loop time = " << wtime << "\n";*/

	/*} while (F_max - F_min >= eps);*/
	} while (F_max - F_min >= eps && loop_count < LOOP_MAX);
};

void crs::print_solution()
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


