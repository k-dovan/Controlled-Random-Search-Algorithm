#include <vector>

class crs_omp
{
private:
	// n-no. dimensions of the function, N-no. sampled points
	int n, N;
	const double eps = 0.001;
	// min and max values of the domain  D
	const double D_max = 39.9999999999999999;
	// use it to form N from n
	const int multiplier = 50;
	int numb_of_evals = 0;
	// sample data
	std::vector<std::vector<double>> P;
	// evaluation values of sample's vectors
	std::vector<double> F;

	// temporary min, max evaluation values and it's indexes respectively
	double F_min, F_max;
	int idx_min, idx_max;

private:
	double rand_val();
	double evaluate_F(std::vector<double>);
	bool is_in_domain(std::vector<double>);

public:
	// n-no. dimensions of the func
	void init(int n);
	void crs_exec();
	void print_solution();
};

