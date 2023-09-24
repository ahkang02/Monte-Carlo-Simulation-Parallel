#include<stdlib.h>
#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<climits>
#include<list>
#include<algorithm>
#include<numeric>
#include<vector>
#include<iostream>
#include<random>
#include<chrono>
#include<ppl.h>
#include "BSM.hpp"
#include "Discretization_schemes.h"
#include "Payoffs.h"

using namespace std;
using namespace concurrency;
using std::cout;
using std::vector;
using namespace concurrency;

int main() {
	Scheduler::SetDefaultSchedulerPolicy(SchedulerPolicy(2, MinConcurrency, 1, MaxConcurrency, 1));
	clock_t start, end;

	double S0 = 90;
	double K = 110;
	double r = 0.05;
	double v = 0.2;
	double T = 1;
	double act = exp(-r * T);
	double B = 100;
	unsigned num_intervals = 250;

	cout << "\nFor S0 = " << S0 << ", K = " << K << ", r = " << r << ", v = " << v << " and B = " << B << " : ";

	int num_sims = pow(10, 5);

	vector<double> put_euler(num_sims, 0.0);
	vector<double> call_euler(num_sims, 0.0);
	vector<double> put_euler_squared(num_sims, 0.0);
	vector<double> call_euler_squared(num_sims, 0.0);
	vector<double> MC_UIC(num_sims, 0.0);
	vector<double> MC_UIP(num_sims, 0.0);
	vector<double> MC_UOC(num_sims, 0.0);
	vector<double> MC_UOP(num_sims, 0.0);;
	vector<double> MC_DIC(num_sims, 0.0);
	vector<double> MC_DIP(num_sims, 0.0);
	vector<double> MC_DOC(num_sims, 0.0);
	vector<double> MC_DOP(num_sims, 0.0);
	double final_put_euler_sum = 0.0;
	double final_put_euler_squared_sum = 0.0;
	double final_call_euler_sum = 0.0;
	double final_call_euler_squared_sum = 0.0;
	double final_MC_UIC_sum = 0.0;
	double final_MC_UIP_sum = 0.0;
	double final_MC_UOC_sum = 0.0;
	double final_MC_UOP_sum = 0.0;
	double final_MC_DIC_sum = 0.0;
	double final_MC_DIP_sum = 0.0;
	double final_MC_DOC_sum = 0.0;
	double final_MC_DOP_sum = 0.0;

	start = clock();

	parallel_for(0, num_sims, [&](int n) {
		std::vector<double> euler_spot_prices(num_intervals, S0);
		GBM_EULER(euler_spot_prices, r, v, T);

		put_euler[n] += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
		call_euler[n] += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
		put_euler_squared[n] += (pow(Put(euler_spot_prices[num_intervals - 1], K), 2)) / (static_cast<double>(num_sims) - 1);
		call_euler_squared[n] += (pow(Call(euler_spot_prices[num_intervals - 1], K), 2)) / (static_cast<double>(num_sims) - 1);

		if (*max_element(euler_spot_prices.begin(), euler_spot_prices.end()) >= B) {
			MC_UIC[n] += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_UIP[n] += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_UOC[n] += 0; MC_UOP[n] += 0;
		}
		else {
			MC_UIC[n] += 0; MC_UIP[n] += 0;
			MC_UOC[n] += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_UOP[n] += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
		}

		if (*min_element(euler_spot_prices.begin(), euler_spot_prices.end()) <= B) {
			MC_DIC[n] += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_DIP[n] += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);;
			MC_DOC[n] += 0; MC_DOP[n] += 0;
		}
		else {
			MC_DOC[n] += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_DOP[n] += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_DIC[n] += 0; MC_DIP[n] += 0;
		}


		});

	for (int n = 0; n < num_sims; n++) {
		final_put_euler_sum += put_euler[n];
		final_put_euler_squared_sum += put_euler_squared[n];
		final_call_euler_sum += call_euler[n];
		final_call_euler_squared_sum += call_euler_squared[n];
		final_MC_UIC_sum += MC_UIC[n];
		final_MC_UIP_sum += MC_UIP[n];
		final_MC_UOC_sum += MC_UOC[n];
		final_MC_UOP_sum += MC_UOP[n];
		final_MC_DIC_sum += MC_DIC[n];
		final_MC_DIP_sum += MC_DIP[n];
		final_MC_DOC_sum += MC_DOC[n];
		final_MC_DOP_sum += MC_DOP[n];
	};


	end = clock();

	double put_MC_std = pow((final_put_euler_squared_sum - pow(final_put_euler_sum, 2)), 0.5);
	double call_MC_std = pow((final_call_euler_squared_sum - pow(final_call_euler_sum, 2)), 0.5);

	cout << "\n\n\nStandard MC with nbr of simulations = " << num_sims << " :";

	cout << "\n\nMC Put price GBM Euler = " << final_put_euler_sum << " with MC std dev = " << put_MC_std * 100 << "%";
	cout << " and MC std error = " << put_MC_std / (static_cast<double>(num_sims)) * 100 << "%";

	cout << "\nMC Call price GBM Euler = " << final_call_euler_sum << " with MC std dev = " << call_MC_std * 100 << "%";
	cout << " and MC std error = " << call_MC_std / (static_cast<double>(num_sims)) * 100 << "%";

	cout << "\n\n\nMC Up-and-Out Put price Euler = " << final_MC_UOP_sum;
	cout << "\nMC Up-and-Out Call price Euler = " << final_MC_UOC_sum;

	cout << "\n\nMC Up-and-In Put price Euler = " << final_MC_UIP_sum;
	cout << "\nMC Up-and-In Call price Euler = " << final_MC_UIC_sum;

	cout << "\n\nMC Down-and-Out Put price Euler = " << final_MC_DOP_sum;
	cout << "\nMC Down-and-Out Call price Euler = " << final_MC_DOC_sum;

	cout << "\n\nMC Down-and-In Put price Euler = " << final_MC_DIP_sum;
	cout << "\nMC Down-and-In Call price Euler = " << final_MC_DIC_sum;

	cout << "\n\nMC CPU time = " << end - start << " sec";


	cout << "\n";

	return 0;
}
