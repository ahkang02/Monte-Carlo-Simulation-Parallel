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
#include "BSM.hpp"
#include "Discretization_schemes.h"
#include "Payoffs.h"

using namespace std;

int main() {

	clock_t start, end;

	double S0 = 90;
	double K = 110;
	double r = 0.05;
	double v = 0.2;
	double T = 1;
	double act = exp(-r * T);
	double B = 100;
	unsigned num_intervals = 250;

	std::cout << "\nFor S0 = " << S0 << ", K = " << K << ", r = " << r << ", v = " << v << " and B = " << B << " : ";

	vector<double> euler_spot_prices(num_intervals, S0);

	double put_euler = 0; double call_euler = 0;
	double put_euler_squared = 0; double call_euler_squared = 0;
	double MC_UIC = 0; double MC_UIP = 0; double MC_UOC = 0; double MC_UOP = 0;
	double MC_DIC = 0; double MC_DIP = 0; double MC_DOC = 0; double MC_DOP = 0;

	int num_sims = pow(10, 5);

	start = clock();

	for (int n = 0; n < num_sims; n++) {

		GBM_EULER(euler_spot_prices, r, v, T);

		put_euler += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
		call_euler += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
		put_euler_squared += (pow(Put(euler_spot_prices[num_intervals - 1], K), 2)) / (static_cast<double>(num_sims) - 1);
		call_euler_squared += (pow(Call(euler_spot_prices[num_intervals - 1], K), 2)) / (static_cast<double>(num_sims) - 1);

		if (*max_element(euler_spot_prices.begin(), euler_spot_prices.end()) >= B) {
			MC_UIC += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_UIP += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_UOC += 0; MC_UOP += 0;
		}
		else {
			MC_UIC += 0; MC_UIP += 0;
			MC_UOC += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_UOP += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
		}

		if (*min_element(euler_spot_prices.begin(), euler_spot_prices.end()) <= B) {
			MC_DIC += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_DIP += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);;
			MC_DOC += 0; MC_DOP += 0;
		}
		else {
			MC_DOC += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_DOP += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_DIC += 0; MC_DIP += 0;
		}

	}

	end = clock();

	double put_MC_std = pow((put_euler_squared - pow(put_euler, 2)), 0.5);
	double call_MC_std = pow((call_euler_squared - pow(call_euler, 2)), 0.5);

	std::cout << "\n\n\nStandard MC with nbr of simulations = " << num_sims << " :";

	std::cout << "\n\nMC Put price GBM Euler = " << put_euler << " with MC std dev = " << put_MC_std * 100 << "%";
	std::cout << " and MC std error = " << put_MC_std / (static_cast<double>(num_sims)) * 100 << "%";

	std::cout << "\nMC Call price GBM Euler = " << call_euler << " with MC std dev = " << call_MC_std * 100 << "%";
	std::cout << " and MC std error = " << call_MC_std / (static_cast<double>(num_sims)) * 100 << "%";

	std::cout << "\n\n\nMC Up-and-Out Put price Euler = " << MC_UOP;
	std::cout << "\nMC Up-and-Out Call price Euler = " << MC_UOC;

	std::cout << "\n\nMC Up-and-In Put price Euler = " << MC_UIP;
	std::cout << "\nMC Up-and-In Call price Euler = " << MC_UIC;

	std::cout << "\n\nMC Down-and-Out Put price Euler = " << MC_DOP;
	std::cout << "\nMC Down-and-Out Call price Euler = " << MC_DOC;

	std::cout << "\n\nMC Down-and-In Put price Euler = " << MC_DIP;
	std::cout << "\nMC Down-and-In Call price Euler = " << MC_DIC;

	std::cout << "\n\nMC CPU time = " << end - start << " sec";

	std::cout << "\n";

	return 0;
}
