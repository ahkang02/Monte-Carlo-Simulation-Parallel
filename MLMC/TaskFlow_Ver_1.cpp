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
#include "taskflow/taskflow.hpp"
#include "taskflow/algorithm/reduce.hpp"
#include "taskflow/algorithm/for_each.hpp"

using namespace std;

int main() {

	tf::Taskflow taskflow;
	tf::Executor executor(2);

	clock_t start, end;

	double S0 = 90;
	double K = 110;
	double r = 0.05;
	double v = 0.2;
	double T = 1;
	double act = exp(-r * T);
	double B = 100;
	unsigned num_intervals = 250;
	int num_sims = pow(10, 5);

	//in order to use taskflow.reduce, need declare
	//create a vector of size num_sims with all elements initialized to 0
	std::vector<double> put_euler(num_sims, 0.0);
	std::vector<double> call_euler(num_sims, 0.0);
	std::vector<double> put_euler_squared(num_sims, 0.0);
	std::vector<double> call_euler_squared(num_sims, 0.0);
	std::vector<double> MC_UIC(num_sims, 0.0);
	std::vector<double> MC_UIP(num_sims, 0.0);
	std::vector<double> MC_UOC(num_sims, 0.0);
	std::vector<double> MC_UOP(num_sims, 0.0);;
	std::vector<double> MC_DIC(num_sims, 0.0);
	std::vector<double> MC_DIP(num_sims, 0.0);
	std::vector<double> MC_DOC(num_sims, 0.0);
	std::vector<double> MC_DOP(num_sims, 0.0);

	start = clock();

	tf::Task pf = taskflow.for_each_index(0, num_sims, 1, [&](int n) {
		std::vector<double> euler_spot_prices(num_intervals, S0);
		GBM_EULER(euler_spot_prices, r, v, T);

		put_euler[n] = act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
		call_euler[n] = act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
		put_euler_squared[n] = (pow(Put(euler_spot_prices[num_intervals - 1], K), 2)) / (static_cast<double>(num_sims) - 1);
		call_euler_squared[n] = (pow(Call(euler_spot_prices[num_intervals - 1], K), 2)) / (static_cast<double>(num_sims) - 1);

		if (*max_element(euler_spot_prices.begin(), euler_spot_prices.end()) >= B) {
			MC_UIC[n] = act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_UIP[n] = act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_UOC[n] = 0;
			MC_UOP[n] = 0;
		}
		else {
			MC_UIC[n] = 0;
			MC_UIP[n] = 0;
			MC_UOC[n] = act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_UOP[n] = act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
		}

		if (*min_element(euler_spot_prices.begin(), euler_spot_prices.end()) <= B) {
			MC_DIC[n] = act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_DIP[n] = act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);;
			MC_DOC[n] = 0;
			MC_DOP[n] = 0;
		}
		else {
			MC_DOC[n] = act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_DOP[n] = act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
			MC_DIC[n] = 0;
			MC_DIP[n] = 0;
		}

		});

	//run the task and wait for them to run finish
	//use wait_for_all to synchronize the taskflow
	executor.run(taskflow);
	executor.wait_for_all();

	//sum up all element in the vector (from vector 0 to vector num_sims)
	double final_put_euler_sum = std::accumulate(put_euler.begin(), put_euler.end(), 0.0);
	double final_put_euler_squared_sum = std::accumulate(put_euler_squared.begin(), put_euler_squared.end(), 0.0);
	double final_call_euler_sum = std::accumulate(call_euler.begin(), call_euler.end(), 0.0);
	double final_call_euler_squared_sum = std::accumulate(call_euler_squared.begin(), call_euler_squared.end(), 0.0);
	double final_MC_UIC_sum = std::accumulate(MC_UIC.begin(), MC_UIC.end(), 0.0);
	double final_MC_UIP_sum = std::accumulate(MC_UIP.begin(), MC_UIP.end(), 0.0);
	double final_MC_UOC_sum = std::accumulate(MC_UOC.begin(), MC_UOC.end(), 0.0);
	double final_MC_UOP_sum = std::accumulate(MC_UOP.begin(), MC_UOP.end(), 0.0);
	double final_MC_DIC_sum = std::accumulate(MC_DIC.begin(), MC_DIC.end(), 0.0);
	double final_MC_DIP_sum = std::accumulate(MC_DIP.begin(), MC_DIP.end(), 0.0);
	double final_MC_DOC_sum = std::accumulate(MC_DOC.begin(), MC_DOC.end(), 0.0);
	double final_MC_DOP_sum = std::accumulate(MC_DOP.begin(), MC_DOP.end(), 0.0);

	end = clock();

	double put_MC_std = pow((final_put_euler_squared_sum - pow(final_put_euler_sum, 2)), 0.5);
	double call_MC_std = pow((final_call_euler_squared_sum - pow(final_call_euler_sum, 2)), 0.5);

	std::cout << "\n\n\nStandard MC with nbr of simulations = " << num_sims << " :";

	std::cout << "\n\nMC Put price GBM Euler = " << final_put_euler_sum << " with MC std dev = " << put_MC_std * 100 << "%";
	std::cout << " and MC std error = " << put_MC_std / (static_cast<double>(num_sims)) * 100 << "%";

	std::cout << "\nMC Call price GBM Euler = " << final_call_euler_sum << " with MC std dev = " << call_MC_std * 100 << "%";
	std::cout << " and MC std error = " << call_MC_std / (static_cast<double>(num_sims)) * 100 << "%";

	std::cout << "\n\n\nMC Up-and-Out Put price Euler = " << final_MC_UOP_sum;
	std::cout << "\nMC Up-and-Out Call price Euler = " << final_MC_UOC_sum;

	std::cout << "\n\nMC Up-and-In Put price Euler = " << final_MC_UIP_sum;
	std::cout << "\nMC Up-and-In Call price Euler = " << final_MC_UIC_sum;

	std::cout << "\n\nMC Down-and-Out Put price Euler = " << final_MC_DOP_sum;
	std::cout << "\nMC Down-and-Out Call price Euler = " << final_MC_DOC_sum;

	std::cout << "\n\nMC Down-and-In Put price Euler = " << final_MC_DIP_sum;
	std::cout << "\nMC Down-and-In Call price Euler = " << final_MC_DIC_sum;

	std::cout << "\n\nMC CPU time = " << end - start << " sec";

	std::cout << "\n";

	return 0;
}
