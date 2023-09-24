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

	std::cout << "\n\nBSM Put price = " << BSM_Put(S0, r, v, T, K);
	std::cout << "\nBSM Call price = " << BSM_Call(S0, r, v, T, K);

	std::cout << "\n\nBSM Binary Put price = " << BSM_BinPut(S0, r, v, T, K);
	std::cout << "\nBSM Binary Call price = " << BSM_BinCall(S0, r, v, T, K);

	std::cout << "\n\nBSM Up-and-Out Put price = " << BSM_UOP(S0, r, v, T, K, B);
	std::cout << "\nBSM Up-and-Out Call price = " << BSM_UOC(S0, r, v, T, K, B);

	std::cout << "\n\nBSM Up-and-In Put price = " << BSM_UIP(S0, r, v, T, K, B);
	std::cout << "\nBSM Up-and-In Call price = " << BSM_UIC(S0, r, v, T, K, B);

	std::cout << "\n\nBSM Down-and-Out Put price = " << BSM_DOP(S0, r, v, T, K, B);
	std::cout << "\nBSM Down-and-Out Call price = " << BSM_DOC(S0, r, v, T, K, B);

	std::cout << "\n\nBSM Down-and-In Put price = " << BSM_DIP(S0, r, v, T, K, B);
	std::cout << "\nBSM Down-and-In Call price = " << BSM_DIC(S0, r, v, T, K, B);

	vector<double> euler_spot_prices(num_intervals, S0);
	vector<double> milstein_spot_prices(num_intervals, S0);

	double put_euler = 0; double call_euler = 0;
	double put_milstein = 0; double call_milstein = 0;
	double put_euler_squared = 0; double call_euler_squared = 0;
	double proba_crossing = 0;
	double BB_MC_UOP = 0; double BB_MC_UIP = 0;	double BB_MC_DOP = 0; double BB_MC_DIP = 0;
	double BB_MC_UOC = 0; double BB_MC_UIC = 0;	double BB_MC_DOC = 0; double BB_MC_DIC = 0;
	double MC_UIC = 0; double MC_UIP = 0; double MC_UOC = 0; double MC_UOP = 0;
	double MC_DIC = 0; double MC_DIP = 0; double MC_DOC = 0; double MC_DOP = 0;

	double put_heston = 0; double call_heston = 0; ///
	double V0 = 0.225;
	vector<double> vol_paths(num_intervals, V0);
	double kappa = 1;
	double theta = 0.25;
	double xi = 0.02;
	double rho = 0.5;

	int num_sims = pow(10, 5);

	start = clock();

	for (int n = 0; n < num_sims; n++) {

		GBM_EULER(euler_spot_prices, r, v, T);
		//GBM_MILSTEIN(milstein_spot_prices, r, v, T);

		put_euler += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
		call_euler += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);
		put_euler_squared += (pow(Put(euler_spot_prices[num_intervals - 1], K), 2)) / (static_cast<double>(num_sims) - 1);
		call_euler_squared += (pow(Call(euler_spot_prices[num_intervals - 1], K), 2)) / (static_cast<double>(num_sims) - 1);

		/*		HESTON_EULER(euler_spot_prices, vol_paths, r, v, T, kappa, theta, xi, rho);
				put_heston += act * Put(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);//
				call_heston += act * Call(euler_spot_prices[num_intervals - 1], K) / static_cast<double>(num_sims);//
		*/
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

		// Brownian Bridge method

		vector<double> proba_up(num_intervals - 1, 1);
		vector<double> proba_down(num_intervals - 1, 1);
		double proba_up_cross = 1;
		double proba_down_cross = 1;

		Brownian_Bridge(euler_spot_prices, proba_up, proba_down, r, v, T, B, proba_up_cross, proba_down_cross);

		// proba_crossing = (std::accumulate(proba.begin(), proba.end(), proba[0], std::multiplies<double>()));

		BB_MC_UOP += act * Put(euler_spot_prices[num_intervals - 1], K) * (proba_up_cross) / static_cast<double>(num_sims);
		BB_MC_UOC += act * Call(euler_spot_prices[num_intervals - 1], K) * (proba_up_cross) / static_cast<double>(num_sims);

		BB_MC_UIP += act * Put(euler_spot_prices[num_intervals - 1], K) * (1 - proba_up_cross) / static_cast<double>(num_sims);
		BB_MC_UIC += act * Call(euler_spot_prices[num_intervals - 1], K) * (1 - proba_up_cross) / static_cast<double>(num_sims);

		BB_MC_DOP += act * Put(euler_spot_prices[num_intervals - 1], K) * (proba_down_cross) / static_cast<double>(num_sims);
		BB_MC_DOC += act * Call(euler_spot_prices[num_intervals - 1], K) * (proba_down_cross) / static_cast<double>(num_sims);

		BB_MC_DIP += act * Put(euler_spot_prices[num_intervals - 1], K) * (1 - proba_down_cross) / static_cast<double>(num_sims);
		BB_MC_DIC += act * Call(euler_spot_prices[num_intervals - 1], K) * (1 - proba_down_cross) / static_cast<double>(num_sims);

	}

	end = clock();

	double put_MC_std = pow((put_euler_squared - pow(put_euler, 2)), 0.5);
	double call_MC_std = pow((call_euler_squared - pow(call_euler, 2)), 0.5);

	std::cout << "\n\n\nStandard MC with nbr of simulations = " << num_sims << " :";

	std::cout << "\n\nMC Put price GBM Euler = " << put_euler << " with MC std dev = " << put_MC_std * 100 << "%";
	std::cout << " and MC std error = " << put_MC_std / (static_cast<double>(num_sims)) * 100 << "%";

	std::cout << "\nMC Call price GBM Euler = " << call_euler << " with MC std dev = " << call_MC_std * 100 << "%";
	std::cout << " and MC std error = " << call_MC_std / (static_cast<double>(num_sims)) * 100 << "%";

	/*	cout << "\n\nMC Put price Heston Euler = " << put_heston;
		cout << "\nMC Call price Heston Euler = " << call_heston;
	*/
	std::cout << "\n\n\nMC Up-and-Out Put price Euler = " << MC_UOP;
	std::cout << "\nMC Up-and-Out Call price Euler = " << MC_UOC;

	std::cout << "\n\nMC Up-and-In Put price Euler = " << MC_UIP;
	std::cout << "\nMC Up-and-In Call price Euler = " << MC_UIC;

	std::cout << "\n\nMC Down-and-Out Put price Euler = " << MC_DOP;
	std::cout << "\nMC Down-and-Out Call price Euler = " << MC_DOC;

	std::cout << "\n\nMC Down-and-In Put price Euler = " << MC_DIP;
	std::cout << "\nMC Down-and-In Call price Euler = " << MC_DIC;

	// Brownian Bridge

	std::cout << "\n\n\nBrownian Bridge MC :";
	std::cout << "\n\nMC BB Up-and-Out Put price Euler = " << BB_MC_UOP;
	std::cout << "\nMC BB Up-and-Out Call price Euler = " << BB_MC_UOC;

	std::cout << "\n\nMC BB Up-and-In Put price Euler = " << BB_MC_UIP;
	std::cout << "\nMC BB Up-and-In Call price Euler = " << BB_MC_UIC;

	std::cout << "\n\nMC BB Down-and-Out Put price Euler = " << BB_MC_DOP;
	std::cout << "\nMC BB Down-and-Out Call price Euler = " << BB_MC_DOC;

	std::cout << "\n\nMC BB Down-and-In Put price Euler = " << BB_MC_DIP;
	std::cout << "\nMC BB Down-and-In Call price Euler = " << BB_MC_DIC;

	std::cout << "\n\nMC CPU time = " << end - start << " sec";

	//	std::cout << "\nMC Put price with Milstein = " << put_milstein;
	//	std::cout << "\nPrice Call option with Milstein = " << call_milstein;
	//	std::cout << "\n";




	// MLMC

	int L = 4; // nbr of levels
	int M = 5; //Root of the refinements

	double mu = 0;
	double mu0 = 0;
	double sigma = 0;
	double sigma0 = 0;

	double N0 = pow(10, 4);
	double N = pow(10, 4);

	vector<double> X(N0, S0);
	double BB_MLMC_UOC_0 = 0; double BB_MLMC_UOP_0 = 0;
	double BB_MLMC_UOC_1 = 0; double BB_MLMC_UOP_1 = 0;

	start = clock();

	for (int l = 1; l < L; l++) {

		double nl = std::pow(M, l) / M;
		double Nl = int(N * (pow(pow(M, l - 1), -0.5) + pow(pow(M, l), -0.5)) / sqrt(pow(M, l - 1) + pow(M, l)));

		if (l == 1) {
			for (int n = 0; n < N0; n++) {
				GBM_EULER(X, r, v, T);

				/*vector<double> proba_up(N0 - 1, 1);
				vector<double> proba_down(N0 - 1, 1);
				double proba_up_cross = 1;
				double proba_down_cross = 1;*/

				//Brownian_Bridge(X, proba_up, proba_down, r, v, T, B, proba_up_cross, proba_down_cross);

				mu0 += act * Put(X.back(), K) / N0;
				sigma0 += (pow(Put(X.back(), K), 2)) / N0;

				//BB_MLMC_UOP_0 += act * Put(X[N0 - 1], K) * (proba_up_cross) / N0;
				//BB_MLMC_UOC_0 += act * Call(X[N0 - 1], K) * (proba_up_cross) / N0;
			}
		}
		else {
			for (int n = 0; n < static_cast<int>(Nl); n++) {

				vector<double> Xf(Nl, S0);
				vector<double> Xc(Nl / M, S0);
				vector<double> gauss(Nl);

				GBM_EULER_MLMC(Xf, Xc, nl, M, gauss, r, v, T); ///

				mu += act * (Put(Xf.back(), K) - Put(Xc.back(), K)) / Nl;
				sigma += (pow(Put(Xf.back(), K) - Put(Xc.back(), K), 2)) / (Nl - 1);

				//	vector<double> proba_up(N_l[l] - 1, 1); vector<double> proba_down(N_l[l] - 1, 1);
				//	double proba_up_cross = 1; double proba_down_cross = 1;
				//	Brownian_Bridge_v2(X_f, X_c, proba_up, proba_down, N_l[l], r, v, T, B, proba_up_cross, proba_down_cross);
				//	BB_MLMC_UOC_1 = exp(-r * T / N_l[l]) * (Put(X_f[N_l[l] - 1], K) - Put(X_c[N_l[l - 1] - 1], K)) * proba_up_cross / N_l[l];
			}
		}
	}

	end = clock();

	double put_MLMC_std = pow(sigma - pow(mu, 2) + sigma0 - pow(mu0, 2), 0.5);

	std::cout << "\n\n\nMLMC with nbr of levels = " << L << ":";
	std::cout << "\n\nMLMC Put price Euler = " << mu + mu0 << " with MLMC std dev = " << put_MLMC_std * 100 << "%";
	std::cout << " and MLMC std error = " << put_MLMC_std / (std::pow(2.0, L)) * 100 << "%";

	std::cout << "\n\nmu0 = " << mu0;
	std::cout << "\nmu = " << mu;
	std::cout << "\n\nsigma0 = " << sigma0;
	std::cout << "\nsigma = " << sigma;

	std::cout << "\n\nMLMC BB U-O-P 0 = " << BB_MLMC_UOP_0;
	std::cout << "\n\nMLMC BB U-O-P 1 = " << BB_MLMC_UOP_1;
	std::cout << "\nMLMC BB U-O-C = " << BB_MLMC_UOC_0;

	std::cout << "\n\nMLMC CPU time = " << end - start << " sec";


	std::cout << "\n";

	return 0;
}
