#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>

using namespace std;

double gaussian_box_muller() {
    double x = 0.0; double y = 0.0; double euclid_sq = 0.0;

    do {
        x = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
        y = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
        euclid_sq = x * x + y * y;
    } while (euclid_sq >= 1.0);

    return x * sqrt(-2 * log(euclid_sq) / euclid_sq);
}


void GBM_EULER(std::vector<double>& spot_prices, const double& r, const double& v, const double& T) { 

    double dt = T / static_cast<double>(spot_prices.size());
    double drift = exp(dt * (r - 0.5 * v * v));
    double vol = sqrt(v * v * dt);

    for (int i = 1; i < spot_prices.size(); i++) {
        double gauss_bm = gaussian_box_muller();
        spot_prices[i] = spot_prices[i - 1] * drift * exp(vol * gauss_bm);
    }
}



void HESTON_EULER(std::vector<double>& spot_prices, std::vector<double>& vol_paths, const double& r, const double& v, const double& T,
    const double& kappa, const double& theta, const double& xi, const double& rho) {

    double dt = T / static_cast<double>(spot_prices.size());

    for (int i = 1; i < spot_prices.size(); i++) {

        double gauss_bm1 = gaussian_box_muller();
        double gauss_bm2 = gaussian_box_muller();
        gauss_bm2 = rho * gauss_bm1 + sqrt(1 - rho * rho) * gauss_bm2;

        double drift = exp(dt * (r - 0.5 * vol_paths[i - 1] * vol_paths[i - 1]));
        double vol = sqrt(vol_paths[i - 1] * vol_paths[i - 1] * dt);

        spot_prices[i] = spot_prices[i - 1] * drift * exp(vol * gauss_bm1);
        vol_paths[i] = vol_paths[i - 1] + kappa * (theta - vol_paths[i - 1]) * dt + xi * sqrt(max(vol_paths[i - 1],0.0) * dt) * gauss_bm2;
    }
}



void GBM_EULER_MLMC(std::vector<double>& Xf, std::vector<double>& Xc, const double& time_steps, const double& M, std::vector<double>& Gauss,
    const double& r, const double& v, const double& T) {

    double dt_Xf = T / static_cast<double>(Xf.size());
    double dt_Xc = T / static_cast<double>(Xc.size());

    for (int t = 1; t < Xf.size(); t++) {
        Gauss[t-1] = gaussian_box_muller();
        Xf[t] = Xf[t - 1] * exp((r - 0.5 * v * v) * dt_Xf + v * sqrt(dt_Xf) * Gauss[t - 1]);
    }

    for (int t = 1; t < Xc.size(); t++) {
        Xc[t] = Xc[t - 1] * exp((r - 0.5 * v * v) * dt_Xc + v * sqrt(dt_Xc) * (/*Gauss[t + M] +*/ Gauss[t - 1])); // different time steps but same brownian
    }
}



void GBM_MILSTEIN(std::vector<double>& spot_prices, const double& r, const double& v, const double& T) {

    double dt = T / static_cast<double>(spot_prices.size());

    for (int t = 1; t < spot_prices.size(); t++) {
        double gauss_bm = gaussian_box_muller();
        spot_prices[t] = spot_prices[t - 1] + r * dt + v * (sqrt(dt)*gauss_bm) + 0.5 * v * v * (pow(sqrt(dt) * gauss_bm, 2) - dt);
    }
}



void Brownian_Bridge(std::vector<double>& spot_prices, 
                    std::vector<double>& p_up, std::vector<double>& p_down, 
                    const double& r, const double& v, const double& T, double& B, 
                    double& p_up_cross, double& p_down_cross) {

    double dt = T / static_cast<double>(spot_prices.size());
    double drift = exp(dt * (r - 0.5 * v * v));
    double vol = sqrt(v * v * dt);

    for (int i = 0; i < spot_prices.size() - 1; i++) {
        double gauss_bm = gaussian_box_muller();
        spot_prices[i+1] = spot_prices[i] * drift * exp(vol * gauss_bm);
        p_up[i] = exp((-2 * (std::max(spot_prices[i] - B, 0.0)) * (std::max(spot_prices[i + 1] - B, 0.0))) / v * v * dt * spot_prices[i] * spot_prices[i]);
        p_down[i] = exp((-2 * (std::max(B - spot_prices[i], 0.0)) * (std::max(B - spot_prices[i + 1], 0.0))) / v * v * dt * spot_prices[i] * spot_prices[i]);
        p_up_cross *= p_up[i];
        p_down_cross *= p_down[i];
    }
}



void Brownian_Bridge_v2(std::vector<double>& X_f, std::vector<double>& X_c,
                        std::vector<double>& p_up, std::vector<double>& p_down, 
                        const double& time_steps, const double& r, const double& v, const double& T, double& B, 
                        double& p_up_cross, double& p_down_cross) {

    for (int i = 1; i < time_steps; i++) {
        double gauss_bm = gaussian_box_muller();
        X_f[i] = X_f[i - 1] + r * (T / time_steps) + v * (sqrt(T / time_steps) * gauss_bm);
        X_c[i] = X_c[i - 1] + r * (T / (time_steps - 1)) + v * (sqrt(T / (time_steps - 1)) * gauss_bm); // different time steps
        p_up[i] = exp((-2 * (std::max(X_f[i] - B, 0.0)) * (std::max(X_f[i + 1] - B, 0.0))) / v * v * (T / time_steps) * X_f[i] * X_f[i]);
        p_down[i] = exp((-2 * (std::max(B - X_f[i], 0.0)) * (std::max(B - X_f[i + 1], 0.0))) / v * v * (T / time_steps) * X_f[i] * X_f[i]);
        p_up_cross *= p_up[i];
        p_down_cross *= p_down[i];
    }
}