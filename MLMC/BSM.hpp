#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>

using boost::math::normal;

double BSM_Call(const double& S0, const double& r, const double& v, const double& T, const double& K) {

	double d1 = (1 / v * sqrt(T)) * (log(S0 / K) + (r + v * v / 2) * T);
	double d2 = d1 - v * sqrt(T);
	normal s;
	double CallBS = S0 * cdf(s,d1) - K * cdf(s, d2) * exp(-r*T);
	return CallBS;
}

double BSM_Put(const double& S0, const double& r, const double& v, const double& T, const double& K) {

	double d1 = (1 / v * sqrt(T)) * (log(S0 / K) + (r - v * v / 2) * T);
	double d2 = d1 - v * sqrt(T);
	normal s;
	double PutBS = K * cdf(s, -d2) * exp(-r * T) - S0 * cdf(s, -d1);
	return PutBS;
}

double BSM_BinCall(double& S0, const double& r, const double& v, const double& T, const double& K) {

	double d1 = (1 / v * sqrt(T)) * (log(S0 / K) + (r - v * v / 2) * T);
	double d2 = d1 - v * sqrt(T);
	normal s;
	double BinCall = exp(-r * T) * cdf(s, d1);
	return BinCall;
}

double BSM_BinPut(double& S0, const double& r, const double& v, const double& T, const double& K) {

	double d1 = (1 / v * sqrt(T)) * (log(S0 / K) + (r - v * v / 2) * T);
	double d2 = d1 - v * sqrt(T);
	normal s;
	double BinPut = exp(-r * T) * cdf(s, -d2);
	return BinPut;
}


// Binary Knock-out options

double BSM_BinDIC(double& S0, const double& r, const double& v, const double& T, const double& K, double& B) {

	double gamma = 1 - 2 * r / v * v;
	double BSM_BinDIC;
	// Regular: B =< K
	if (B <= K) {
		if (B >= S0) { BSM_BinDIC = BSM_BinCall(S0, r, v, T, K); }
		else { BSM_BinDIC = pow(S0 / B, gamma) * BSM_BinCall(B, r, v, T, K * S0 / B); }
	}
	else {
	// Reverse: B > K
		if (B >= S0) { BSM_BinDIC = BSM_BinCall(S0, r, v, T, K); }
		else { BSM_BinDIC = BSM_BinCall(S0, r, v, T, K) - exp(-r*T) + BSM_BinPut(S0, r, v, T, B) + pow(S0 / B, gamma) * BSM_BinCall(B, r, v, T, S0); }
	}
	return BSM_BinDIC;
}

double BSM_BinDIP(double& S0, const double& r, const double& v, const double& T, const double& K, double& B) {

	double gamma = 1 - 2 * r / v * v;
	double BSM_BinDIP;
	// Regular: B =< K
	if (B <= K) {
		if (B >= S0) { BSM_BinDIP = BSM_BinPut(S0, r, v, T, K); }
		else { BSM_BinDIP = pow(S0 / B, gamma) * BSM_BinPut(B, r, v, T, K * S0 / B); }
	}
	else {
		// Reverse: B > K
		if (B >= S0) { BSM_BinDIP = BSM_BinPut(S0, r, v, T, K); }
		else {
			BSM_BinDIP = BSM_BinPut(S0, r, v, T, K) - exp(-r * T) + BSM_BinCall(S0, r, v, T, B) + pow(S0 / B, gamma) * BSM_BinPut(B, r, v, T, S0);
		}
	}
	return BSM_BinDIP;
}

double BSM_BinUIC(double& S0, const double& r, const double& v, const double& T, const double& K, double& B) {

	double gamma = 1 - 2 * r / v * v;
	double BSM_BinUIC;
	// Regular: B =< K
	if (B <= K) {
		if (B <= S0) { BSM_BinUIC = BSM_BinCall(S0, r, v, T, K); }
		else { BSM_BinUIC = pow(S0 / B, gamma ) * BSM_BinCall(B, r, v, T, K * S0 / B); }
	}
	// Reverse: B > K
	else {
		if (B <= S0) { BSM_BinUIC = BSM_BinCall(S0, r, v, T, K); }
		else {
			BSM_BinUIC = BSM_BinCall(S0, r, v, T, K) - exp(-r * T) + BSM_BinPut(S0, r, v, T, B) + pow(S0 / B, gamma) * BSM_BinCall(B, r, v, T, S0);
		}
	}
	return BSM_BinUIC;
}


double BSM_BinUIP(double& S0, const double& r, const double& v, const double& T, const double& K, double& B) {

	double gamma = 1 - 2 * r / v * v;
	double BSM_BinUIP;
	// Regular: B =< K
	if (B <= K) {
		if (B <= S0) { BSM_BinUIP = BSM_BinPut(S0, r, v, T, K); }
		else { BSM_BinUIP = pow(S0 / B, gamma) * BSM_BinPut(B, r, v, T, K * S0 / B); }
	}
	// Reverse: B > K
	else {
		if (B <= S0) { BSM_BinUIP = BSM_BinPut(S0, r, v, T, K); }
		else {
			BSM_BinUIP = BSM_BinPut(S0, r, v, T, K) - exp(-r * T) + BSM_BinCall(S0, r, v, T, B) + pow(S0 / B, gamma) * BSM_BinPut(B, r, v, T, S0);
		}
	}
	return BSM_BinUIP;
}



// Knock-out options

double BSM_DIC(double& S0, const double& r, const double& v, const double& T, const double& K, double& B) {

	double gamma = 1 - 2 * r / v * v;
	double BSM_DIC;
	// Regular: B =< K
	if (B <= K) {
		if (B >= S0) { BSM_DIC = BSM_Call(S0, r, v, T, K); }
		else { BSM_DIC = pow(S0 / B, gamma - 1) * BSM_Call(B, r, v, T, K * S0 / B); }
	}
	else {
	// Reverse: B > K
		if (S0 > B) {
			BSM_DIC = BSM_Call(S0, r, v, T, K) - BSM_Call(S0, r, v, T, B) - (B - K) * BSM_BinCall(S0, r, v, T, B) + pow(S0 / B, gamma - 1) * BSM_Call(B, r, v, T, S0) + (B - K) * BSM_BinDIC(S0, r, v, T, B, B); }
		else { BSM_DIC = BSM_Call(S0, r, v, T, K); }
	}
	return BSM_DIC;
}

double BSM_DIP(double& S0, const double& r, const double& v, const double& T, const double& K, double& B) {

	double gamma = 1 - 2 * r / v * v;
	double BSM_DIP;
	// Regular: B < K
	if (B <= K) {
		if (B >= S0) { BSM_DIP = BSM_Put(S0, r, v, T, K); }
		else { BSM_DIP = pow(S0 / B, gamma - 1) * BSM_Put(B, r, v, T, S0 * K / B); } ///
	}
	else {
	// Reverse: B > K
		if (B >= S0) { BSM_DIP = BSM_Put(S0, r, v, T, K) - BSM_Put(S0, r, v, T, B) - (B - K) * BSM_BinPut(S0, r, v, T, B) + pow(S0 / B, gamma - 1) * BSM_Put(S0, r, v, T, B) + (B - K) * pow(S0 / B, gamma - 1) * BSM_BinPut(S0, r, v, T, B) /*BSM_BinDIP(S0, r, v, T, B, B)*/; }
		else { BSM_DIP = BSM_Put(S0, r, v, T, K); }
	}
	return BSM_DIP;
}

double BSM_DOP(double& S0, const double& r, const double& v, const double& T, const double& K, double& B) {
	double gamma = 1 - 2 * r / v * v;
	double BSM_DOP;
	BSM_DOP = BSM_Put(S0, r, v, T, K) - BSM_DIP(S0, r, v, T, K, B);
	return BSM_DOP;
}


double BSM_DOC(double& S0, const double& r, const double& v, const double& T, const double& K, double& B) {
	double gamma = 1 - 2 * r / v * v;
	double BSM_DOC;
	BSM_DOC = BSM_Call(S0, r, v, T, K) - BSM_DIC(S0, r, v, T, K, B);
	return BSM_DOC;
}


double BSM_UIC(double& S0, const double& r, const double& v, const double& T, const double& K, double& B) {
	double gamma = 1 - 2 * r / v * v;
	double BSM_UIC;
	if (B <= K){ // Regular: B <= K
		if (B > S0) { BSM_UIC = BSM_Call(S0, r, v, T, K); } ///
		else { BSM_UIC = pow(S0 / B, gamma - 1) * BSM_Call(S0, r, v, T, K * S0 / B); }
	}
	else { // Reverse: B > K
		if (B <= S0) { BSM_UIC = BSM_Call(S0, r, v, T, K); }
		else { 
			BSM_UIC = BSM_Call(S0, r, v, T, K) - BSM_Call(S0, r, v, T, B) - (B - K) * BSM_BinCall(S0, r, v, T, B) + pow(S0 / B, gamma - 1) * BSM_Call(S0, r, v, T, B) + (B - K) * pow(S0 / B, gamma) * BSM_BinCall(S0, r, v, T, B); }
	}
	return BSM_UIC;
}

double BSM_UOC(double& S0, const double& r, const double& v, const double& T, const double& K, double& B) {
	double gamma = 1 - 2 * r / v * v;
	double BSM_UOC;
	BSM_UOC = BSM_Call(S0, r, v, T, K) - BSM_UIC(S0, r, v, T, K, B);
	return BSM_UOC;
}


double BSM_UIP(double& S0, const double& r, const double& v, const double& T, const double& K, double& B) {

	double gamma = 1 - 2 * r / v * v;
	double BSM_UIP;
	// Regular: B < K
	if (B <= K) {
		if (B <= S0) { BSM_UIP = BSM_Put(S0, r, v, T, K); }
		else { BSM_UIP = pow(S0 / B, gamma - 1) * BSM_Put(B, r, v, T, K * S0 / B); }
	}
	else {
	// Reverse: B > K
		if (B <= S0) { BSM_UIP = BSM_Put(S0, r, v, T, K); }
		else {
			BSM_UIP = BSM_Put(S0, r, v, T, K) - BSM_Put(S0, r, v, T, B)  - (B - K) * BSM_BinPut(S0, r, v, T, B) + pow(S0 / B, gamma - 1) * BSM_Put(S0, r, v, T, B) + (B - K) * pow(S0 / B, gamma) * BSM_BinPut(S0, r, v, T, B) /*BSM_BinUIP(S0, r, v, T, B, B)*/;
		}
	}
	return BSM_UIP;
}

double BSM_UOP(double& S0, const double& r, const double& v, const double& T, const double& K, double& B) {
	double gamma = 1 - 2 * r / v * v;
	double BSM_UOP;
	BSM_UOP = BSM_Put(S0, r, v, T, K) - BSM_UIP(S0, r, v, T, K, B);
	return BSM_UOP;
}




// One and only function for barrier

//double BSM_Barrier(int& direction, int& activation, int& type, double& S0, const double& r, const double& v, const double& T, const double& K, double& B) {
//
//	double gamma = 1 - 2 * r / v * v;
//	double BSM_UOC;
//
//	if (type = 1) 
//	// Regular: B < K (payoff à la barrière est nul)
//	if (B < S0) { BSM_UOC = 0; } // BSM_Call(S0, r, v, T, K); 
//	else { BSM_UOC = BSM_Call(S0, r, v, T, K) - BSM_UIC(S0, r, v, T, K, B); }
//	return BSM_UOC;
//
//}
