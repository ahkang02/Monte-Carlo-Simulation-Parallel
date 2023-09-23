#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>

using namespace std;

double Call(const double& x, const double& K) {
	return std::max(x - K, 0.0);
}

double Put(const double& x, const double& K) {
	return std::max(K - x, 0.0);
}
