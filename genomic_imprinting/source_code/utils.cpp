#define _CRT_SECURE_NO_WARNINGS
#include "utils.h"
#include "array"


#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <string>
#include <vector>
#include <cmath>
#include "utils.h"
#include <numeric>
#include <cassert>


std::string datetime()
{
	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);

	std::ostringstream oss;
	oss << std::put_time(&tm, "%Y_%m_%d_%H_%M_%S");
	std::string str = oss.str();

	return(str);
}


double calcDerivZero(double& R, double& r, double& c, double& k, double& z_0, std::ofstream &df)
{
	// The algorithm implemented here is a basic bisection algorithm; the general principle works by
	// evaluating the equation at points x = A, x = B, and x = C = (A+B)/2 (the midpoint between A and B).
	// Check then if f(C) is sufficiently close to 0, in which case stop and return C. Otherwise, if f(A)
	// and f(C) have the same sign, repeat the procedure by updating A = C and then C = (A+B)/2. Otherwise,
	// set B = C and again update C = (A+B)/2.

	// Define range to look in
	double z_l = z_0; // minimum
	double z_h = R; // maximum
	double z_m = (z_l + z_h) / 2; // midpoint

	// Calculate derivative at minimum, maximum, and midpoint.
	double dwdz_l = Derivative(R, r, c, k, z_0, z_l);
	double dwdz_h = Derivative(R, r, c, k, z_0, z_h);
	double dwdz_m = Derivative(R, r, c, k, z_0, z_m);
	
	// Track number of iterations procedure has been running for.	
	int iteration = 0;
	bool monotonic = false;
	while (std::abs(dwdz_m) > epsilon && iteration < max_iteration)
	{
		if (0 < dwdz_l * dwdz_h)
		{
				z_m = z_h;
				monotonic = true;
				iteration = max_iteration + 1;
		} else {
			// If dwdz is positive for the low value...
			if (dwdz_l > 0)
			{ 
				// Check if it is not positive for all values within the range by checking the high
				// value. If so, the function is monotonically increasing.
			
				// ... and the midpoint value.
				if (dwdz_m > 0)
				{
					// Use current midpoint as new low
					z_l = z_m;

					// recalculate new midpoint using new low
					z_m = z_l + z_h / 2;
				}

				// ... but the midpoint value is negative
				else {

					// Use current midpoint as new high 
					z_h = z_m;

					// recalculate new midpoint
					z_m = (z_l + z_h) / 2;
				}

			}
			else {
				// Check if it is not negative for all values within the range by checking the high
				// value. If so, the function is monotonically decreasing.
				if (dwdz_h <= 0)
				{
					z_m = z_h;
					monotonic = true;
					break;
				}

				// Logic from above if-statement is inverted, but works in similar ways
				if (dwdz_m < 0)
				{
					z_l = z_m;
					z_m = z_l + z_h / 2;

				}
				else {
					z_h = z_m;
					z_m = (z_l + z_h) / 2;
				}
			}
		}

		// Recalculate derivatives
		dwdz_l = Derivative(R, r, c, k, z_0, z_l);
		dwdz_h = Derivative(R, r, c, k, z_0, z_h);
		dwdz_m = Derivative(R, r, c, k, z_0, z_m);
		++iteration;
	}

	double r_pat = 0.0;
	double r_mat = 0.5;
	double r_org = 0.25;

	double w_pat = calcEquilFitness(R, r_pat, c, k, z_0, z_m);
	double w_mat = calcEquilFitness(R, r_mat, c, k, z_0, z_m);
	double w_org = calcEquilFitness(R, r_org, c, k, z_0, z_m);
	bool i_max = iteration == max_iteration ? true : false;
	df << iteration << "\t"<< i_max << "\t" << monotonic << "\t" << R << "\t" << r << "\t" << c << "\t" << k << 
		"\t" << z_0 << "\t" << z_m << "\t" << w_pat << "\t" << w_pat << "\t" << w_pat << std::endl;
	return(z_m, dwdz_m);
}

double calcEquilFitness(double& R, double& r, double& c, double& k, double& z_0, double& z_hat)
{
	double w = k * (1 - exp(-c * (z_hat - z_0))) + r * R * k * (1 - exp(-c * (z_hat - z_0))) / z_hat;
	return(w);
}


double Derivative(double& R, double& r, double& c, double& k, double& z_0, double z)
{
	double dwdz = c * k * exp(-c * (z - z_0)) + R * k * r * ((c * z + 1) * exp(-c * (z - z_0)) - 1) / pow(z, 2);
	return(dwdz);
}