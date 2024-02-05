#include "utils.h"
#include <vector>

int main()
{
	// Define parameter value ranges
	std::vector<double> R_values;
	for (double dR = R_min; dR <= R_max;)
	{
		R_values.push_back(dR);
		dR = dR + R_step_size;
	}
	std::vector<double> z_0_values;

	for (double dz_0 = z_0_min; dz_0 <= z_0_max;)
	{
		z_0_values.push_back(dz_0);
		dz_0 = dz_0 + z_0_step_size;
	}

	std::vector<double> r_values = { 0, 0.25, 0.5 };
	std::vector<double> c_values = { 0.25, 0.5, 1.0, 1.5, 2.0};
	std::vector<double> k_values = { 1.0 };

	// Set up data file
	std::ofstream df("data.txt");

	// Add header with variables names to datafile.
	df << "iteration" << "\t" << "i_max" << "\t" << "monotonic" << "\t" << "R" << "\t" << "r" << "\t" << "c" << "\t" << "k" << 
		"\t" << "z_0" << "\t" << "z_m" << "\t" << "w_pat" << "\t" << "w_mat" << "\t" << "w_org" << std::endl;

	// Find roots
	for (size_t iR = 0, len = R_values.size(); iR < len; ++iR)
	{		
		double R = R_values[iR];
		std::cout << "Progress = "  << R << "/" << R_max << std::endl;
		for (size_t ir = 0, len = r_values.size(); ir < len; ++ir)
		{
			double r = r_values[ir];
			for (size_t ic = 0, len = c_values.size(); ic < len; ++ic)
			{
				double c = c_values[ic];
				for (size_t ik = 0, len = k_values.size(); ik < len; ++ik)
				{
					double k = k_values[ik];
					for (size_t iz_0 = 0, len = z_0_values.size(); iz_0 < len; ++iz_0)
					{
						double z_0= z_0_values[iz_0];
						
						calcDerivZero(R, r, c, k, z_0, df);
					}
				}
			}
		}
	}

	return(0);
}