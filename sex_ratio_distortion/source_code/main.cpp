#include "SexRatioDistortion.h"
#include "utils.h"
#include <fstream>
#include <iostream>

int main()
{
	std::ofstream df("data.txt");
	df << "\t" << "n" << "\t" << "zI" << "\t" << "zG_X_XX" << "\t" << "zG_X_XY" << "\t" << "zG_Y_XY" << "\t" << "zG_Z_ZZ" << "\t" << "zG_Z_ZW" << "\t" << "zG_W_ZW" << "\t" <<

		"wI_star" << "\t" <<
		"wG_X_XX" << "\t" << "wG_X_XY" << "\t" << "wG_Y_XY" << "\t" <<
		"wG_Z_ZZ" << "\t" << "wG_Z_ZW" << "\t" << "wG_W_ZW" << "\t" << std::endl;

	double n_max = 20.001;
	double n_t = 1.001;
	double interval = 0.001;
	while(n_t <= n_max)
	{ 
		Hamiltonian(n_t, df);

		n_t += interval;
	}

	return(0);
}