#include "SexRatioDistortion.h"


double Fitness(double& n, double& z, double& z_competitor)
{
	double z_hat = ((n - 1) * z_competitor + z) / n;

	double w_daughters = (1 - z);
	double w_sons = (z / z_hat) * (1 - z_hat);

	return(w_daughters + w_sons);

}

double Hamiltonian(double& n, std::ofstream& d)
{

	// Optimal sex ratio (individual)
	double zI = (n - 1) / (2 * n);
	
	// Optimal sex ratio (genic elements)
	double zG_X_XX = zI * (1 - (1 / (4 * n - 1)));
	double zG_X_XY = zI * (1 - (4 * n - 3) / (4 * n - 1));
	double zG_Y_XY = zI * 2;

	double zG_Z_ZZ = zI * (1 + (1 / (4 * n - 1)));
	double zG_Z_ZW = zI * (1 + (4 * n - 3) / (4 * n - 1));
	double zG_W_ZW = 0;

	// Individual fitness given specific optimal sex ratios
	double wI = Fitness(n, zI, zI);
	
	double wG_X_XX = Fitness(n, zG_X_XX, zI);
	double wG_X_XY = Fitness(n, zG_X_XY, zI);
	double wG_Y_XY = Fitness(n, zG_Y_XY, zI);

	double wG_Z_ZZ = Fitness(n, zG_Z_ZZ, zI);
	double wG_Z_ZW = Fitness(n, zG_Z_ZW, zI);
	double wG_W_ZW = Fitness(n, zG_W_ZW, zI);


	d << "\t" << n << "\t" << zI << "\t" <<
		zG_X_XX << "\t" << zG_X_XY << "\t" << zG_Y_XY << "\t" <<
		zG_Z_ZZ << "\t" << zG_Z_ZW << "\t" << zG_W_ZW << "\t" <<

		wI << "\t" <<
		wG_X_XX << "\t" << wG_X_XY << "\t" << wG_Y_XY << "\t" <<
		wG_Z_ZZ << "\t" << wG_Z_ZW << "\t" << wG_W_ZW << "\t" << std::endl;

	return(zI);
}