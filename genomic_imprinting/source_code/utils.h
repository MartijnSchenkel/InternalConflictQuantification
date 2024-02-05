#pragma once
#include "parameters.h"
#include "math.h"
#include <iostream>
#include <fstream>

// Function to get the date, used to write output with the starting date of the simulation to simplify data management.
std::string datetime();

// Fitness functions + associated.
double calcDerivZero(double& R, double& r, double& c, double& k, double& z_0, std::ofstream& df);
double calcEquilFitness(double& R, double& r, double& c, double& k, double& z_0, double& z_hat);
double Derivative(double& R, double& r, double& c, double& k, double& z_0, double z);