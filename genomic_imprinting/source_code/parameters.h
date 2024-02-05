#pragma once
#include <vector>

// Derivative parameters
constexpr double R_min = 4.0;
constexpr double R_max = 20.0;
constexpr double R_step_size = 0.05;


constexpr double z_0_min = 0.5;
constexpr double z_0_max = 4.0;
constexpr double z_0_step_size = 0.01;

// Approximation algorithm parameters
constexpr int max_iteration = 10000; // maximum number of calculations.
constexpr int iter_skip = 1000;
constexpr double epsilon = 0.0000000001; // mimimum level of precision before stopping (unless max. iterations reached)