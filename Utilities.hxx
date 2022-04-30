#pragma once

#define _USE_MATH_DEFINES

#include <cmath>
#include <random>


constexpr double INF = std::numeric_limits<double>::infinity();
constexpr double PI  = M_PI;


std::mt19937_64 & GetPRNG();

double RandomReal();

double MaxwellBoltzmannSpeed(double M, double T);

std::pair<double,double> MaxwellBoltzmannVelocity(double M, double T);
