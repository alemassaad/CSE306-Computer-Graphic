#ifndef RANDOM_H
#define RANDOM_H

#include <random>

extern std::default_random_engine engine;
extern std::uniform_real_distribution<double> uniform;

double generate_random();

#endif 
