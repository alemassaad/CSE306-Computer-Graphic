#include "random.h"

std::default_random_engine engine;
std::uniform_real_distribution<double> uniform(0.0, 1.0);

double generate_random() {
    return uniform(engine);
}
