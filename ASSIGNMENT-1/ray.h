#ifndef RAY_H
#define RAY_H

#include "vector.h"

class Ray {
public:
    Vector O, dirrection;
    Ray(const Vector& pos, const Vector& dir);
};

#endif 
