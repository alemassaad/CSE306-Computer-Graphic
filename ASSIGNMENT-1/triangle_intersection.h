#ifndef TRIANGLE_INTERSECTION_H
#define TRIANGLE_INTERSECTION_H

#include "vector.h"
#include "ray.h"

bool Triangle_Intersection(const Ray& ray, const Vector& A, const Vector& B, const Vector& C, 
                           double& alpha, double& beta, double& gamma, double& t, Vector& N);

#endif 
