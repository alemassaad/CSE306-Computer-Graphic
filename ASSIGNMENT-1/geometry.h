#pragma once
#include "vector.h"
#include "ray.h"

class Geometry {
public:
    Geometry() {};
    virtual bool intersection(const Ray& ray, Vector &P, Vector &N, double &t,
                              const Vector& A, const Vector& B, const Vector& C, 
                              double &alpha, double &beta, double &gamma, int &id_mesh) = 0;
    Vector albedo;
    bool isMirror, isTrans;
};
