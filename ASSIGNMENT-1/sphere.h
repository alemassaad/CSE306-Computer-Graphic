#ifndef SPHERE_H
#define SPHERE_H

#include "vector.h"
#include "ray.h"
#include "geometry.h"

class Sphere : public Geometry {
public:
    Vector center;
    double r;

    Sphere(const Vector& center, double r, const Vector& color, bool ism = false, bool ist = false);

    bool intersection(const Ray& ray, Vector& P, Vector& N, double& t);
    bool intersection(const Ray& ray, Vector& P, Vector& N, double& t,
                      const Vector& A, const Vector& B, const Vector& C, 
                      double& alpha, double& beta, double& gamma, int& id_mesh) override;
};

#endif 
