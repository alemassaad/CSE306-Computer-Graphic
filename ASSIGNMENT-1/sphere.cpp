#include "sphere.h"
#include <cmath>

Sphere::Sphere(const Vector& center, double r, const Vector& color, bool ism, bool ist)
    : center(center), r(r) {
    this->albedo = color;
    this->isMirror = ism;
    this->isTrans = ist;
}

bool Sphere::intersection(const Ray& ray, Vector& P, Vector& N, double& t) {
    Vector OC = ray.O - center;
    double b = dot(ray.dirrection, OC);
    double c = dot(OC, OC) - r * r;
    double delta = b * b - c;

    if (delta < 0) {
        return false;
    }

    double sqrtdelta = std::sqrt(delta);
    double t1 = -b - sqrtdelta;
    double t2 = -b + sqrtdelta;

    if (t2 < 0) {
        return false;
    }

    t = (t1 > 0) ? t1 : t2;
    P = ray.O + t * ray.dirrection;
    N = P - center;
    N.normalize();
    return true;
}

bool Sphere::intersection(const Ray& ray, Vector& P, Vector& N, double& t,
                          const Vector& A, const Vector& B, const Vector& C, 
                          double& alpha, double& beta, double& gamma, int& id_mesh) {
    return intersection(ray, P, N, t);
}
