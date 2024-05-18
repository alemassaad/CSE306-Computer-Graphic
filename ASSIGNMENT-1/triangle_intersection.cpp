#include "triangle_intersection.h"
#include <cmath>

bool Triangle_Intersection(const Ray& ray, const Vector& A, const Vector& B, const Vector& C, 
                           double& alpha, double& beta, double& gamma, double& t, Vector& N) {
    Vector e1 = B - A;
    Vector e2 = C - A;
    N = cross(e1, e2);
    double b1 = dot(e2, cross(A - ray.O, ray.dirrection));
    double y1 = dot(e1, cross(A - ray.O, ray.dirrection));
    double dv = dot(ray.dirrection, N);

    beta = b1 / dv;
    gamma = -y1 / dv;
    alpha = 1 - beta - gamma;
    t = dot(A - ray.O, N) / dot(ray.dirrection, N);

    if (alpha < 0 || alpha > 1 || beta < 0 || beta > 1 || gamma < 0 || gamma > 1 || t < 0) {
        return false;
    }
    return true;
}
