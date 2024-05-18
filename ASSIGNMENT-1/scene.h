#pragma once
#include <vector>
#include "ray.h"
#include "vector.h"
#include "sphere.h"
#include "geometry.h"
#include "random_cos.h"

class Scene {
public:
    Scene();
    void addObject(Geometry* s);
    bool intersection(const Ray& ray, Vector& P, Vector& N, double& t, int& ID, int& id_mesh);
    Vector getColor(const Ray& ray, int bounce);

    std::vector<Geometry*> objects;
};
