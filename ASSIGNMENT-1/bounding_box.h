#pragma once
#include "vector.h"
#include "ray.h"

class BoundingBox {
public:
    Vector min, max;
    BoundingBox(Vector min = Vector(0, 0, 0), Vector max = Vector(0, 0, 0));
    bool intersect(const Ray &ray);
};
