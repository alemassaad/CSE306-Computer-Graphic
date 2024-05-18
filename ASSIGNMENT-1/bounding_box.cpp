#include "bounding_box.h"
#include <algorithm>

BoundingBox::BoundingBox(Vector min, Vector max) : min(min), max(max) {}

bool BoundingBox::intersect(const Ray &ray) {
    double t0[3], t1[3];
    for (int i = 0; i < 3; i++) {
        double x, y;
        x = std::min((min[i] - ray.O[i]) / ray.dirrection[i], (max[i] - ray.O[i]) / ray.dirrection[i]);
        y = std::max((min[i] - ray.O[i]) / ray.dirrection[i], (max[i] - ray.O[i]) / ray.dirrection[i]);
        t0[i] = x;
        t1[i] = y;
    }

    double maxt0, mint1;
    maxt0 = std::max(t0[0], std::max(t0[1], t0[2]));
    mint1 = std::min(t1[0], std::min(t1[1], t1[2]));

    if (mint1 > maxt0 && maxt0 > 0) {
        return true;
    }

    return false;
}
