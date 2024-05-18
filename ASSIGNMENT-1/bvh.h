#pragma once
#include "bounding_box.h"

class BVH {
public:
    BVH();
    ~BVH();
    int head, tail;
    BoundingBox box;
    BVH *left, *right;
};
