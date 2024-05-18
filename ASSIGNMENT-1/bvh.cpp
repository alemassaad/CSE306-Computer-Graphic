#include "bvh.h"

BVH::BVH() : head(-1), tail(-1), left(nullptr), right(nullptr) {}

BVH::~BVH() {
    delete left;
    delete right;
}
