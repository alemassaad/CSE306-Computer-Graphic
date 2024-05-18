#pragma once
#include "vector.h"
#include "ray.h"
#include "geometry.h"
#include "triangle_indices.h"
#include "bounding_box.h"
#include "bvh.h"
#include <vector>
#include <list>

class TriangleMesh : public Geometry {
public:
    ~TriangleMesh();
    TriangleMesh(const Vector color);

    bool intersection(const Ray& ray, Vector &P, Vector &N, double &t,
                      const Vector& A, const Vector& B, const Vector& C, 
                      double &alpha, double &beta, double &gamma, int &id_mesh);

    void Compute_min_max(Vector &min, Vector &max, int l, int r);
    void recursive_call(BVH *H, int l, int r);
    void factor_move(double s, const Vector& t);
    void readOBJ(const char* obj);

    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs;
    std::vector<Vector> vertexcolors;
    BVH root;
};
