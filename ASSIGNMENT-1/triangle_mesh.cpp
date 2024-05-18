#include "triangle_mesh.h"
#include "triangle_intersection.h"
#include <algorithm>
#include <list>

TriangleMesh::~TriangleMesh() {}

TriangleMesh::TriangleMesh(const Vector color) {
    this->albedo = color;
    this->isMirror = false;
    this->isTrans = false;
}

bool TriangleMesh::intersection(const Ray& ray, Vector& P, Vector& N, double& t,
                                const Vector& A, const Vector& B, const Vector& C, 
                                double& alpha, double& beta, double& gamma, int& id_mesh) {
    if (!root.box.intersect(ray)) {
        return false;
    }

    std::list<BVH*> nodes_to_visit;
    nodes_to_visit.push_front(&root);
    while (!nodes_to_visit.empty()) {
        BVH* curNode = nodes_to_visit.back();
        nodes_to_visit.pop_back();

        if (curNode->left != nullptr && curNode->right != nullptr) {
            if (curNode->left->box.intersect(ray)) {
                nodes_to_visit.push_back(curNode->left);
            }

            if (curNode->right->box.intersect(ray)) {
                nodes_to_visit.push_back(curNode->right);
            }
        }
        if ((curNode->left == nullptr) || (curNode->right == nullptr)) {
            t = 3E10;
            bool f_inter = false;
            int id = -1;
            for (int i = curNode->head; i < curNode->tail; i++) {
                Vector localN;
                double localt, localalpha, localbeta, localgamma;
                if (Triangle_Intersection(ray, vertices[indices[i].vtxi], vertices[indices[i].vtxj], vertices[indices[i].vtxk], localalpha, localbeta, localgamma, localt, localN)) {
                    if (localt < t) {
                        t = localt;
                        P = ray.O + ray.dirrection * t;
                        N = localN;
                        id = i;
                        f_inter = true;
                    }
                }
            }
            id_mesh = id;
            return f_inter;
        }
    }
    return false;
}

void TriangleMesh::Compute_min_max(Vector& min, Vector& max, int l, int r) {
    min = vertices[indices[l].vtxi];
    max = vertices[indices[l].vtxi];
    for (int i = l; i < r; i++) {
        min[0] = std::min(min[0], vertices[indices[i].vtxi][0]);
        min[1] = std::min(min[1], vertices[indices[i].vtxi][1]);
        min[2] = std::min(min[2], vertices[indices[i].vtxi][2]);
        max[0] = std::max(max[0], vertices[indices[i].vtxi][0]);
        max[1] = std::max(max[1], vertices[indices[i].vtxi][1]);
        max[2] = std::max(max[2], vertices[indices[i].vtxi][2]);

        min[0] = std::min(min[0], vertices[indices[i].vtxj][0]);
        min[1] = std::min(min[1], vertices[indices[i].vtxj][1]);
        min[2] = std::min(min[2], vertices[indices[i].vtxj][2]);
        max[0] = std::max(max[0], vertices[indices[i].vtxj][0]);
        max[1] = std::max(max[1], vertices[indices[i].vtxj][1]);
        max[2] = std::max(max[2], vertices[indices[i].vtxj][2]);

        min[0] = std::min(min[0], vertices[indices[i].vtxk][0]);
        min[1] = std::min(min[1], vertices[indices[i].vtxk][1]);
        min[2] = std::min(min[2], vertices[indices[i].vtxk][2]);
        max[0] = std::max(max[0], vertices[indices[i].vtxk][0]);
        max[1] = std::max(max[1], vertices[indices[i].vtxk][1]);
        max[2] = std::max(max[2], vertices[indices[i].vtxk][2]);
    }
}

void TriangleMesh::recursive_call(BVH* H, int l, int r) {
    Vector min, max;
    Compute_min_max(min, max, l, r);
    H->box = BoundingBox(min, max);
    H->head = l;
    H->tail = r;
    H->left = nullptr;
    H->right = nullptr;
    Vector diag = max - min;
    int diag_max = 0;
    if (diag[1] > diag[0]) {
        diag_max = 1;
    }
    if (diag[2] > diag[diag_max]) {
        diag_max = 2;
    }

    int pivot_index = l;
    Vector middle = diag * 0.5 + min;
    double middle_axis = middle[diag_max];
    for (int i = l; i < r; i++) {
        Vector barycenter = (vertices[indices[i].vtxi] + vertices[indices[i].vtxj] + vertices[indices[i].vtxk]) / 3;
        if (barycenter[diag_max] < middle_axis) {
            std::swap(indices[i], indices[pivot_index]);
            pivot_index++;
        }
    }

    if (pivot_index - l <= 5 || r - pivot_index <= 5 || r - l <= 10) {
        return;
    }
    H->left = new BVH;
    H->right = new BVH;
    recursive_call(H->left, l, pivot_index);
    recursive_call(H->right, pivot_index, r);
}

void TriangleMesh::factor_move(double s, const Vector& t) {
    for (int i = 0; i < vertices.size(); i++) {
        vertices[i] = vertices[i] * s + t;
    }
}

void TriangleMesh::readOBJ(const char* obj) {
    char matfile[255];
    char grp[255];

    FILE* f;
    f = fopen(obj, "r");
    int curGroup = -1;
    while (!feof(f)) {
        char line[255];
        if (!fgets(line, 255, f)) break;

        std::string linetrim(line);
        linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
        strcpy(line, linetrim.c_str());

        if (line[0] == 'u' && line[1] == 's') {
            sscanf(line, "usemtl %[^\n]\n", grp);
            curGroup++;
        }

        if (line[0] == 'v' && line[1] == ' ') {
            Vector vec;

            Vector col;
            if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                col[0] = std::min(1., std::max(0., col[0]));
                col[1] = std::min(1., std::max(0., col[1]));
                col[2] = std::min(1., std::max(0., col[2]));

                vertices.push_back(vec);
                vertexcolors.push_back(col);

            } else {
                sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                vertices.push_back(vec);
            }
        }
        if (line[0] == 'v' && line[1] == 'n') {
            Vector vec;
            sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
            normals.push_back(vec);
        }
        if (line[0] == 'v' && line[1] == 't') {
            Vector vec;
            sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
            uvs.push_back(vec);
        }
        if (line[0] == 'f') {
            TriangleIndices t;
            int i0, i1, i2, i3;
            int j0, j1, j2, j3;
            int k0, k1, k2, k3;
            int nn;
            t.group = curGroup;

            char* consumedline = line + 1;
            int offset;

            nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
            if (nn == 9) {
                if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                if (j0 < 0) t.uvi = uvs.size() + j0; else t.uvi = j0 - 1;
                if (j1 < 0) t.uvj = uvs.size() + j1; else t.uvj = j1 - 1;
                if (j2 < 0) t.uvk = uvs.size() + j2; else t.uvk = j2 - 1;
                if (k0 < 0) t.ni = normals.size() + k0; else t.ni = k0 - 1;
                if (k1 < 0) t.nj = normals.size() + k1; else t.nj = k1 - 1;
                if (k2 < 0) t.nk = normals.size() + k2; else t.nk = k2 - 1;
                indices.push_back(t);
            } else {
                nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                if (nn == 6) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else t.uvk = j2 - 1;
                    indices.push_back(t);
                } else {
                    nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                    if (nn == 3) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        indices.push_back(t);
                    } else {
                        nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (k0 < 0) t.ni = normals.size() + k0; else t.ni = k0 - 1;
                        if (k1 < 0) t.nj = normals.size() + k1; else t.nj = k1 - 1;
                        if (k2 < 0) t.nk = normals.size() + k2; else t.nk = k2 - 1;
                        indices.push_back(t);
                    }
                }
            }

            consumedline = consumedline + offset;

            while (true) {
                if (consumedline[0] == '\n') break;
                if (consumedline[0] == '\0') break;
                nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                TriangleIndices t2;
                t2.group = curGroup;
                if (nn == 3) {
                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else t2.vtxi = i0 - 1;
                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else t2.vtxj = i2 - 1;
                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else t2.vtxk = i3 - 1;
                    if (j0 < 0) t2.uvi = uvs.size() + j0; else t2.uvi = j0 - 1;
                    if (j2 < 0) t2.uvj = uvs.size() + j2; else t2.uvj = j2 - 1;
                    if (j3 < 0) t2.uvk = uvs.size() + j3; else t2.uvk = j3 - 1;
                    if (k0 < 0) t2.ni = normals.size() + k0; else t2.ni = k0 - 1;
                    if (k2 < 0) t2.nj = normals.size() + k2; else t2.nj = k2 - 1;
                    if (k3 < 0) t2.nk = normals.size() + k3; else t2.nk = k3 - 1;
                    indices.push_back(t2);
                    consumedline = consumedline + offset;
                    i2 = i3;
                    j2 = j3;
                    k2 = k3;
                } else {
                    nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                    if (nn == 2) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else t2.uvk = j3 - 1;
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        indices.push_back(t2);
                    } else {
                        nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else t2.vtxk = i3 - 1;
                            if (k0 < 0) t2.ni = normals.size() + k0; else t2.ni = k0 - 1;
                            if (k2 < 0) t2.nj = normals.size() + k2; else t2.nj = k2 - 1;
                            if (k3 < 0) t2.nk = normals.size() + k3; else t2.nk = k3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            k2 = k3;
                            indices.push_back(t2);
                        } else {
                            nn = sscanf(consumedline, "%u%n", &i3, &offset);
                            if (nn == 1) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else t2.vtxk = i3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                indices.push_back(t2);
                            } else {
                                consumedline = consumedline + 1;
                            }
                        }
                    }
                }
            }
        }
    }
    fclose(f);
}

