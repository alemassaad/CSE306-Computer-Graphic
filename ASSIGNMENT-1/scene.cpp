#include "scene.h"
#include <cmath>
#include <algorithm>
#include "random_cos.h"

Scene::Scene() {}

void Scene::addObject(Geometry* s) {
    objects.push_back(s);
}

bool Scene::intersection(const Ray& ray, Vector& P, Vector& N, double& t, int& ID, int& id_mesh) {
    t = 3E10;
    bool found_intersection = false;

    for (int i = 0; i < objects.size(); i++) {
        Vector localP, localN;
        double localt;

        Vector A, B, C;
        double alpha, beta, gamma;
        id_mesh = -1;
        int localid_mesh;

        bool inter = objects[i]->intersection(ray, localP, localN, localt, A, B, C, alpha, beta, gamma, localid_mesh);
        if (inter && localt > 0) {
            found_intersection = true;
            if (localt < t) {
                ID = i;
                t = localt;
                P = localP;
                N = localN;
                id_mesh = localid_mesh;
            }
        }
    }

    return found_intersection;
}

Vector Scene::getColor(const Ray& ray, int bounce) {
    if (bounce <= 0) {
        return Vector(0, 0, 0);
    }

    Vector P, N;
    int ID;
    double t;
    int id_mesh;
    if (!intersection(ray, P, N, t, ID, id_mesh)) {
        return Vector(0, 0, 0);
    }

    Geometry* obj = objects[ID];
    if (obj->isMirror) {
        Vector R = ray.dirrection - 2 * dot(ray.dirrection, N) * N;
        return getColor(Ray(P + 0.001 * N, R), bounce - 1);
    }

    if (obj->isTrans) {
        Vector R = ray.dirrection - 2 * dot(ray.dirrection, N) * N;
        double n1 = 1, n2 = 1.4;
        Vector N_trans = N;

        if (dot(ray.dirrection, N) > 0) {
            std::swap(n1, n2);
            N_trans = -N_trans;
        }

        double cosI = -dot(N_trans, ray.dirrection);
        double sinT2 = sqr(n1 / n2) * (1 - sqr(cosI));
        if (sinT2 > 1) {
            return getColor(Ray(P + 0.001 * N, R), bounce - 1);
        }

        Vector T = (n1 / n2) * (ray.dirrection + cosI * N_trans) - N_trans * std::sqrt(1 - sinT2);
        return getColor(Ray(P + 0.001 * T, T), bounce - 1);
    }

    Vector L(-10, 20, 40);
    Vector lightVec = L - P;
    double distLight2 = lightVec.norm2();
    lightVec.normalize();

    Ray lightRay(P + 0.001 * N, lightVec);
    Vector Plight, Nlight;
    int IDlight;
    double tlight;

    double shadow = 1.0;
    if (intersection(lightRay, Plight, Nlight, tlight, IDlight, id_mesh) && tlight * tlight < distLight2) {
        shadow = 0.0;
    }

    double I = 2E10;
    Vector L0 = shadow * I / (distLight2 * 4 * M_PI) * obj->albedo / M_PI * std::max(0.0, dot(lightVec, N));
    Ray ray_random(P, random_cos(N)); 
    Vector color = L0 + (obj->albedo * getColor(ray_random, bounce - 1));
    return color;
}
