#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include "stb_image_write.h"
#include "stb_image.h"
#include <math.h>
#include <iostream>
#include <chrono>
#include <list>
#include <random>
#include "vector.h"
#include "ray.h"
#include "sphere.h"
#include "scene.h"
#include "triangle_mesh.h"
#include "random_cos.h"

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0, 1);

std::vector<int> ob_list;

void boxMuller(double stdev, double &x, double &y) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * stdev;
    y = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2) * stdev;
}

int main() {
    auto start = std::chrono::steady_clock::now();
    int W = 512;
    int H = 512;

    double fov = 60 * M_PI / 180;
    Vector Q(0, 0, 55);

    Sphere S1(Vector(-20, 0, 0), 10, Vector(0.7, 0.3, 0.1));
    Sphere S2(Vector(0, 0, 0), 10, Vector(0.7, 0.3, 0.1), true);
    Sphere S3(Vector(20, 0, 0), 10, Vector(0.7, 0.3, 0.1), false, true);

    Sphere Sfloor(Vector(0, -1000, 0), 990, Vector(0.1, 0.5, 0.1));
    Sphere Sceilling(Vector(0, 1000, 0), 940, Vector(0.2, 0.4, 0.2));
    Sphere Sleft(Vector(-1000, 0, 0), 940, Vector(0.1, 0.3, 0.8));
    Sphere Sright(Vector(1000, 0, 0), 940, Vector(0.7, 0.2, 0.5));
    Sphere Sback(Vector(0, 0, 1000), 940, Vector(0.1, 0.4, 0.2));
    Sphere Sfront(Vector(0, 0, -1000), 940, Vector(0.9, 0.3, 0.5));

    Scene scene;
    scene.addObject(&Sfloor);
    scene.addObject(&Sceilling);
    scene.addObject(&Sleft);
    scene.addObject(&Sright);
    scene.addObject(&Sback);
    scene.addObject(&Sfront);

    TriangleMesh mesh(Vector(0.3, 0.3, 0.3));
    mesh.readOBJ("cat.obj");
    mesh.factor_move(0.6, Vector(0, -10, 0));
    mesh.recursive_call(&mesh.root, 0, mesh.indices.size());

    std::cout << "Number of triangles in the mesh: " << mesh.indices.size() / 3 << std::endl;
    std::cout << "Number of vertices in the mesh: " << mesh.vertices.size() << std::endl;    scene.addObject(&mesh);

    ob_list.push_back(1);
    ob_list.push_back(1);
    ob_list.push_back(1);
    ob_list.push_back(1);
    ob_list.push_back(1);
    ob_list.push_back(1);
    ob_list.push_back(2);

    std::vector<unsigned char> image(W * H * 3, 0);
    int nb_paths = 20;
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color(0., 0., 0.);
            for (int k = 0; k < nb_paths; k++) {
                double x, y;
                boxMuller(1, x, y);
                x = i;
                y = j;
                Vector u(y - W / 2 + 0.5, H / 2 - x - 0.5, -W / (2 * tan(fov / 2)));
                u.normalize();
                Ray ray(Q, u);
                color = color + scene.getColor(ray, 5);
            }
            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0] / nb_paths, 1. / 2.2));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1] / nb_paths, 1. / 2.2));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2] / nb_paths, 1. / 2.2));
        }
    }
    stbi_write_png("image_rmi.png", W, H, 3, &image[0], 0);

    auto finish = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "Time elapsed: " << elapsed << " microseconds = " << elapsed / 1e6 << " seconds" << std::endl;

    return 0;
}
