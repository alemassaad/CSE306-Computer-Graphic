#define _CRT_SECURE_NO_WARNINGS 1
#include "scene.h"
#include "stb_image_write.h"
#include <vector>
#include <iostream>
#include <cmath>

int main() {
    int W = 512;
    int H = 512;
    double fov = 60 * M_PI / 180;
    Vector Q(0, 0, 55);

    Sphere S1(Vector(-20, 0, 0), 10, Vector(1.0, 0.0, 0.0));  // Red
    Sphere S2(Vector(0, 0, 0), 10, Vector(0.0, 1.0, 0.0), true);  // Green, mirror
    Sphere S3(Vector(20, 0, 0), 10, Vector(0.0, 0.0, 1.0), false, true);  // Blue, transparent

    Sphere Sfloor(Vector(0, -1000, 0), 990, Vector(0.9, 0.5, 0.1));  // Orange
    Sphere Sceilling(Vector(0, 1000, 0), 940, Vector(0.5, 0.1, 0.9));  // Purple
    Sphere Sleft(Vector(-1000, 0, 0), 940, Vector(0.1, 0.9, 0.5));  // Teal
    Sphere Sright(Vector(1000, 0, 0), 940, Vector(0.9, 0.1, 0.3));  // Pink
    Sphere Sback(Vector(0, 0, 1000), 940, Vector(0.1, 0.3, 0.9));  // Light Blue
    Sphere Sfront(Vector(0, 0, -1000), 940, Vector(0.3, 0.9, 0.1));  // Lime Green

    Scene scene;
    scene.addObject(&S1);
    scene.addObject(&S2);
    scene.addObject(&S3);
    scene.addObject(&Sfloor);
    scene.addObject(&Sceilling);
    scene.addObject(&Sleft);
    scene.addObject(&Sright);
    scene.addObject(&Sback);
    scene.addObject(&Sfront);

    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector u(j - W / 2 + 0.5, H / 2 - i - 0.5, -W / (2 * tan(fov / 2)));
            u.normalize();
            Ray ray(Q, u);
            Vector color = scene.getColor(ray, 5);

            image[(i * W + j) * 3 + 0] = std::min(255., std::pow(color[0], 0.45));
            image[(i * W + j) * 3 + 1] = std::min(255., std::pow(color[1], 0.45));
            image[(i * W + j) * 3 + 2] = std::min(255., std::pow(color[2], 0.45));
        }
    }
    stbi_write_png("image_surfaces.png", W, H, 3, &image[0], 0);

    return 0;
}
