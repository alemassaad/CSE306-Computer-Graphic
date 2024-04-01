#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <limits>
#include <cmath>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};
 
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class Ray {
public:
    Vector origin;
    Vector direction;
    Ray(const Vector& o, const Vector& d) : origin(o), direction(d) {
        direction.normalize(); 
    }
};

class Sphere {
public:
    Vector center;
    double radius;
    Vector color;

    Sphere(const Vector& c, double r, const Vector& col) : center(c), radius(r), color(col) {}

    bool intersect(const Ray& ray, double& t) const {
        Vector O_minus_C = ray.origin - center;
        //since direction is a unit vector, a==1
        double b = 2 * dot(ray.direction, O_minus_C);
        double c = dot(O_minus_C, O_minus_C) - radius * radius;
        double delta = b * b - 4 * c;
        
        if (delta < 0) return false; //no intersection
        
        double sqrtDelta = sqrt(delta);
        double t0 = (-b - sqrtDelta) / 2;
        double t1 = (-b + sqrtDelta) / 2;

        if (t0 < 0) t0 = t1; //if t0 is negative, use t1 instead
        if (t0 < 0) return false; //both t0 and t1 are negative, no intersection

        t = t0; //t0 is the closest intersection
        return true;
    }
};

class Scene {
public:
    std::vector<Sphere> spheres;
    void add(const Sphere& sphere) {
        spheres.push_back(sphere);
    }
    
    bool intersect(const Ray& ray, Vector& hit, Vector& N, const Sphere*& closest_sphere) const {
        double min_t = std::numeric_limits<double>::max();
        bool found = false;
        closest_sphere = nullptr;

        for (const auto& sphere : spheres) {
            double t = 0;
            if (sphere.intersect(ray, t) && t < min_t) {
                min_t = t;
                found = true;
                closest_sphere = &sphere;
            }
        }

        if (found) {
            hit = ray.origin + ray.direction * min_t;
            N = (hit - closest_sphere->center) / closest_sphere->radius;
            N.normalize();
        }

        return found;
    }
};

Vector computeColor(const Vector& hit, const Vector& N, const Vector& materialColor, const Vector& lightPos, double intensity) {    
    Vector L = lightPos - hit;
    double distance2 = L.norm2();
    L.normalize();
    double lambertian = std::max(dot(N, L), 0.0);
    Vector color = (intensity / (4 * M_PI * distance2)) * materialColor * lambertian / M_PI;
    color[0] = std::min(std::max(color[0], 0.0), 255.0);
    color[1] = std::min(std::max(color[1], 0.0), 255.0);
    color[2] = std::min(std::max(color[2], 0.0), 255.0);
    return color;
};

int main() {
    int W = 512;
    int H = 512;
    double fov = M_PI / 3.;
    std::vector<unsigned char> image(W * H * 3, 0);
    Vector camera(0, 0, 0);
    Scene scene;
    
    // Add a red sphere to the scene
    Vector sphereColor(255, 0, 0); // Red color
    scene.add(Sphere(Vector(0, 0, 55), 10, sphereColor));
    
    double scale = tan(fov * 0.5);
    double imageAspectRatio = W / (double)H;
    Vector lightPos(50, 50, 50);
    double intensity = 1000;
    
    for (int j = 0; j < H; ++j) {
        for (int i = 0; i < W; ++i) {
            double x = (2 * (i + 0.5) / (double)W - 1) * imageAspectRatio * scale;
            double y = (1 - 2 * (j + 0.5) / (double)H) * scale;
            Vector pixel(x, y, -1);
            pixel.normalize();
            
            Ray ray(camera, pixel);
            Vector hit, N;
            const Sphere* closest_sphere = nullptr;
            if (scene.intersect(ray, hit, N, closest_sphere)) {
                Vector color = computeColor(hit, N, closest_sphere->color, lightPos, intensity);
                image[(j * W + i) * 3 + 0] = (unsigned char)(color[0]);
                image[(j * W + i) * 3 + 1] = (unsigned char)(color[1]);
                image[(j * W + i) * 3 + 2] = (unsigned char)(color[2]);
            }
        }
    }
    
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    return 0;
}
