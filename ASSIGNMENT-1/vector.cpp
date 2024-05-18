#include "vector.h"
#include <cmath>

Vector::Vector(double x, double y, double z) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
}

double Vector::norm2() const {
    return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
}

double Vector::norm() const {
    return std::sqrt(norm2());
}

void Vector::normalize() {
    double n = norm();
    if (n > 0) {  // Avoid division by zero
        double inv_n = 1.0 / n;
        data[0] *= inv_n;
        data[1] *= inv_n;
        data[2] *= inv_n;
    }
}

double Vector::operator[](int i) const { return data[i]; }
double& Vector::operator[](int i) { return data[i]; }

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator-(const Vector& b) {
    return Vector(-b[0], -b[1], -b[2]);
}

Vector operator*(const double a, const Vector& b) {
    return Vector(a * b[0], a * b[1], a * b[2]);
}

Vector operator*(const Vector& a, const double b) {
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}

Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector operator/(const Vector& a, const double b) {
    double inv_b = 1.0 / b;
    return Vector(a[0] * inv_b, a[1] * inv_b, a[2] * inv_b);
}

double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

double sqr(double x) {
    return x * x;
}
