#include "Vector.h"
#include <cmath>

Vector::Vector(double x, double y, double z) {
    components[0] = x;
    components[1] = y;
    components[2] = z;
}

double Vector::lengthSquared() const {
    return components[0] * components[0] + components[1] * components[1] + components[2] * components[2];
}

double Vector::length() const {
    return sqrt(lengthSquared());
}

void Vector::normalize() {
    double len = length();
    for (int i = 0; i < 3; ++i) {
        components[i] /= len;
    }
}

double Vector::operator[](int index) const {
    return components[index];
}

double& Vector::operator[](int index) {
    return components[index];
}

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector operator-(const Vector& b) {
    return Vector(-b[0], -b[1], -b[2]);
}

Vector operator*(const double scalar, const Vector& vec) {
    return Vector(scalar * vec[0], scalar * vec[1], scalar * vec[2]);
}

Vector operator*(const Vector& vec, const double scalar) {
    return Vector(vec[0] * scalar, vec[1] * scalar, vec[2] * scalar);
}

Vector operator/(const Vector& vec, const double scalar) {
    return Vector(vec[0] / scalar, vec[1] / scalar, vec[2] / scalar);
}

double dotProduct(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector crossProduct(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
