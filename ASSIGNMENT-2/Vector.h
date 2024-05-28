#ifndef VECTOR_H
#define VECTOR_H

class Vector {
public:
    Vector(double x = 0, double y = 0, double z = 0);
    double lengthSquared() const;
    double length() const;
    void normalize();
    double operator[](int index) const;
    double& operator[](int index);
private:
    double components[3];
};

Vector operator+(const Vector& a, const Vector& b);
Vector operator-(const Vector& a, const Vector& b);
Vector operator-(const Vector& b);
Vector operator*(const double scalar, const Vector& vec);
Vector operator*(const Vector& vec, const double scalar);
Vector operator/(const Vector& vec, const double scalar);
double dotProduct(const Vector& a, const Vector& b);
Vector crossProduct(const Vector& a, const Vector& b);

#endif // VECTOR_H
