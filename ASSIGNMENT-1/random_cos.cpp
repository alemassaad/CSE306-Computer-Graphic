#include "random_cos.h"
#include "random.h"

Vector random_cos(const Vector &N) {
    double r1 = generate_random();
    double r2 = generate_random();

    double x = cos(2 * M_PI * r1) * sqrt(1 - r2);
    double y = sin(2 * M_PI * r1) * sqrt(1 - r2);
    double z = sqrt(r2);

    Vector T1, T2;
    double min_num = std::min({N[0], N[1], N[2]});

    if (min_num == N[0]) {
        T1 = Vector(0, N[2], -N[1]);
    } else if (min_num == N[1]) {
        T1 = Vector(N[2], 0, -N[0]);
    } else {
        T1 = Vector(N[1], -N[0], 0);
    }

    T2 = cross(N, T1);
    T1.normalize();
    T2.normalize();
    return x * T1 + y * T2 + z * N;
}
