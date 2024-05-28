#ifndef POLYGON_H
#define POLYGON_H

#include <vector>
#include "Vector.h"

class Polygon {
public:
    std::vector<Vector> vertices;
};

Polygon clipPolygon(const Polygon& poly, const Vector& referencePoint, const Vector& edgeStart, const Vector& edgeEnd);
Polygon createVoronoiCell(const Vector* points, int totalPoints, int index);

#endif // POLYGON_H
