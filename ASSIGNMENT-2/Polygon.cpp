#include "Polygon.h"
#include <algorithm>

namespace {
    void findIntersection(const Vector& A, const Vector& B, const Vector& lineStart, const Vector& lineEnd, Vector& intersection, bool& flag) {
        Vector normal = Vector(lineEnd[1] - lineStart[1], -lineEnd[0] + lineStart[0]);
        Vector u = lineStart;
        double t = dotProduct(u - A, normal) / dotProduct(B - A, normal);
        intersection = A + t * (B - A);
        flag = true;
        if (((dotProduct(u - B, normal) <= 0) && (dotProduct(u - A, normal) > 0)) || ((dotProduct(u - B, normal) > 0) && (dotProduct(u - A, normal) <= 0))) {
            flag = false;
        }
    }

    Vector rotate90Degrees(const Vector& p, const Vector& origin) {
        double x = p[0];
        double y = p[1];
        double dx = origin[0];
        double dy = origin[1];
        double xx = -(y - dy) + dx;
        double yy = (x - dx) + dy;
        return Vector(xx, yy, 0);
    }
}

Polygon clipPolygon(const Polygon& poly, const Vector& referencePoint, const Vector& edgeStart, const Vector& edgeEnd) {
    Polygon result;
    int numVertices = poly.vertices.size();
    std::vector<Vector> positions(numVertices);
    std::vector<bool> inside(numVertices);

    for (int i = 0; i < numVertices; ++i) {
        bool flag;
        Vector intersectPoint;
        findIntersection(poly.vertices[i], referencePoint, edgeStart, edgeEnd, intersectPoint, flag);
        inside[i] = flag;

        int next = (i + 1) % numVertices;
        findIntersection(poly.vertices[i], poly.vertices[next], edgeStart, edgeEnd, intersectPoint, flag);
        positions[i] = intersectPoint;
    }

    for (int i = 0; i < numVertices; ++i) {
        if (inside[i]) {
            result.vertices.push_back(poly.vertices[i]);
        }
        if (inside[i] != inside[(i + 1) % numVertices]) {
            result.vertices.push_back(positions[i]);
        }
    }

    return result;
}

Polygon createVoronoiCell(const Vector* points, int totalPoints, int index) {
    Polygon cell;
    cell.vertices.push_back(Vector(0, 0, 0));
    cell.vertices.push_back(Vector(0, 1, 0));
    cell.vertices.push_back(Vector(1, 1, 0));
    cell.vertices.push_back(Vector(1, 0, 0));

    for (int i = 0; i < totalPoints; ++i) {
        if (i == index) continue;

        Vector midPoint = (points[i] + points[index]) * 0.5;
        Vector perpendicular = rotate90Degrees(points[i], midPoint);
        cell = clipPolygon(cell, points[index], midPoint, perpendicular);
    }

    return cell;
}
