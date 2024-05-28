#include "Voronoi.h"

std::vector<Polygon> computeVoronoiDiagram(int numPoints, Vector* points) {
    std::vector<Polygon> diagram;
    for (int i = 0; i < numPoints; ++i) {
        Polygon cell = createVoronoiCell(points, numPoints, i);
        diagram.push_back(cell);
    }
    return diagram;
}
