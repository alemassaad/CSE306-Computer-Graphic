#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include "Vector.h"
#include "Polygon.h"
#include "Voronoi.h"
#include "SVG.h"

int main() {
    srand(time(0)); // use current time as our seed for randomness generator

    int numPoints = 400;
    std::vector<Vector> points(numPoints);
    for (int i = 0; i < numPoints; ++i) {
        points[i][0] = rand() / static_cast<double>(RAND_MAX);
        points[i][1] = rand() / static_cast<double>(RAND_MAX);
        points[i][2] = 0;
    }

    auto diagram = computeVoronoiDiagram(numPoints, points.data());
    save_svg(diagram, "diagram.svg");

    return 0;
}
