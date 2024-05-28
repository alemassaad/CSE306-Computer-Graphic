#ifndef VORONOI_H
#define VORONOI_H

#include <vector>
#include "Polygon.h"

std::vector<Polygon> computeVoronoiDiagram(int numPoints, Vector* points);

#endif // VORONOI_H
