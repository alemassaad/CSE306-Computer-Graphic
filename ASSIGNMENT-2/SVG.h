#ifndef SVG_H
#define SVG_H

#include <vector>
#include <string>
#include "Polygon.h"

void save_svg(const std::vector<Polygon>& polygons, const std::string& filename, const std::string& fillcol = "none");
void save_svg_animated(const std::vector<Polygon>& polygons, const std::string& filename, int frameid, int nbframes);

#endif // SVG_H
