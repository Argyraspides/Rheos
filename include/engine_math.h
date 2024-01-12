#pragma once

#include "cartesian.h"

namespace Math
{
    static Point getNormal2D(Point p);
    static double dotProd(const Point &p1, const Point &p2);
    static double dist(const Point &p1, const Point &p2);

    // Origin of the 3D cartesian space (x, y, and z = 0)
    static Point origin = {0, 0, 0};
    static Point defaultPt = origin;


    // Calculates the distance between two points
    static double dist(const Point &p1, const Point &p2)
    {
        return sqrt(
            pow((p1.x - p2.x), 2) +
            pow((p1.y - p2.y), 2) +
            pow((p1.z - p2.z), 2));
    }

}
