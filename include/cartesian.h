#pragma once
#include <cmath>

struct Point
{
    float x = 0, y = 0, z = 0;
    Point(float x = 0, float y = 0, float z = 0) : x(x), y(y), z(z) {}

    Point operator+(const Point &point) const;
    Point operator-(const Point &point) const;
    Point operator/(const float &num) const;
    Point operator*(const float &num) const;

    void operator=(const Point &point);

    bool operator!=(const Point &p) const;
    bool operator==(const Point &p) const;

    void normalize();

    float magnitude() const;
};