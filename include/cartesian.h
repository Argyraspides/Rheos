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


struct iPoint
{
    int x = 0, y = 0, z = 0;
    iPoint(int x = 0, int y = 0, int z = 0) : x(x), y(y), z(z) {}

    iPoint operator+(const iPoint &iPoint) const;
    iPoint operator-(const iPoint &iPoint) const;
    iPoint operator/(const int &num) const;
    iPoint operator*(const int &num) const;

    void operator=(const iPoint &iPoint);

    bool operator!=(const iPoint &p) const;
    bool operator==(const iPoint &p) const;

    void normalize();

    int magnitude() const;
};