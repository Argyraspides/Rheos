#include "cartesian.h"

Point Point::operator+(const Point &point) const
{
    return {
        x + point.x,
        y + point.y,
        z + point.z};
}

Point Point::operator-(const Point &point) const
{
    return {
        x - point.x,
        y - point.y,
        z - point.z};
}

void Point::operator=(const Point &point)
{
    x = point.x;
    y = point.y;
    z = point.z;
}

Point Point::operator/(const float &num) const
{
    return {
        x / num,
        y / num,
        z / num};
}

Point Point::operator*(const float &num) const
{
    return {
        x * num,
        y * num,
        z * num};
}

bool Point::operator!=(const Point &p) const
{
    return (x != p.x && y != p.y && z != p.z);
}

bool Point::operator==(const Point &p) const
{
    return (x == p.x && y == p.y && z == p.z);
}

void Point::normalize()
{
    float len = sqrt(x * x + y * y + z * z);
    if (len == 0)
        return;

    x /= len;
    y /= len;
    z /= len;
}

float Point::magnitude() const
{
    return sqrt(x * x + y * y + z * z);
}
