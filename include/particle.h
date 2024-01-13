#pragma once
#include "cartesian.h"

class Particle
{
public:
    Particle(Point pos, Point vel) : pos(pos), vel(vel){};
    Particle(Point pos, Point predPos, Point vel) : pos(pos), predPos(predPos), vel(vel){};
    Point pos; // Position of the particle
    Point vel; // Velocity of the particle
    Point predPos;

    static const int densityIdx = 0, pressureIdx = 1;
    float prop[2] = {0.0f, 0.0f}; // 0 = density, 1 = pressure

    static float mass;  // Mass of the particle
    static float smRad; // Radial influence a particle has (smoothing radius)
    static float rad;   // Size of the particle
    static float restingDensity;
    static float pressureScaler;
};
