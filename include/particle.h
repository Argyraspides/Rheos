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
    float density;
    static float mass;  // Mass of the particle
    static float smRad; // Radial influence a particle has (smoothing radius)
    static float rad;   // Size of the particle
    static float restingDensity;
    static float pressureScaler;
};
