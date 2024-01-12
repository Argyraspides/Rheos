#include "model.h"
#include "engine_math.h"
#include "omp.h"
#include <chrono>
#include <iostream>
#include <algorithm>
#include <cmath>

Model::Model()
{
    int particleCount = 4000;
    int maxColumns = 50;

    float separation = Particle::rad * 2;

    int width = maxColumns * separation;
    int height = (particleCount / maxColumns) * separation;

    int xCenterOffset = (xBounds - width) / 2;
    int yCenterOffset = (yBounds - height) / 2;

    m_particles.reserve(particleCount);

    // Initialize all the particles in a grid
    for (int i = 0; i < particleCount; i++)
    {
        float xPos = (i % maxColumns) * separation + xCenterOffset;
        float yPos = (i / maxColumns) * separation + yCenterOffset;
        m_particles.push_back(Particle({xPos, yPos}, {0, 0}));
    }
}

Model::~Model()
{
}

// Begins the engine
void Model::run()
{

    auto start = std::chrono::high_resolution_clock::now();

    while (m_isRunning)
    {
        auto end = std::chrono::high_resolution_clock::now();
        double duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        if (duration >= ENGINE_TIME_STEP * 1000.0f)
        {
            m_time += ENGINE_TIME_STEP;
            update();
            start = std::chrono::high_resolution_clock::now();
        }
    }
}

void Model::update()
{

#pragma omp parallel for
    for (Particle &p : m_particles)
    {
        p.vel.y += m_gravity * ENGINE_TIME_STEP;
        p.prop[Particle::densityIdx] = getDensity(p.pos);
    }

#pragma omp parallel for
    for (Particle &p : m_particles)
    {
        Point pressureForce = getPressureGradient(p);
        Point pressureAcceleration = pressureForce / p.prop[Particle::densityIdx];
        p.vel = p.vel + pressureAcceleration * ENGINE_TIME_STEP;
    }

#pragma omp parallel for
    for (Particle &p : m_particles)
    {
        p.pos = p.pos + p.vel * ENGINE_TIME_STEP;
    }

    handleWallCollisions();
}

// Calculates the influence a particle has at a point 'dist' away from the particles center,
// given that the particle has a radius of influence 'rad'/
float Model::smoothingKernel(float rad, float dist)
{
    // Smoothing kernel function:
    // y = (r - x)^2 / (πr^4 / 6)
    if (dist > rad)
        return 0;
    float vol = (M_PI * pow(rad, 4.0f)) / 6.0f;
    float delta = rad - dist;
    return (delta * delta) / vol;
}

// Get the density of the fluid at a point 'p'
float Model::getDensity(const Point &p)
{
    float density = 0.0f;

    for (const Particle &particle : m_particles)
    {
        float dist = Math::dist(particle.pos, p);
        if (dist > Particle::smRad)
            continue;
        density += Particle::mass * smoothingKernel(Particle::smRad, dist);
    }

    return density;
}

// Converts density to pressure
float Model::densityToPressure(float density)
{
    float dX = density - Particle::restingDensity;
    float pressure = dX * Particle::pressureScaler;
    return pressure;
}

// Gets an arbitrary scalar property (density, pressure, etc)
float Model::getScalarProperty(const Point &p, int propertyIdx)
{
    float A = 0.0f;

    for (int i = 0; i < m_particles.size(); i++)
    {
        float dist = Math::dist(m_particles[i].pos, p);
        float influence = smoothingKernel(dist, Particle::smRad);
        float density = m_particles[i].prop[propertyIdx];
        A += m_particles[i].prop[propertyIdx] * Particle::mass / density * influence;
    }

    return A;
}

// Gets the derivative of the smoothed kernel function
float Model::getSmoothedKernelDerivative(float rad, float dist)
{
    if (dist >= rad)
        return 0;

    float scale = 12.0f / (pow(rad, 4) * M_PI);
    return (dist - rad) * scale;
}

// Gets the net gradient at a particular point in space, of a particular scalar property
Point Model::getNetGradient(const Point &v, const int &propertyIndex)
{
    Point totalGradient = {0, 0};
    for (int i = 0; i < m_particles.size(); i++)
    {
        Point direction = m_particles[i].pos - v;
        // If two particles are directly on top of each other, or we are comparing a particle to itself,
        // we don't want divide by zero errors. In this case we can just take a random direction
        if (direction == Point(0, 0))
            direction = {1, 1};

        float dist = direction.magnitude();
        direction.normalize();

        float gradient = getSmoothedKernelDerivative(Particle::smRad, dist);
        // Get the current density stored in the particle
        float density = m_particles[i].prop[Particle::densityIdx];

        Point mult = direction * m_particles[i].prop[propertyIndex] * gradient * Particle::mass;
        totalGradient = totalGradient + mult / density;
    }
    return totalGradient;
}

// Gets the pressure gradient of a particle
Point Model::getPressureGradient(Particle &pC)
{
    Point pressureGrad = {0.0f, 0.0f, 0.0f};
    for (const Particle &p : m_particles)
    {
        // Get the distance from the current particle to the point in space we want to find the gradient of
        float dist = (p.pos - pC.pos).magnitude();
        Point dir;

        if (dist == 0)
        {
            dir = {static_cast<float>(rand()) / 1.0f, static_cast<float>(rand()) / 1.0f};
            dist = dir.magnitude();
        }
        else
        {
            // The actual direction vector
            dir = (p.pos - pC.pos) / dist;
        }

        // The rate of change of the property at this point
        float slope = getSmoothedKernelDerivative(Particle::smRad, dist);
        // Density of the current particle
        float density = p.prop[Particle::densityIdx];

        float sharedP = getSharedPressure(density, pC.prop[Particle::densityIdx]);
        pressureGrad = pressureGrad - dir * sharedP * slope * Particle::mass / p.prop[Particle::densityIdx];
    }
    return pressureGrad;
}

// Gets the average pressure of two density regions
float Model::getSharedPressure(float d1, float d2)
{
    float pA = densityToPressure(d1);
    float pB = densityToPressure(d2);
    return (pA + pB) / 2.0f;
}

// Handles particle wall collisions
void Model::handleWallCollisions()
{
    static float collisionDamping = .75f;
    for (Particle &particle : m_particles)
    {
        if (xBounds - particle.pos.x <= Particle::rad)
        {
            particle.pos.x = xBounds - Particle::rad;
            particle.vel.x *= -1.0f * collisionDamping;
        }
        if (particle.pos.x <= Particle::rad)
        {
            particle.pos.x = Particle::rad;
            particle.vel.x *= -1.0f * collisionDamping;
        }
        if (yBounds - particle.pos.y <= Particle::rad)
        {
            particle.pos.y = yBounds - Particle::rad;
            particle.vel.y *= -1.0f * collisionDamping;
        }
        if (particle.pos.y <= Particle::rad)
        {
            particle.pos.y = Particle::rad;
            particle.vel.y *= -1.0f * collisionDamping;
        }
    }
}
