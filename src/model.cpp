#include "model.h"
#include "engine_math.h"
#include <chrono>
#include <iostream>
#include <algorithm>

Model::Model()
{
    int particleCount = 3000;
    int maxColumns = 70;

    float separation = Particle::rad * 2;
    m_gravity = 300;

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
        m_particles.push_back(Particle({xPos, yPos}, {xPos, yPos}, {0, 0}));
    }

    initializeGrid();
}

Model::~Model()
{
}

// Begins the engine
void Model::run()
{

    std::srand(static_cast<unsigned>(std::time(nullptr)));

    auto start = std::chrono::high_resolution_clock::now();

    while (m_isRunning)
    {
        auto end = std::chrono::high_resolution_clock::now();
        double duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        if (duration >= ENGINE_TIME_STEP * 10.0f)
        {
            m_time += ENGINE_TIME_STEP;
            update();
            start = std::chrono::high_resolution_clock::now();
        }
    }
}

void Model::update()
{

    for (Particle &p : m_particles)
    {
        p.vel.y += m_gravity * ENGINE_TIME_STEP;
        p.predPos = p.pos + p.vel * (1.0f / 500.0f);
    }

    handlePredWallCollisions();
    handleWallCollisions();

    updateGrid();

    for (Particle &p : m_particles)
    {
        p.density = getDensity(p.predPos);
    }

    for (Particle &p : m_particles)
    {
        Point pressureForce = getPressureGradient(p);
        Point pressureAcceleration = pressureForce / p.density;
        if (p.density > 0)
            p.vel = p.vel + pressureAcceleration * ENGINE_TIME_STEP;
    }

    for (Particle &p : m_particles)
    {
        p.pos = p.pos + p.vel * ENGINE_TIME_STEP;
    }
}

// Calculates the influence a particle has at a point 'dist' away from the particles center,
// given that the particle has a radius of influence 'rad'/
float Model::smoothingKernel(float rad, float dist)
{
    // Smoothing kernel function:
    // y = (r - x)^2 / (πr^4 / 6)
    if (dist > rad)
        return 0;
    float vol = (M_PI * rad * rad * rad * rad) * 0.1666666666666f;
    float delta = rad - dist;
    return (delta * delta) / vol;
}

// Get the density at a particle 'p'
float Model::getDensity(const Point &p)
{
    float density = 0.0f;
    std::vector<iPoint> inRangeCells = getInRangeCells(p);
    for (int c = 0; c < inRangeCells.size(); c++)
    {
        int iterations = m_particleGridSizes[inRangeCells[c].x][inRangeCells[c].y];
        for (int i = 0; i < iterations; i++)
        {
            float dist = Math::dist(p, m_particleGrid[inRangeCells[c].x][inRangeCells[c].y][i]->pos);
            density += Particle::mass * smoothingKernel(Particle::smRad, dist);
        }
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

// Gets the derivative of the smoothed kernel function
float Model::getSmoothedKernelDerivative(float rad, float dist)
{
    if (dist >= rad)
        return 0;

    float scale = 12.0f / (rad * rad * rad * rad * M_PI);
    return (dist - rad) * scale;
}

// Gets the pressure gradient of a particle
Point Model::getPressureGradient(Particle &pC)
{
    Point pressureGrad = {0.0f, 0.0f, 0.0f};

    std::vector<iPoint> inRangeCells = getInRangeCells(pC.predPos);
    for (int c = 0; c < inRangeCells.size(); c++)
    {
        int iterations = m_particleGridSizes[inRangeCells[c].x][inRangeCells[c].y];
        for (int i = 0; i < iterations; i++)
        {
            // Get the distance from the current particle to the point in space we want to find the gradient of
            Particle *_pC = m_particleGrid[inRangeCells[c].x][inRangeCells[c].y][i];
            Point dir = (_pC->predPos - pC.predPos);
            float dist = dir.magnitude();

            // If the particles are directly on top of each other we can just choose a random direction
            if (dist == 0)
            {
                dir = {static_cast<float>(std::rand()) / 1.0f + 0.01f, static_cast<float>(std::rand()) / 1.0f + 0.01f};
                dist = dir.magnitude();
            }

            dir.normalize();

            // The rate of change of the property at this point
            float slope = getSmoothedKernelDerivative(Particle::smRad, dist);
            // Density of the current particle
            float density = _pC->density;

            float sharedP = getSharedPressure(density, pC.density);

            if (density > 0)
                pressureGrad = pressureGrad - dir * sharedP * slope * Particle::mass / density;
        }
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

void Model::handlePredWallCollisions()
{
    for (Particle &particle : m_particles)
    {
        if (xBounds - particle.predPos.x <= Particle::rad)
        {
            particle.predPos.x = xBounds - Particle::rad;
        }
        if (particle.predPos.x <= Particle::rad)
        {
            particle.predPos.x = Particle::rad;
        }
        if (yBounds - particle.predPos.y <= Particle::rad)
        {
            particle.predPos.y = yBounds - Particle::rad;
        }
        if (particle.predPos.y <= Particle::rad)
        {
            particle.predPos.y = Particle::rad;
        }
    }
}

iPoint Model::getGridCoo(const Point &loc)
{
    return iPoint(loc.y / Particle::smRad, loc.x / Particle::smRad, 0);
}

std::vector<iPoint> Model::getInRangeCells(const Point &p)
{

    static int rows = yBounds / Particle::smRad;
    static int cols = xBounds / Particle::smRad;

    static std::vector<iPoint> cells(9, {0, 0, 0});
    iPoint ctr = getGridCoo(p);
    iPoint begin = {Math::positiveMod(ctr.x - 1, rows), Math::positiveMod(ctr.y - 1, cols)};

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cells[j + i * 3] = {(begin.x + i) % rows, (begin.y + j) % cols};
        }
    }

    return cells;
}

void Model::applyForce(Particle &p, Point dir)
{

}

void Model::updateGrid()
{

    static int rows = yBounds / Particle::smRad;
    static int cols = xBounds / Particle::smRad;

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            m_particleGridSizes[i][j] = 0;
        }
    }

    for (int i = 0; i < m_particles.size(); i++)
    {
        iPoint gc = getGridCoo(m_particles[i].predPos);
        Particle *ptr = &m_particles[i];
        m_particleGrid[gc.x][gc.y][m_particleGridSizes[gc.x][gc.y]++] = ptr;
    }
}

void Model::initializeGrid()
{
    static int rows = yBounds / Particle::smRad;
    static int cols = xBounds / Particle::smRad;

    for (int i = 0; i < rows; i++)
    {
        m_particleGridSizes.push_back(std::vector<int>(cols, 0));
    }

    for (int i = 0; i < rows; i++)
    {
        m_particleGrid.push_back(std::vector<std::vector<Particle *>>());
        for (int j = 0; j < cols; j++)
        {
            m_particleGrid[i].push_back(std::vector<Particle *>());
            for (int k = 0; k < m_particles.size(); k++)
            {
                m_particleGrid[i][j].push_back(nullptr);
            }
        }
    }
}
