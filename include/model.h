#pragma once
#include <pthread.h>
#include "particle.h"
#include "cartesian.h"
#include <vector>

class Model
{
public:
    Model();
    ~Model();
    void run(); // BEGINS THE MODEL
    void update();

    // Calculates how much influence a particle has at a point 'dist' units away from its center, given
    // its circle of influence is 'rad'.
    float smoothingKernel(float rad, float dist);

    // Gets the density of the fluid at point 'p'
    float getDensity(const Point &p);

    // Gets a scalar property of the fluid at point 'p'
    float getScalarProperty(const Point &p, int propertyIdx);

    // Returns the derivative of the smoothingKernel() function
    float getSmoothedKernelDerivative(float rad, float dist);

    // Resolves collisions of the particle against the wall
    void handleWallCollisions();

    // Returns the net gradient at a point in space
    Point getNetGradient(const Point &v, const int &propertyIndex);

    // Returns the pressure gradient at a point in space
    Point getPressureGradient(Particle &pC);

    // Converts density to pressure
    float densityToPressure(float density);

    // Gets the average pressure of two densities
    float getSharedPressure(float d1, float d2);

    // Gets the grid ID of a particle based on its location
    int getGridID(const Point &location);

    // Initializes the particle grid
    void initializeGrid();

    // Updates the particle grid
    void updateGrid();

    std::vector<Particle *> &getInRangeParticles(const Particle &p);

    void updateSpatialLookup();

    iPoint posToCellCoo(const Point &p);

    unsigned int hashCell(iPoint cell);

    unsigned int getKeyFromHash(unsigned int hash);

    float m_gravity = 0;
    std::vector<Particle> m_particles;

    // Emscripten doesn't support std::thread for multithreading, only C-type pthread's. This will essentially be a pointer
    // to the run() function so we can actually pass it into pthread_create() in main.cpp
    static void *threadEntry(void *instance)
    {
        reinterpret_cast<Model *>(instance)->run();
        return nullptr;
    }

    bool m_isRunning = true;
    int ENGINE_POLLING_RATE = 60;
    float ENGINE_TIME_STEP = 1.0f / (float)ENGINE_POLLING_RATE;
    double m_time;

    int xBounds = 1920, yBounds = 1080;
};
