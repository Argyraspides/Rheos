#include "controller.h"
#include "view.h"
#include "model.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <mutex>
#include "engine_math.h"

Controller::Controller(Model *model)
{
    m_model = model;
}

void Controller::ShutModel()
{
    m_model->m_isRunning = false;
}

void Controller::UpdateModel_ChangePressureScaler(float pressure)
{
    Particle::pressureScaler = pressure;
}

void Controller::UpdateModel_ChangeRestingDensity(float density)
{
    Particle::restingDensity = density;
}

void Controller::UpdateModel_ChangeGravity(float gravity)
{
    m_model->m_gravity = gravity;
}

void Controller::UpdateModel_ChangeRadialInfluence(float rad)
{
    Particle::smRad = rad;
}

void Controller::UpdateModel_ApplyFluidForce(const Point &loc, float radInfluence, float forceStrength)
{
    Point force = {0,0,0};
    for (int i = 0; i < m_model->m_particles.size(); i++)
    {
        float d = Math::dist(loc, m_model->m_particles[i].pos);
        if (d < radInfluence)
        {
            Point offset = loc - m_model->m_particles[i].pos;
            Point dir = d < std::numeric_limits<float>::min() ? Math::origin : offset / d;
            float centerT = 1-d/radInfluence;
            force = force + (dir * forceStrength - m_model->m_particles[i].vel) * centerT;
            m_model->applyForce(m_model->m_particles[i], force);
        }
        force = {0,0,0};
    }
}

const float Controller::RetrieveModel_GetPressureScaler()
{
    return Particle::pressureScaler;
}

const float Controller::RetrieveModel_GetRestingDensity()
{
    return Particle::restingDensity;
}

const float Controller::RetrieveModel_GetGravity()
{
    return m_model->m_gravity;
}

const float Controller::RetrieveModel_GetRadialInfluence()
{
    return Particle::smRad;
}

const std::vector<Particle> &Controller::RetrieveModel_GetParticles()
{
    return m_model->m_particles;
}
