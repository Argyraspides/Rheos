#include "controller.h"
#include "view.h"
#include "model.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <mutex>

Controller::Controller(Model *model)
{
    m_model = model;
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

const std::vector<Particle>& Controller::RetrieveModel_GetParticles()
{
    return m_model->m_particles;
}
