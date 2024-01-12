#pragma once
#include "imgui.h"
#include "model.h"

class Controller
{
public:
    Model *m_model;
    Controller(Model *model);

    void UpdateModel_ChangePressureScaler(float pressure);
    void UpdateModel_ChangeRestingDensity(float density);
    void UpdateModel_ChangeGravity(float gravity);
    void UpdateModel_ChangeRadialInfluence(float rad);

    const float RetrieveModel_GetPressureScaler();
    const float RetrieveModel_GetRestingDensity();
    const float RetrieveModel_GetGravity();
    const float RetrieveModel_GetRadialInfluence();

    const std::vector<Particle> &RetrieveModel_GetParticles();
};