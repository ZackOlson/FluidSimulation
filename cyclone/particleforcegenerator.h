//
// Created by arija.hartel on 9/11/2025.
//

#ifndef PARTICLEFORCEGENERATOR_H
#define PARTICLEFORCEGENERATOR_H
#pragma once

#include "core.h"

namespace cyclone {
    // Abstract Base Class
    class ParticleForceGenerator {
    public:
        virtual void updateForce(Particle *particle, real duration) = 0;

    };

} // cyclone

#endif //PARTICLEFORCEGENERATOR_H
