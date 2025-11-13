//
// Created by cocoa on 19/09/2025.
//

#ifndef PARTICLESPRING_H
#define PARTICLESPRING_H
#pragma once
#include "core.h"
#include "particle.h"
#include "particleforcegenerator.h"

namespace cyclone {

    class ParticleSpring : public ParticleForceGenerator {

        // Other half of spring
        Particle *other;

        // Anchor
        //Vector3 *anchor;

        // Spring Constant
        real springConstant;

        // Rest length
        real restLength;

    public:
        ParticleSpring(Particle *other, real springConstant, real restLength) : other(other), springConstant(springConstant), restLength(restLength) {}; // : other(other), springConstant(springConstant), restLength(restLength) {}

        virtual void updateForce(Particle *particle, real duration);

    };

}



#endif //PARTICLESPRING_H
