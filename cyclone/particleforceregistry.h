//
// Created by arija.hartel on 9/11/2025.
//

#ifndef PARTICLEFORCEREGISTRY_H
#define PARTICLEFORCEREGISTRY_H

#pragma once
#include <vector>

#include "cyclone.h"

namespace cyclone {

    class ParticleForceRegistry {
    protected:
        struct ParticleForceRegistration {
            Particle* particle;
            ParticleForceGenerator *fg;
        };

        // registry to update forces, which in turn updates particles
        typedef std::vector<ParticleForceRegistration> Registry;
        Registry registry;

    public:
        // Add force generator to registry
        void add(Particle* particle, ParticleForceGenerator *fg);
        // delete if needed
        void remove(Particle* particle, ParticleForceGenerator *fg);
        // Update all force generators
        void updateForces(real duration);

    };

} // cyclone

#endif //PARTICLEFORCEREGISTRY_H