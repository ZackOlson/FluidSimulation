//
// Created by arija.hartel on 9/11/2025.
//

#include "particleforceregistry.h"

namespace cyclone {
    void ParticleForceRegistry::add(Particle* particle, ParticleForceGenerator *fg) {
        auto *registration = new ParticleForceRegistration();
        registration->particle = particle;
        registration->fg = fg;
        registry.push_back(*registration);

    }

    // Delete particle from list, this might not work
    void ParticleForceRegistry::remove(Particle* particle, ParticleForceGenerator *fg) {
        auto *temp = new ParticleForceRegistration();
        temp->particle = particle;
        temp->fg = fg;
        for (auto it = registry.begin(); it != registry.end(); it++) {
            if (temp->particle == it->particle) {
                registry.erase(it);
                return;
            }
        }

    }

    // Update all force generators
    void ParticleForceRegistry::updateForces(real duration) {
        // go through all in registry
        // use update force function, which updates their particle
        auto i = registry.begin();
        for (; i != registry.end(); i++) {
            i->fg->updateForce(i->particle, duration);
        }
    }
} // cyclone