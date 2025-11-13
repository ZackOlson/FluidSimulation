//
// Created by cocoa on 19/09/2025.
//

#include "particlespring.h"

namespace cyclone {
    void ParticleSpring::updateForce(Particle *particle, real duration) {
        // Calculate Force
        Vector3 force;
        force = particle->getPosition();
        force -= other->getPosition();

        // Calculate Magnitude
        real magnitude = force.magnitude();
        magnitude = real_abs(magnitude - restLength);
        magnitude *= springConstant;

        // Final Force
        force.normalise();
        force *= -magnitude;
        particle->addForce(force);

    }


}