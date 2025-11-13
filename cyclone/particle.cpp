//
// Created by arija.hartel on 9/1/2025.
//

#include <cassert>
#include "particle.h"

using namespace cyclone;

void Particle::integrate(real deltaTime) {
    assert(deltaTime > 0.0);

    // Update Pos
    position.addScaledVector(velocity, deltaTime);

    // Acceleration from Force
    Vector3 resultAcceleration = acceleration;
    resultAcceleration += forceAccum * inverseMass;

    // Velocity
    velocity.addScaledVector(resultAcceleration, deltaTime);

    // Drag
    velocity *= real_pow(damping, deltaTime);

    // Clear force (?)
    clearForce();

}
