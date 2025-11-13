//
// Created by arija.hartel on 9/1/2025.
// Adding onto the Cyclone library
// References the book when questions arise
//

#pragma once
#ifndef PARTICLE_H
#define PARTICLE_H
#include "core.h"

namespace cyclone {
    class Particle {
    public:

        // Helpers
        Vector3 getPosition() { return position; };
        Vector3 getVelocity() { return velocity; };
        Vector3 getAcceleration() { return acceleration; };
        real getDamping() { return damping; };
        real getMass() { return inverseMass; };
        Vector3 getTotalAccum() { return forceAccum; };

        // Setters
        void setPosition(Vector3 newPos) { position = newPos; };
        void setPosition(float x, float y, float z) { position = Vector3(x, y, z); };
        void setVelocity(Vector3 newVel) { velocity = newVel; };
        void setVelocity(float x, float y, float z) { velocity = Vector3(x, y, z); };
        void setAcceleration(Vector3 newAcc) { acceleration = newAcc; };
        void setAcceleration(float x, float y, float z) { acceleration = Vector3(x, y, z); };
        void setDamping(real newDamping) { damping = newDamping; };
        // Forces
        void setMass(real newMass) { inverseMass = newMass; };
        void setTotalAccum(Vector3 newAcc) { forceAccum = newAcc; };
        void addForce(Vector3 force) { forceAccum += force; };
        void clearForce() { forceAccum = Vector3(); };

        // Funky guys
        void integrate(real deltaTime);


        // These are a maybe
        //void launchParticle(real initialSpeed, real angle, real startPosition);
        //void updateParticle(real deltaTime); // might use


    protected:
        // linear position
        Vector3 position;
        // linear velocity
        Vector3 velocity;
        // constant acceleration, like gravity
        Vector3 acceleration;
        real damping;

        // Forces
        real inverseMass;
        Vector3 forceAccum;

    };
}

#endif //PARTICLE_H
