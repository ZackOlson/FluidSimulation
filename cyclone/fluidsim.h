//
// Created by cocoa on 20/11/2025.
//

#ifndef FLUIDSIM_H
#define FLUIDSIM_H
#include "particle.h"
#include "particleforcegenerator.h"

namespace cyclone {
    class FluidSim : public ParticleForceGenerator {
    public:
        FluidSim(real mDensity) : density(mDensity) {}; // Set variables that we need

        virtual void updateForce(Particle* p, real duration);

        // Not sure what variables are needed exactly yet but we have a plan
        void updateParticleDensities(Particle* p, real duration);
        void solveIncompressibilties(Particle* p, real duration);
        void transferVelocities(Particle* p, real duration);
        void handleParticleCollisions(Particle* p, real duration);

    protected:
        // Density of everything
        real density;

    };
}

#endif //FLUIDSIM_H
