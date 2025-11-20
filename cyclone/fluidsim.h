//
// Created by cocoa on 20/11/2025.
// Edited by zack on 20/11/2025
//

#ifndef FLUIDSIM_H
#define FLUIDSIM_H
#include "particle.h"
#include "particleforcegenerator.h"
#include <vector>

namespace cyclone {
    class FluidSim : public ParticleForceGenerator {
    public:
        // constructor
        FluidSim(int gridW, int gridH, int gridD, real cellSize = 4.0f, real density = 1000.0f);

        // simulation step
        void step(std::vector<Particle>& particles, real dt);

        virtual void updateForce(Particle* p, real duration) override;

    protected:
        // Core steps
        void particlesToGrid(const std::vector<Particle>& particles);
        void applyGridForces(real dt);
        void saveGridVelocities();
        void solvePressure(real dt);
        void gridToParticles(std::vector<Particle>& particles, real flipRatio);
        void advectParticles(std::vector<Particle>& particles, real dt);
        void handleParticleCollision(Particle& p);

        inline int idx3(int x, int y, int z) const 
        {
            return x + y * gridWidth + z * gridWidth * gridHeight;
        }
        void clearGrid();

        // Grid data
        int gridWidth, gridHeight, gridDepth;
        real cellSize;
        real invCellSize;

        std::vector<real> u, v, w;
        std::vector<real> uPrev, vPrev, wPrev;
        std::vector<real> mass;
        std::vector<real> pressure;

        // parameters
        real gravity{ -400.0f };
        real density;
        int pressureIterations{ 30 };
        // bounciness
        real boundaryDamping{ -0.5f };
    };
}

#endif //FLUIDSIM_H
