//
// Created by cocoa on 20/11/2025.
// Edited by zack on 20/11/2025.
//

#include "fluidsim.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>

using namespace cyclone;

// using grid mehtod (poorly lol)
FluidSim::FluidSim(int gridW, int gridH, int gridD, real cSize, real dens)
    : gridWidth(gridW), gridHeight(gridH), gridDepth(gridD), cellSize(cSize), density(dens)
{

    invCellSize = 1.0f / cellSize;
    int n = gridWidth * gridHeight * gridDepth;

    u.assign(n, (real)0.0);
    v.assign(n, (real)0.0);
    w.assign(n, (real)0.0);

    uPrev.assign(n, (real)0.0);
    vPrev.assign(n, (real)0.0);
    wPrev.assign(n, (real)0.0);

    mass.assign(n, (real)0.0);
    pressure.assign(n, (real)0.0);
}

void FluidSim::updateForce(Particle* p, real duration) {
    Vector3 force;
}

void FluidSim::clearGrid() {
    int N = gridWidth * gridHeight * gridDepth;

    std::fill(u.begin(), u.end(), 0.0f);
    std::fill(v.begin(), v.end(), 0.0f);
    std::fill(w.begin(), w.end(), 0.0f);

    std::fill(mass.begin(), mass.end(), 0.0f);
}

// Transfers from particles to grid (velocity and mass)
void FluidSim::particlesToGrid(const std::vector<Particle>& particles)
{
    clearGrid();

    for (const auto& p : particles)
    {
        Vector3 pos = p.getPosition();

        float px = pos.x * invCellSize;
        float py = pos.y * invCellSize;
        float pz = pos.z * invCellSize;

        int ix = (int)floor(px);
        int iy = (int)floor(py);
        int iz = (int)floor(pz);

        float fx = px - ix;
        float fy = py - iy;
        float fz = pz - iz;

        float wx[2] = { 1.0f - fx, fx };
        float wy[2] = { 1.0f - fy, fy };
        float wz[2] = { 1.0f - fz, fz };

        Vector3 vel = p.getVelocity();

        for (int k = 0; k < 2; k++)
        {
            int gz = iz + k;
            if (gz < 0 || gz >= gridDepth) continue;

            for (int j = 0; j < 2; j++)
            {
                int gy = iy + j;
                if (gy < 0 || gy >= gridHeight) continue;

                for (int i = 0; i < 2; i++)
                {
                    int gx = ix + i;
                    if (gx < 0 || gx >= gridWidth) continue;

                    float wght = wx[i] * wy[j] * wz[k];
                    int id = idx3(gx, gy, gz);

                    mass[id] += wght;
                    u[id] += wght * vel.x;
                    v[id] += wght * vel.y;
                    w[id] += wght * vel.z;
                }
            }
        }
    }

    // Normalize
    int N = gridWidth * gridHeight * gridDepth;
    for (int i = 0; i < N; i++)
    {
        if (mass[i] > 0.0f) {
            u[i] /= mass[i];
            v[i] /= mass[i];
            w[i] /= mass[i];
        }
    }
}

// Apply forces to the grid (gravity only for the moment, add mouse movement thing later)
void FluidSim::applyGridForces(real dt) {
    
    for (size_t i = 0; i < v.size(); ++i) 
    {
        v[i] += gravity * dt;
    }
}

void FluidSim::saveGridVelocities() {
    uPrev = u;
    vPrev = v;
    wPrev = w;
}

// attempt at a basic Jacobi solver for pressure
// Used this code base as a reference, along with the resources we already looked at 
// https://github.com/blaxill/jacobi/blob/master/examples/fluid.rs
void FluidSim::solvePressure(real dt)
{
    int N = gridWidth * gridHeight * gridDepth;

    std::vector<real> div(N, 0.0f);
    std::vector<real> pnew(N, 0.0f);

    // divergence
    for (int z = 1; z < gridDepth - 1; ++z) 
    {
        for (int y = 1; y < gridHeight - 1; ++y) 
        {
            for (int x = 1; x < gridWidth - 1; ++x) 
            {

                int k = idx3(x, y, z);

                real du_dx = (u[idx3(x + 1, y, z)] - u[idx3(x - 1, y, z)])
                    * 0.5f * invCellSize;

                real dv_dy = (v[idx3(x, y + 1, z)] - v[idx3(x, y - 1, z)])
                    * 0.5f * invCellSize;

                real dw_dz = (w[idx3(x, y, z + 1)] - w[idx3(x, y, z - 1)])
                    * 0.5f * invCellSize;

                div[k] = du_dx + dv_dy + dw_dz;
            }
        }
    }

    float h2inv = 1.0f / (invCellSize * invCellSize);

    // Jacobi iterations
    for (int it = 0; it < pressureIterations; ++it)
    {
        for (int z = 1; z < gridDepth - 1; ++z) 
        {
            for (int y = 1; y < gridHeight - 1; ++y) 
            {
                for (int x = 1; x < gridWidth - 1; ++x) 
                {

                    int k = idx3(x, y, z);

                    real sumNei =
                        pressure[idx3(x - 1, y, z)] +
                        pressure[idx3(x + 1, y, z)] +
                        pressure[idx3(x, y - 1, z)] +
                        pressure[idx3(x, y + 1, z)] +
                        pressure[idx3(x, y, z - 1)] +
                        pressure[idx3(x, y, z + 1)];

                    pnew[k] = (sumNei - div[k] * h2inv) / 6.0f;
                }
            }
        }

        // swap pnew into pressure
        for (int i = 0; i < N; ++i)
            pressure[i] = pnew[i];
    }

    // Subtract pressure
    for (int z = 1; z < gridDepth - 1; ++z) 
    {
        for (int y = 1; y < gridHeight - 1; ++y) 
        {
            for (int x = 1; x < gridWidth - 1; ++x) 
            {

                int k = idx3(x, y, z);

                real gradPx = (pressure[idx3(x + 1, y, z)] -
                    pressure[idx3(x - 1, y, z)]) * 0.5f * invCellSize;

                real gradPy = (pressure[idx3(x, y + 1, z)] -
                    pressure[idx3(x, y - 1, z)]) * 0.5f * invCellSize;

                real gradPz = (pressure[idx3(x, y, z + 1)] -
                    pressure[idx3(x, y, z - 1)]) * 0.5f * invCellSize;

                u[k] -= gradPx;
                v[k] -= gradPy;
                w[k] -= gradPz;
            }
        }
    }

    // boundary
    for (int z = 0; z < gridDepth; ++z)
        for (int x = 0; x < gridWidth; ++x) 
        {
            v[idx3(x, 0, z)] = 0;     v[idx3(x, gridHeight - 1, z)] = 0;
            u[idx3(x, 0, z)] = 0;     u[idx3(x, gridHeight - 1, z)] = 0;
            w[idx3(x, 0, z)] = 0;     w[idx3(x, gridHeight - 1, z)] = 0;
        }

    for (int z = 0; z < gridDepth; ++z)
        for (int y = 0; y < gridHeight; ++y) 
        {
            u[idx3(0, y, z)] = 0;     u[idx3(gridWidth - 1, y, z)] = 0;
            v[idx3(0, y, z)] = 0;     v[idx3(gridWidth - 1, y, z)] = 0;
            w[idx3(0, y, z)] = 0;     w[idx3(gridWidth - 1, y, z)] = 0;
        }

    for (int y = 0; y < gridHeight; ++y)
        for (int x = 0; x < gridWidth; ++x) 
        {
            w[idx3(x, y, 0)] = 0;     w[idx3(x, y, gridDepth - 1)] = 0;
        }
}

// Transfers from grid to particles (velocity and mass)
void FluidSim::gridToParticles(std::vector<Particle>& particles, real flipRatio)
{
    for (auto& p : particles)
    {
        Vector3 pos = p.getPosition();

        float px = pos.x * invCellSize;
        float py = pos.y * invCellSize;
        float pz = pos.z * invCellSize;

        int ix = (int)floor(px);
        int iy = (int)floor(py);
        int iz = (int)floor(pz);

        float fx = px - ix;
        float fy = py - iy;
        float fz = pz - iz;

        float wx[2] = { 1 - fx, fx };
        float wy[2] = { 1 - fy, fy };
        float wz[2] = { 1 - fz, fz };

        Vector3 pic(0, 0, 0);
        Vector3 flipDelta(0, 0, 0);

        for (int k = 0; k < 2; k++) 
        {
            int gz = iz + k;
            if (gz < 0 || gz >= gridDepth) continue;

            for (int j = 0; j < 2; j++) 
            {
                int gy = iy + j;
                if (gy < 0 || gy >= gridHeight) continue;

                for (int i = 0; i < 2; i++) 
                {
                    int gx = ix + i;
                    if (gx < 0 || gx >= gridWidth) continue;

                    float wght = wx[i] * wy[j] * wz[k];
                    int id = idx3(gx, gy, gz);

                    pic.x += u[id] * wght;
                    pic.y += v[id] * wght;
                    pic.z += w[id] * wght;

                    flipDelta.x += (u[id] - uPrev[id]) * wght;
                    flipDelta.y += (v[id] - vPrev[id]) * wght;
                    flipDelta.z += (w[id] - wPrev[id]) * wght;
                }
            }
        }

        Vector3 newVelocity = pic * (1.0f - flipRatio) + (p.getVelocity() + flipDelta) * flipRatio;

        p.setVelocity(newVelocity);
    }
}

void FluidSim::advectParticles(std::vector<Particle>& particles, real dt) {
    /*
    for (auto& p : particles) 
    {

        // Euler advection
        Vector3 pos = p.getPosition();
        Vector3 vel = p.getVelocity();
        pos += vel * dt;
        p.setPosition(pos);
        handleParticleCollision(p);
    } */



    real stiffness = 1.0f;

    // Trying to see if particles can check distance
    for (int i = 0; i < particles.size(); i++) {
        Vector3 tempForce = Vector3(0, 0, 0);
        for (int j = 0; j < particles.size(); j++) {
            if (j != i) { // skip self

                Vector3 direction = Vector3(particles[i].getPosition() - particles[j].getPosition());
                real distance = direction.magnitude();

                float restSize = 1.0f;
                if (distance < restSize) {
                    tempForce += direction.unit() * stiffness * (restSize - distance);
                }

            }
        }

        // Give force!
        Vector3 force = tempForce * dt;
        Vector3 pos = particles[i].getPosition();
        pos += force * gravity * dt;
        particles[i].setPosition(pos);
        handleParticleCollision(particles[i]);


    }

}

// Handles particle collisions (I think, although im not sure if i got it working properly)
void FluidSim::handleParticleCollision(Particle& p)
{
    Vector3 pos = p.getPosition();
    Vector3 vel = p.getVelocity();

    real minX = 1.0f;
    real minY = 1.0f;
    real minZ = 1.0f;

    real maxX = gridWidth * cellSize - 2.0f;
    real maxY = gridHeight * cellSize - 2.0f;
    real maxZ = gridDepth * cellSize - 2.0f;

    if (pos.x < minX) { pos.x = minX; vel.x *= boundaryDamping; }
    if (pos.y < minY) { pos.y = minY; vel.y *= boundaryDamping; }
    if (pos.z < minZ) { pos.z = minZ; vel.z *= boundaryDamping; }

    if (pos.x > maxX) { pos.x = maxX; vel.x *= boundaryDamping; }
    if (pos.y > maxY) { pos.y = maxY; vel.y *= boundaryDamping; }
    if (pos.z > maxZ) { pos.z = maxZ; vel.z *= boundaryDamping; }

    p.setPosition(pos);
    p.setVelocity(vel);
}

void FluidSim::step(std::vector<Particle>& particles, real dt) {

    // Convert particle info to grid
    particlesToGrid(particles);

    // Apply forces to grid
    //applyGridForces(dt); TEST

    // Save grid
    saveGridVelocities();

    // Pressure solve
    solvePressure(dt);

    // convert grid info to particles (PIC/FLIP combo, mostly flip tho)
    real flipRatio = 0.95f;
    gridToParticles(particles, flipRatio);

    // Advect particles
    advectParticles(particles, dt);
}

