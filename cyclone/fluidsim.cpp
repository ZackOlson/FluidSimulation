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

    const real h = cellSize;
    const real invh = invCellSize;
    const real invh2 = invh * invh;


    std::vector<real> div(N, 0.0f);

    for (int z = 1; z < gridDepth - 1; ++z)
        for (int y = 1; y < gridHeight - 1; ++y)
            for (int x = 1; x < gridWidth - 1; ++x)
            {
                int k = idx3(x, y, z);

                real du_dx = (u[idx3(x + 1, y, z)] - u[idx3(x - 1, y, z)]) * 0.5f * invh;
                real dv_dy = (v[idx3(x, y + 1, z)] - v[idx3(x, y - 1, z)]) * 0.5f * invh;
                real dw_dz = (w[idx3(x, y, z + 1)] - w[idx3(x, y, z - 1)]) * 0.5f * invh;

                div[k] = (du_dx + dv_dy + dw_dz) / dt;
            }


    std::fill(pressure.begin(), pressure.end(), 0.0f);


    std::vector<real> r = div;
    std::vector<real> z(N, 0.0f);
    std::vector<real> p(N, 0.0f);
    std::vector<real> Ap(N, 0.0f);

    const real M = (h * h) / 6.0f;
    for (int i = 0; i < N; ++i) z[i] = r[i] * M;
    p = z;

    real rzOld = 0.0f;
    for (int i = 0; i < N; ++i) rzOld += r[i] * z[i];


    const int maxIt = 60;
    for (int it = 0; it < maxIt; ++it)
    {

        for (int zc = 1; zc < gridDepth - 1; ++zc)
            for (int yc = 1; yc < gridHeight - 1; ++yc)
                for (int xc = 1; xc < gridWidth - 1; ++xc)
                {
                    int k = idx3(xc, yc, zc);

                    auto sample = [&](int X, int Y, int Z)
                        {
                            // If out of bounds, return center
                            if (X < 0 || X >= gridWidth ||
                                Y < 0 || Y >= gridHeight ||
                                Z < 0 || Z >= gridDepth)
                                return p[k];
                            return p[idx3(X, Y, Z)];
                        };

                    real pl = sample(xc - 1, yc, zc);
                    real pr = sample(xc + 1, yc, zc);
                    real pb = sample(xc, yc - 1, zc);
                    real pt = sample(xc, yc + 1, zc);
                    real pd = sample(xc, yc, zc - 1);
                    real pu = sample(xc, yc, zc + 1);

                    Ap[k] = invh2 * (pl + pr + pb + pt + pd + pu - 6.0f * p[k]);
                }

        // Get alpha
        real denom = 0;
        for (int i = 0; i < N; ++i) denom += p[i] * Ap[i];
        if (fabs(denom) < 1e-20f) break;

        real alpha = rzOld / denom;

        for (int i = 0; i < N; ++i) pressure[i] += alpha * p[i];

        for (int i = 0; i < N; ++i) r[i] -= alpha * Ap[i];

        // Check convergence
        real rr = 0;
        for (int i = 0; i < N; i++) rr += r[i] * r[i];
        if (rr < 1e-10f) break;

        for (int i = 0; i < N; ++i) z[i] = r[i] * M;

        real rzNew = 0.0f;
        for (int i = 0; i < N; ++i) rzNew += r[i] * z[i];

        real beta = rzNew / rzOld;
        rzOld = rzNew;

        for (int i = 0; i < N; ++i) p[i] = z[i] + beta * p[i];
    }

    // X walls
    for (int y = 0; y < gridHeight; ++y)
        for (int z = 0; z < gridDepth; ++z)
        {
            pressure[idx3(0, y, z)] = pressure[idx3(1, y, z)];
            pressure[idx3(gridWidth - 1, y, z)] = pressure[idx3(gridWidth - 2, y, z)];
        }
    // Y walls
    for (int x = 0; x < gridWidth; ++x)
        for (int z = 0; z < gridDepth; ++z)
        {
            pressure[idx3(x, 0, z)] = pressure[idx3(x, 1, z)];
            pressure[idx3(x, gridHeight - 1, z)] = pressure[idx3(x, gridHeight - 2, z)];
        }
    // Z walls
    for (int x = 0; x < gridWidth; ++x)
        for (int y = 0; y < gridHeight; ++y)
        {
            pressure[idx3(x, y, 0)] = pressure[idx3(x, y, 1)];
            pressure[idx3(x, y, gridDepth - 1)] = pressure[idx3(x, y, gridDepth - 2)];
        }

    for (int zc = 1; zc < gridDepth - 1; ++zc)
        for (int yc = 1; yc < gridHeight - 1; ++yc)
            for (int xc = 1; xc < gridWidth - 1; ++xc)
            {
                int k = idx3(xc, yc, zc);
                u[k] -= (pressure[idx3(xc + 1, yc, zc)] - pressure[idx3(xc - 1, yc, zc)]) * 0.5f * invh;
                v[k] -= (pressure[idx3(xc, yc + 1, zc)] - pressure[idx3(xc, yc - 1, zc)]) * 0.5f * invh;
                w[k] -= (pressure[idx3(xc, yc, zc + 1)] - pressure[idx3(xc, yc, zc - 1)]) * 0.5f * invh;
            }

    for (int z = 0; z < gridDepth; ++z)
        for (int x = 0; x < gridWidth; ++x) {
            v[idx3(x, 0, z)] = v[idx3(x, gridHeight - 1, z)] = 0;
            u[idx3(x, 0, z)] = u[idx3(x, gridHeight - 1, z)] = 0;
            w[idx3(x, 0, z)] = w[idx3(x, gridHeight - 1, z)] = 0;
        }
    for (int z = 0; z < gridDepth; ++z)
        for (int y = 0; y < gridHeight; ++y) {
            u[idx3(0, y, z)] = u[idx3(gridWidth - 1, y, z)] = 0;
            v[idx3(0, y, z)] = v[idx3(gridWidth - 1, y, z)] = 0;
            w[idx3(0, y, z)] = w[idx3(gridWidth - 1, y, z)] = 0;
        }
    for (int y = 0; y < gridHeight; ++y)
        for (int x = 0; x < gridWidth; ++x) {
            w[idx3(x, y, 0)] = w[idx3(x, y, gridDepth - 1)] = 0;
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
    // Build cell lists
    int nCells = gridWidth * gridHeight * gridDepth;

    std::vector<std::vector<int>> cellParticles;
    cellParticles.resize(nCells);

    auto cellIndexFromCoords = [&](int x, int y, int z)->int {
        return x + y * gridWidth + z * gridWidth * gridHeight;
        };

    // clamp cell coords
    auto clampCell = [&](int& c, int maxC) {
        if (c < 0) c = 0;
        else if (c >= maxC) c = maxC - 1;
        };

    // Insert particles into cells
    const int P = (int)particles.size();
    for (int i = 0; i < P; ++i) {
        Vector3 pos = particles[i].getPosition();

        int cx = (int)floor(pos.x * invCellSize);
        int cy = (int)floor(pos.y * invCellSize);
        int cz = (int)floor(pos.z * invCellSize);

        clampCell(cx, gridWidth);
        clampCell(cy, gridHeight);
        clampCell(cz, gridDepth);

        int idx = cellIndexFromCoords(cx, cy, cz);
        cellParticles[idx].push_back(i);
    }

    // For each particle, test only neighboring cells
    const real restSize = 1.0f;
    const real stiffness = 1.0f;
    const real invDt = (dt > 0.0f) ? 1.0f / dt : 0.0f;

    for (int i = 0; i < P; ++i) {
        Vector3 pos_i = particles[i].getPosition();

        int cx = (int)floor(pos_i.x * invCellSize);
        int cy = (int)floor(pos_i.y * invCellSize);
        int cz = (int)floor(pos_i.z * invCellSize);

        // neighbor cell bounds (clamped)
        int minx = cx - 1, miny = cy - 1, minz = cz - 1;
        int maxx = cx + 1, maxy = cy + 1, maxz = cz + 1;
        if (minx < 0) minx = 0;
        if (miny < 0) miny = 0;
        if (minz < 0) minz = 0;
        if (maxx >= gridWidth) maxx = gridWidth - 1;
        if (maxy >= gridHeight) maxy = gridHeight - 1;
        if (maxz >= gridDepth) maxz = gridDepth - 1;

        Vector3 tempForce(0, 0, 0);

        for (int z = minz; z <= maxz; ++z) {
            for (int y = miny; y <= maxy; ++y) {
                for (int x = minx; x <= maxx; ++x) {
                    int cellIdx = cellIndexFromCoords(x, y, z);
                    const auto& list = cellParticles[cellIdx];
                    for (int id : list) {
                        if (id == i) continue;
                        Vector3 pos_j = particles[id].getPosition();
                        Vector3 dir = pos_i - pos_j;
                        real dist = dir.magnitude();
                        if (dist <= 0.0001f) continue;
                        if (dist < restSize) {
                            tempForce += dir.unit() * stiffness * (restSize - dist);
                        }
                    }
                }
            }
        }


        Vector3 vel = particles[i].getVelocity();

        Vector3 accel = tempForce * (1.0f / particles[i].getMass()); // mass is 1.0f by your init, but keep general
        vel += accel * dt;

        vel.y += gravity * dt;

        // Update position from velocity
        Vector3 newPos = particles[i].getPosition() + vel * dt;

        particles[i].setVelocity(vel);
        particles[i].setPosition(newPos);

        // handle walls
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
    real flipRatio = 0.05f;
    gridToParticles(particles, flipRatio);

    // Advect particles
    advectParticles(particles, dt);
}

