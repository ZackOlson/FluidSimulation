/*
 *  Final: Fluid Simulation
 *  GPR-350-01 Fall '25
 *
 *  Author: Arija Hartel (@cocoatehcat) & Zachary Olson (@ZackOlson)
 *  Due: 12/08/2025
 *
 *  Certificate of Authenticity: We did not use AI and all of our resources are linked in the
 *  "Works Used" ReadMe
 */

#include <iostream>
#include <vector>

#include <raylib.h>
#include <raymath.h> //Needed for most Raylib math functions
#include <cyclone/cyclone.h>
#include <imgui.h>
#include <rlImGui.h>
//#include <bits/stl_algo.h>

#include <cyclone/particle.h>
#include <ew/camera.h>
#include <cyclone/fluidsim.h>

struct AppSettings {
    int numParticles;
    float spawnRad;
    int gridWidth;
    int gridHeight;
    int gridDepth;
    float cellSize;
    float density;

}_appSettings;

int main() {
    // Initialization
   //--------------------------------------------------------------------------------------
    SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE);
    InitWindow(1280, 720, "Final Demo");

    rlImGuiSetup(true);

    Camera3D camera = { 0 };
    camera.position = {0.0f, 10.0f, 15.0f };   // Camera position
    camera.target = { 0.0f, 0.0f, 0.0f };       // Camera looking at point
    camera.up = { 0.0f, 1.0f, 0.0f };           // Camera up vector
    camera.fovy = 45.0f;                        // Camera vertical field-of-view in degrees
    camera.projection = CAMERA_PERSPECTIVE;     // Camera projection type - perspective vs orthographic

    //Vector used by the Cyclone physics engine
    cyclone::Vector3 spherePosition(0.0f, 0.0f, 0.0f);
    Model sphereModel = LoadModelFromMesh(GenMeshSphere(1, 12, 6)); // looks like a beach ball so I'm keeping it

    // Billboard texture for particles
    Image circleImg = GenImageColor(32, 32, BLANK);
    ImageDrawCircle(&circleImg, 16, 16, 16, WHITE);
    Texture2D circleTexture = LoadTextureFromImage(circleImg);
    UnloadImage(circleImg);

    bool sphereMove = false;
    bool applied = false;

    // FLUID
    // Number of particles and where (random numbers at the moment)
    const int NUM_PARTICLES = 5000;
    const float SPAWN_RADIUS = 8.0f;

    std::vector<cyclone::Particle> fluidParticles;

    fluidParticles.resize(NUM_PARTICLES);

    // Spawn particles in a cube
    for (int i = 0; i < NUM_PARTICLES; i++) {
        float x = ((float)rand() / RAND_MAX - 0.5f) * SPAWN_RADIUS;
        float y = ((float)rand() / RAND_MAX) * SPAWN_RADIUS + 1.0f;
        float z = ((float)rand() / RAND_MAX - 0.5f) * SPAWN_RADIUS;

        // Correct position
        x += 5.0f;
        y += 5.0f;
        z += 5.0f;

        fluidParticles[i].setPosition(x, y, z);
        fluidParticles[i].setVelocity(0, 0, 0);
        fluidParticles[i].setAcceleration(0, -9.8f, 0);
        fluidParticles[i].setDamping(0.98f);
        fluidParticles[i].setMass(1.0f);
    }

    int gridW = 30;
    int gridH = 30;
    int gridD = 30;
    float cellSize = 0.7f;
    float density = 1000.0f;

    _appSettings.numParticles = NUM_PARTICLES;
    _appSettings.spawnRad = SPAWN_RADIUS;
    _appSettings.gridWidth = gridW;
    _appSettings.gridHeight = gridH;
    _appSettings.gridDepth = gridD;
    _appSettings.cellSize = cellSize;
    _appSettings.density = density;

    cyclone::FluidSim fluidSim(gridW, gridH, gridD, cellSize, density);


    while (!WindowShouldClose()) {
        //Seconds between previous frame and this one
        float deltaTime = GetFrameTime();
        
        //Input
        if (IsMouseButtonPressed(MOUSE_RIGHT_BUTTON)) {
            DisableCursor();
        }
        if (IsMouseButtonReleased(MOUSE_RIGHT_BUTTON)) {
            EnableCursor();
        }
        //Only allow movement if the cursor is hidden
        if (IsCursorHidden()) {
            ew::UpdateFlyCamera(&camera, deltaTime);
        }

        // Running?
        if (sphereMove) {
            fluidSim.step(fluidParticles, deltaTime);
        }


        //Drawing
        BeginDrawing();

        //3D mode draws objects in right handed vector space
        BeginMode3D(camera);
        ClearBackground(SKYBLUE);

        //Position is zeroed because Model.transform is handling the position for us
        for (auto& p : fluidParticles) {
            Vector3 pos = { p.getPosition().x, p.getPosition().y, p.getPosition().z };

            float size = 0.2f; // diameter of the particle
            DrawBillboard(camera, circleTexture, pos, size, BLUE);
        }

        DrawGrid(100, 1.0f);

        EndMode3D();


        DrawText("Fluid Sim!", 64, 16, 32, WHITE);
        DrawFPS(GetScreenWidth() - 128, 16);

        //ImGUI
        //All ImGui calls must be between rlGuiBegin() and rlImGuiEnd()
        // ImGUI returns bool checking if it has changed, can collect those and refresh at the end to be cool
        // DragFloat3, SliderFloat, Button for Start
        // CollapsingHeader to nest things
        rlImGuiBegin();
        ImGui::Begin("SETTINGS");
        // Particles
        ImGui::SliderInt("Number of Particles", &_appSettings.numParticles, 10, 5000);
        // Spawn Radius
        ImGui::SliderFloat("Spawn Radius", &_appSettings.spawnRad, 1.0f, 10.0f);
        // Width
        ImGui::SliderInt("Grid Width", &_appSettings.gridWidth, 10, 50);
        // Height
        ImGui::SliderInt("Grid Height", &_appSettings.gridHeight, 10, 50);
        // Depth
        ImGui::SliderInt("Grid Depth", &_appSettings.gridDepth, 10, 50);
        // Cell Size
        ImGui::SliderFloat("Cell Size", &_appSettings.cellSize, 0.5f, 5.0f);
        // Density
        ImGui::SliderFloat("Density", &_appSettings.density, 10.0f, 5000.0f);
        if (ImGui::Button("Apply", ImVec2(100, 20))) {
            // recreate particle vector
            fluidParticles.clear();
            fluidParticles.resize(_appSettings.numParticles);

            // spawn particles
            for (int i = 0; i < _appSettings.numParticles; ++i) {
                float x = ((float)rand() / RAND_MAX - 0.5f) * _appSettings.spawnRad;
                float y = ((float)rand() / RAND_MAX) * _appSettings.spawnRad * 0.5f + 1.0f; // lower spawn Y
                float z = ((float)rand() / RAND_MAX - 0.5f) * _appSettings.spawnRad;

                // place
                fluidParticles[i].setPosition(x, y, z);
                fluidParticles[i].setVelocity(0, 0, 0);
                fluidParticles[i].setAcceleration(0, -9.8f, 0);
                fluidParticles[i].setDamping(0.98f);
                fluidParticles[i].setMass(1.0f);
            }

            // Recreate fluid sim
            fluidSim = cyclone::FluidSim(_appSettings.gridWidth,
                _appSettings.gridHeight,
                _appSettings.gridDepth,
                _appSettings.cellSize,
                _appSettings.density);

            applied = false;
        }
        if (ImGui::Button("Start!", ImVec2(100, 20))) {
            sphereMove = true; // start sim
        }
        // Stop
        if (ImGui::Button("Stop!", ImVec2(100, 20))) {
            sphereMove = false; // pause
        }

        ImGui::End();

        rlImGuiEnd(); // end of on window

        EndDrawing();
    }
    UnloadTexture(circleTexture);
    UnloadModel(sphereModel);
    CloseWindow();
}