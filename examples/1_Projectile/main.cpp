/*
 *  Assignment 1: Projectile
 *  GPR-350-01 Fall '25
 *
 *  Author: Arija Hartel (@cocoatehcat)
 *  Due: 09/08/2025
 *
 *  Certificate of Authenticity: I did not use AI or any help except for the book and base Google,
 *  and enlisting Victor Diab to look at my code when my particle trail was being grumpy.
 */

#include <iostream>
#include <vector>

#include <raylib.h>
#include <raymath.h> //Needed for most Raylib math functions
#include <cyclone/cyclone.h>
#include <imgui.h>
#include <rlImGui.h>
#include <bits/stl_algo.h>

#include <cyclone/particle.h>
#include <ew/camera.h>

struct AppSettings {
    float spinSpeed = 2.0f; // Radians per second
    Vector3 position; // Position of the cube

    float angle;
    float speed;
    Vector3 acceleration;
    float damping;
    Vector2 initialVelocity;
}_appSettings;

int main() {
    // Initialization
   //--------------------------------------------------------------------------------------
    SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE);
    InitWindow(1280, 720, "Projectile Demo");

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

    // setting values
    Vector3 sphereAcceleration = { 0.0f, -9.8f, 0.0f };
    float sphereSpeed = 10.0f;
    float sphereDamping = 0.9f;

    bool sphereMove = false;

    // list for particle trail
    std::vector<cyclone::Particle> particleTrail;
    auto mainProj = new cyclone::Particle();

    // Set default, change to 0
    spherePosition.x = 0;
    spherePosition.y = 3;
    spherePosition.z = 0;

    _appSettings.position.y = spherePosition.y;

    _appSettings.acceleration = sphereAcceleration;
    _appSettings.speed = sphereSpeed;
    _appSettings.damping = sphereDamping;

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

        // Move?
        if (sphereMove) {
            // integrate?
            mainProj->integrate((float)deltaTime);
            spherePosition = mainProj->getPosition();

            // Check to see if sphere should stop moving
            if (spherePosition.y < 0.0f) {
            sphereMove = false;
            }
        }

        //Model.transform allows us to modify the model (local->world) matrix directly
        //Concatenation order is left to right
        //Note that Raylib does come with optional operator overloads for C++!
        sphereModel.transform = MatrixTranslate(spherePosition.x, spherePosition.y, spherePosition.z);

        //Drawing
        BeginDrawing();

        //3D mode draws objects in right handed vector space
        BeginMode3D(camera);
        ClearBackground(SKYBLUE);

        //Position is zeroed because Model.transform is handling the position for us
        DrawModel(sphereModel, { 0,0,0 }, 1.0f, RED);
        DrawModelWires(sphereModel, { 0,0,0 }, 1.0f, BLACK);

        // for loop of all the baby spheres in collection
        if (sphereMove) {
            for (auto i = 0; i < particleTrail.size(); i++) {
                sphereModel.transform = MatrixTranslate(particleTrail[i].getPosition().x, particleTrail[i].getPosition().y, particleTrail[i].getPosition().z);
                DrawModel(sphereModel, { 0,0,0 }, 0.25f, GRAY);
            }
        }

        DrawGrid(100, 1.0f);

        EndMode3D();

        //Outside of 3D mode, we can draw 2D things directly to the screen in pixel coordinates
        //DrawCircle(32, 32, 16.f, BLUE);
        //DrawText("RayLib Text!", 64, 16, 32, WHITE);
        DrawFPS(GetScreenWidth() - 128, 16);

        //ImGUI
        //All ImGui calls must be between rlGuiBegin() and rlImGuiEnd()
        // ImGUI returns bool checking if it has changed, can collect those and refresh at the end to be cool
        // DragFloat3, SliderFloat, Button for Start
        // CollapsingHeader to nest things
        rlImGuiBegin();
        ImGui::Begin("SETTINGS");
        // Initial Position:
        ImGui::InputFloat3( "Position", (float*)&_appSettings.position);
        // Angle
        ImGui::SliderFloat("Angle", &_appSettings.angle, 0.0f, 90.0f);
        // Speed
        ImGui::SliderFloat("Speed", &_appSettings.speed, 0.0f, 60.0f);
        // Acceleration
        ImGui::InputFloat3("Acceleration", (float*)&_appSettings.acceleration);
        // Damping
        ImGui::SliderFloat("Damping", &_appSettings.damping, 0.0f, 1.0f);
        if (ImGui::Button("Start!", ImVec2(100, 20))) {
            sphereMove = true; // create projectile, add coalition to queue
            // Change to radians
            float angleInRadians = _appSettings.angle * (PI / 180.0f);

            // Initial Stuff
            cyclone::real initialVelocityX = _appSettings.speed * cos(angleInRadians);
            cyclone::real initialVelocityY = _appSettings.speed * sin(angleInRadians);
            _appSettings.initialVelocity.x = initialVelocityX;
            _appSettings.initialVelocity.y = initialVelocityY;

            // Main guy
            auto projectile = new cyclone::Particle();

            // Position
            projectile->setPosition(_appSettings.position.x, _appSettings.position.y, _appSettings.position.z);

            // Velocity
            projectile->setVelocity(initialVelocityX, initialVelocityY, 0);

            // Acceleration (currently gravity)
            projectile->setAcceleration(_appSettings.acceleration.x, _appSettings.acceleration.y, _appSettings.acceleration.z);

            // Damping
            projectile->setDamping(_appSettings.damping);

            // The whole point
            mainProj = projectile;

            // Trail
            float fixedTime = 1.0f / 60.0f;

            // Clear
            particleTrail.clear();

            // Multiplying by 4 since the scale is lowered to make ball smaller
            Vector3 fakePosition;
            fakePosition.x = _appSettings.position.x * 4.0f;
            fakePosition.y = _appSettings.position.y * 4.0f;
            fakePosition.z = _appSettings.position.z * 4.0f;

            Vector3 fakeVelocity;
            fakeVelocity.x = initialVelocityX * 4.0f;
            fakeVelocity.y = initialVelocityY * 4.0f;
            fakeVelocity.z = 0 * 4.0f;

            // Calculating how many balls are needed dynamically
            // Got this from Google & help from Victor
            float landTime = (-fakeVelocity.y - sqrt(fakeVelocity.y * fakeVelocity.y - (2 * _appSettings.acceleration.y * fakePosition.y))) / _appSettings.acceleration.y;

            for (float t = 0; t <= landTime; t += deltaTime) {

                // Setting the bases
                auto fakeProjectile = new cyclone::Particle();
                fakeProjectile->setPosition(fakePosition.x, fakePosition.y, fakePosition.z);
                fakeProjectile->setVelocity(fakeVelocity.x, fakeVelocity.y, fakeVelocity.z);
                fakeProjectile->setAcceleration(_appSettings.acceleration.x * 4.0f, _appSettings.acceleration.y * 4.0f, _appSettings.acceleration.z * 4.0f);
                fakeProjectile->setDamping(_appSettings.damping);

                fakeProjectile->integrate((float)fixedTime);

                // Setting variables for next particle
                fakePosition.x = fakeProjectile->getPosition().x;
                fakePosition.y = fakeProjectile->getPosition().y;
                fakePosition.z = fakeProjectile->getPosition().z;

                fakeVelocity.x = fakeProjectile->getVelocity().x;
                fakeVelocity.y = fakeProjectile->getVelocity().y;
                fakeVelocity.z = fakeProjectile->getVelocity().z;

                if (fakeProjectile->getPosition().y < 0.0f) {
                    continue;
                }

                particleTrail.push_back(*fakeProjectile);
            }
        }

        ImGui::End();

        rlImGuiEnd(); // end of on window

        EndDrawing();
    }
    UnloadModel(sphereModel);
    CloseWindow();
}