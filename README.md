# Fluid Simulation Project

## Overview

This project implements a 2D fluid simulation using Unreal Engine. The simulation is based on the Navier-Stokes equations and provides a visual representation of fluid dynamics.

## Classes

### AFluidGrid

The `AFluidGrid` class is the core of the fluid simulation. It initializes properties such as density, velocity, and diffusion, and updates the simulation at each tick. The visual representation of the simulation is rendered onto a `UTextureRenderTarget2D` and displayed on a plane mesh.

#### Key Properties
- **Size**: The resolution of the simulation grid.
- **AreaSize**: The size of the area affected by the simulation.
- **AffectedDensity**: The density added to the grid when affected.
- **AffectedVelocity**: The velocity added to the grid when affected.
- **Dt**: The time step for the simulation.
- **Diffusion**: The diffusion rate of the fluid.
- **Viscosity**: The viscosity of the fluid.
- **TurbulenceScale**: The scale of the turbulence effect.
- **TurbulenceSpeed**: The speed of the turbulence effect.

#### Key Methods
- `InitializeRenderTarget()`: Initializes the render target for the simulation.
- `BeginPlay()`: Called when the game starts or when the actor is spawned. Initializes the render target and material instance.
- `Tick(float DeltaSeconds)`: Called every frame to update the simulation. Handles input and updates the fluid properties.
- `HandleInput()`: Handles user input to manipulate the simulation.
- `RenderDensity()`: Renders the density field onto the render target.
- `RenderVelocity()`: Renders the velocity field.
- `FadeDensity()`: Gradually fades the density field over time.
- `LineTraceAndColor()`: Performs a line trace to detect mouse clicks and updates the simulation accordingly.
- `UpdateRenderTarget()`: Updates the render target with the current simulation data.
- `GetSmoothGradientColor(float Intensity)`: Returns a color based on the intensity of the fluid properties.
- `AddDensity(int32 x, int32 y, float amount)`: Adds density to a specific grid cell.
- `AddVelocity(int32 x, int32 y, float amountX, float amountY)`: Adds velocity to a specific grid cell.
- `StepSimulation()`: Performs a single step of the fluid simulation, updating density and velocity fields.
- `AddRandomCentralVelocity(float magnitude)`: Adds a random velocity to the center of the grid.
- `Diffuse(int32 b, TArray<float>& x, TArray<float>& x0, float diff, float dt)`: Diffuses the fluid properties.
- `Advect(int32 b, TArray<float>& d, TArray<float>& d0, TArray<float>& velocX, TArray<float>& velocY, float dt)`: Advects the fluid properties based on velocity.
- `Project(TArray<float>& velocX, TArray<float>& velocY, TArray<float>& p, TArray<float>& div)`: Projects the velocity field to ensure incompressibility.
- `LinearSolve(int32 b, TArray<float>& x, TArray<float>& x0, float a, float c)`: Solves linear systems for diffusion and projection steps.
- `SetBoundary(int32 b, TArray<float>& x)`: Sets the boundary conditions for the fluid properties.
- `IX(int32 x, int32 y) const`: Converts 2D grid coordinates to a 1D array index.

## Features

- **Real-time Fluid Simulation**: Updates and renders the fluid simulation in real-time.
- **User Interaction**: Allows users to interact with the simulation through mouse clicks.
- **Adjustable Parameters**: Parameters like density, velocity, diffusion, and viscosity can be adjusted to see different fluid behaviors.
- **Visual Representation**: Uses a gradient color scheme to visualize the density and velocity of the fluid.

## Development Process

1. **Initialize Properties**: Set up the grid resolution, diffusion rates, and other fluid properties.
2. **Initialize Render Target**: Create and configure the render target for visualizing the simulation.
3. **Handle User Input**: Implement functions to detect and respond to user interactions.
4. **Simulation Steps**: Develop the core simulation loop to update fluid properties using advection, diffusion, and projection methods.
5. **Render Simulation**: Implement rendering functions to visualize the fluid properties on a 2D plane.

## Credits

This project is inspired by and based on the fluid simulation techniques discussed in the following resources:
- [Fluid Simulation for Dummies](https://mikeash.com/pyblog/fluid-simulation-for-dummies.html) by Mike Ash
- [Real-Time Fluid Dynamics for Games](https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/GDC03.pdf) by Jos Stam

Special thanks to these authors for their detailed explanations and insights into fluid dynamics simulation.
[YouTube Video](https://youtu.be/W2d810dHGM8?si=tdzvJ25-kiccODCy)

![GIF](https://media.giphy.com/media/uuVpe02tsJV21hmTqF/giphy.gif)

![image](https://github.com/user-attachments/assets/6a1a610a-328f-44d3-922c-1089eb5c82b2)
![image](https://github.com/user-attachments/assets/f41c514c-f1d5-457a-a170-bfedd688186b)
![image](https://github.com/user-attachments/assets/b55749ad-6201-4c31-a8ad-fd1e4d586f19)
![image](https://github.com/user-attachments/assets/5517c1c6-7792-4335-b634-e3dd731d6999)
![image](https://github.com/user-attachments/assets/02b7ec54-5ba3-40cf-b03f-c3722510b1bc)
![image](https://github.com/user-attachments/assets/801d6bff-7d93-42a9-9077-81400f88625d)



