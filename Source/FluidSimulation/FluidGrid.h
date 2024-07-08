#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "FluidGrid.generated.h"

class UProceduralMeshComponent;

UCLASS()
class FLUIDSIMULATION_API AFluidGrid : public AActor
{
	GENERATED_BODY()

public:
	AFluidGrid();

protected:
	virtual void BeginPlay() override;

public:
	virtual void Tick(float DeltaTime) override;

private:
	// Size of the simulation grid
	UPROPERTY(EditAnywhere)
	int32 Size;

	// Number of iterations for the solver
	UPROPERTY(EditAnywhere)
	int32 Iteration = 10;

	// Time step for the simulation
	UPROPERTY(EditAnywhere)
	float Dt;

	// Diffusion rate for the fluid
	UPROPERTY(EditAnywhere)
	float Diffusion;

	// Viscosity of the fluid
	UPROPERTY(EditAnywhere)
	float Viscosity;

	// Arrays to store fluid properties
	TArray<float> S;          // Source array for density
	TArray<float> Density;    // Density array

	// Velocity components
	TArray<float> Vx;
	TArray<float> Vy;
	TArray<float> Vz;

	// Previous velocity components
	TArray<float> Vx0;
	TArray<float> Vy0;
	TArray<float> Vz0;

	// Vorticity array to store curl of the velocity field
	TArray<FVector> Vorticity;

	// Max vorticity value to control the strength of the turbulence
	float MaxVorticity = 2.0f;

	// Scale of noise added to the velocity field
	float NoiseScale = 0.1f;

	UPROPERTY(EditAnywhere)
	UProceduralMeshComponent* ProceduralMesh;

	// Function to add density at a specific grid cell
	void AddDensity(int32 x, int32 y, int32 z, float amount);

	// Function to add velocity at a specific grid cell
	void AddVelocity(int32 x, int32 y, int32 z, float amountX, float amountY, float amountZ);

	// Convert 3D coordinates to 1D index
	int32 IX(int32 x, int32 y, int32 z);

	// Diffuse function to spread out the density and velocity
	void diffuse(int32 b, TArray<float>& x, TArray<float>& x0, float diff, float dt);

	// Solve linear equations
	void LinearSolve(int32 b, TArray<float>& x, const TArray<float>& x0, float a, float c);

	// Handle boundary conditions
	void SetBoundary(int32 b, TArray<float>& x);

	// Project function to enforce incompressibility
	void project(TArray<float>& velocX, TArray<float>& velocY, TArray<float>& velocZ, TArray<float>& p, TArray<float>& div);

	// Advection function to move fluid based on velocity
	void Advect(int32 b, TArray<float>& d, const TArray<float>& d0, const TArray<float>& velocX, const TArray<float>& velocY, const TArray<float>& velocZ, float dt);

	void GenerateGrid();
	void StepSimulation();
	void UpdateMesh();

private:
	FTimerHandle TimerHandle;
};
