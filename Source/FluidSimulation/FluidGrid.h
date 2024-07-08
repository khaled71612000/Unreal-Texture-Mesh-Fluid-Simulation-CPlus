#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "FluidGrid.generated.h"

class UTexture2D;
class UMaterialInstanceDynamic;
class UStaticMeshComponent;

UCLASS()
class FLUIDSIMULATION_API AFluidGrid : public AActor
{
	GENERATED_BODY()

public:
	//Initializes the fluid grid's size, time step, diffusion, viscosity, and property arrays. Creates the plane component for visualization.
	AFluidGrid();

protected:

	//Schedules the simulation step to run periodically.
	virtual void BeginPlay() override;

	void UpdateTexture();
private:
	//Dimension of the fluid grid.
	UPROPERTY(EditAnywhere)
	int32 Size = 100;
	//Time step for the simulation.
	UPROPERTY(EditAnywhere)
	float Dt = 0.1f;
	//Diffusion rate of the fluid.
	UPROPERTY(EditAnywhere)
	float Diffusion = 0.1f;
	//Viscosity of the fluid.
	UPROPERTY(EditAnywhere)
	float Viscosity = 0.1f;
	//Arrays to store fluid properties like density and velocity in x, y, z directions.
	TArray<float> Density, Vx, Vy, Vz;
	// Texture to visualize the fluid simulation.
	UPROPERTY(VisibleAnywhere)
	UTexture2D* DynamicTexture;
	//Material instance to apply the dynamic texture.
	UPROPERTY(VisibleAnywhere)
	UMaterialInstanceDynamic* DynamicMaterialInstance;
	// Static mesh component to display the fluid grid.
	UPROPERTY(VisibleAnywhere)
	UStaticMeshComponent* PlaneComponent;
	//Material to base the dynamic material on.
	UPROPERTY(EditAnywhere, Category = "Fluid Simulation")
	UMaterial* BaseMaterial;

	//Creates a transient texture to visualize the fluid simulation.
	void InitializeDynamicTexture();

	//Adds a specified amount of density at a given grid position.
	void AddDensity(int32 x, int32 y, float amount);
	void AddVelocity(int32 x, int32 y, float amountX, float amountY);

	/*Main simulation step that updates the fluid properties :
	Adds density and velocity for testing.
		Diffuses velocities.
		Projects velocities to ensure incompressibility.
		Advects velocities and densities.*/

	void StepSimulation();
	//Diffuses the fluid properties using linear solve.
	void Diffuse(int32 b, TArray<float>& x, TArray<float>& x0, float diff, float dt);
	//Moves fluid properties through the grid based on velocity.
	void Advect(int32 b, TArray<float>& d, TArray<float>& d0, TArray<float>& velocX, TArray<float>& velocY, float dt);
	//Ensures fluid velocity field is incompressible.
	void Project(TArray<float>& velocX, TArray<float>& velocY, TArray<float>& p, TArray<float>& div);
	void LinearSolve(int32 b, TArray<float>& x, TArray<float>& x0, float a, float c);
	void SetBoundary(int32 b, TArray<float>& x);
	//Converts 2D grid coordinates to a 1D index.
	int32 IX(int32 x, int32 y) const { return x + y * Size; }
	FTimerHandle TimerHandle;
};