#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Engine/TextureRenderTarget2D.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "FluidGrid.generated.h"

UCLASS()
class FLUIDSIMULATION_API AFluidGrid : public AActor
{
	GENERATED_BODY()

public:
	AFluidGrid();

protected:
	virtual void BeginPlay() override;
	virtual void Tick(float DeltaSeconds) override;

	void HandleInput();

private:
	UPROPERTY(EditAnywhere)
	int32 Size = 50;

	UPROPERTY(EditAnywhere)
	float Dt = 0.1f;

	UPROPERTY(EditAnywhere)
	float Diffusion = 1.0f;

	UPROPERTY(EditAnywhere)
	float Viscosity = 0.1f;

	TArray<float> Density, Vx, Vy, Vz;

	UPROPERTY(VisibleAnywhere)
	UTextureRenderTarget2D* RenderTarget;

	UPROPERTY(VisibleAnywhere)
	UMaterialInstanceDynamic* DynamicMaterialInstance;

	UPROPERTY(VisibleAnywhere)
	UStaticMeshComponent* PlaneComponent;

	UPROPERTY(EditAnywhere, Category = "Fluid Simulation")
	UMaterial* BaseMaterial;

	void InitializeRenderTarget();
	void UpdateRenderTarget();
	void AddDensity(int32 x, int32 y, float amount);
	void AddVelocity(int32 x, int32 y, float amountX, float amountY);
	void StepSimulation();
	void Diffuse(int32 b, TArray<float>& x, TArray<float>& x0, float diff, float dt);
	void Advect(int32 b, TArray<float>& d, TArray<float>& d0, TArray<float>& velocX, TArray<float>& velocY, float dt);
	void Project(TArray<float>& velocX, TArray<float>& velocY, TArray<float>& p, TArray<float>& div);
	void LinearSolve(int32 b, TArray<float>& x, TArray<float>& x0, float a, float c);
	void SetBoundary(int32 b, TArray<float>& x);
	int32 IX(int32 x, int32 y) const { return x + y * Size; }
};
