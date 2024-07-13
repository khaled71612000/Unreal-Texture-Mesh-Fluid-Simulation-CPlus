#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "Engine/TextureRenderTarget2D.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "Components/BoxComponent.h"
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
	void LineTraceAndColor();

private:
	UPROPERTY(EditAnywhere)
	int32 Size = 128;

	UPROPERTY(EditAnywhere, Category = "Fluid Simulation")
	int32 AreaSize = 60;  // Adjust this value to make the affected area bigger

	UPROPERTY(EditAnywhere, Category = "Fluid Simulation")
	float AffectedDensity = 10;  // Amount of density to add

	UPROPERTY(EditAnywhere, Category = "Fluid Simulation")
	float AffectedVelocity = 2;  // Amount of velocity to add

	UPROPERTY(EditAnywhere)
	float Dt = 100.f;

	UPROPERTY(EditAnywhere)
	float Diffusion = 100.0f;

	UPROPERTY(EditAnywhere)
	float Viscosity = 10.f;

	TArray<float> Density, Vx, Vy, Vz;

	UPROPERTY(VisibleAnywhere)
	UTextureRenderTarget2D* RenderTarget;

	UPROPERTY(VisibleAnywhere)
	UMaterialInstanceDynamic* DynamicMaterialInstance;

	UPROPERTY(VisibleAnywhere)
	UStaticMeshComponent* PlaneComponent;

	UPROPERTY(VisibleAnywhere)
	UBoxComponent* CollisionBox;  // New collision box component

	UPROPERTY(EditAnywhere, Category = "Fluid Simulation")
	UMaterial* BaseMaterial;

	UPROPERTY(EditAnywhere)
	float Scale = 100.0f;  // Adjust this scale factor as needed

	void InitializeRenderTarget();
	void UpdateRenderTarget();

	FColor GetGradientColor(float Intensity);
	void AddDensity(int32 x, int32 y, float amount);
	void AddVelocity(int32 x, int32 y, float amountX, float amountY);
	void Diffuse(int32 b, TArray<float>& x, TArray<float>& x0, float diff, float dt);
	void Advect(int32 b, TArray<float>& d, TArray<float>& d0, TArray<float>& velocX, TArray<float>& velocY, float dt);
	void Project(TArray<float>& velocX, TArray<float>& velocY, TArray<float>& p, TArray<float>& div);
	void LinearSolve(int32 b, TArray<float>& x, TArray<float>& x0, float a, float c);
	void SetBoundary(int32 b, TArray<float>& x);

	int32 IX(int32 x, int32 y) const;
	void StepSimulation();
	void AddRandomCentralVelocity(float magnitude);

	void RenderDensity();
	void RenderVelocity();
	void FadeDensity();

};
