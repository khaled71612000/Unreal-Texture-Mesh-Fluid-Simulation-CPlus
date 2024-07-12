#include "FluidGrid.h"
#include "Engine/World.h"
#include "TimerManager.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "UObject/ConstructorHelpers.h"
#include "Components/StaticMeshComponent.h"
#include "Engine/StaticMesh.h"
#include "RHI.h"
#include "RenderGraphResources.h"
#include "RenderCommandFence.h"
#include "RenderingThread.h"
#include "DrawDebugHelpers.h"
#include "Components/BoxComponent.h"
#include "Engine/StaticMesh.h"

AFluidGrid::AFluidGrid()
{
	PrimaryActorTick.bCanEverTick = true;

	int32 TotalSize = Size * Size;
	Density.SetNumZeroed(TotalSize);
	Vx.SetNumZeroed(TotalSize);
	Vy.SetNumZeroed(TotalSize);
	Vz.SetNumZeroed(TotalSize);

	PlaneComponent = CreateDefaultSubobject<UStaticMeshComponent>(TEXT("PlaneComponent"));
	RootComponent = PlaneComponent;

	static ConstructorHelpers::FObjectFinder<UStaticMesh> PlaneMeshAsset(TEXT("/Game/StarterContent/Shapes/Shape_Plane.Shape_Plane"));
	if (PlaneMeshAsset.Succeeded())
	{
		PlaneComponent->SetStaticMesh(PlaneMeshAsset.Object);
		PlaneComponent->SetWorldScale3D(FVector(10.0f, 10.0f, 1.0f));
		PlaneComponent->SetCollisionEnabled(ECollisionEnabled::QueryAndPhysics);
		PlaneComponent->SetCollisionResponseToAllChannels(ECR_Block);
	}

	RenderTarget = CreateDefaultSubobject<UTextureRenderTarget2D>(TEXT("RenderTarget"));
}

void AFluidGrid::InitializeRenderTarget()
{
	RenderTarget->InitAutoFormat(Size, Size);
	RenderTarget->RenderTargetFormat = ETextureRenderTargetFormat::RTF_RGBA8;
	RenderTarget->ClearColor = FLinearColor::Black;
	RenderTarget->UpdateResource();
}

void AFluidGrid::BeginPlay()
{
	Super::BeginPlay();

	InitializeRenderTarget();

	if (BaseMaterial)
	{
		DynamicMaterialInstance = UMaterialInstanceDynamic::Create(BaseMaterial, this);
		if (DynamicMaterialInstance)
		{
			DynamicMaterialInstance->SetTextureParameterValue(FName("DynamicTexture"), RenderTarget);
			PlaneComponent->SetMaterial(0, DynamicMaterialInstance);
		}
	}
}

void AFluidGrid::Tick(float DeltaSeconds)
{
	Super::Tick(DeltaSeconds);

	HandleInput();
	StepSimulation();
	UpdateRenderTarget();
}

void AFluidGrid::HandleInput()
{
	if (GetWorld()->GetFirstPlayerController()->IsInputKeyDown(EKeys::LeftMouseButton))
	{
		LineTraceAndColor();
	}
}
void AFluidGrid::LineTraceAndColor()
{
	FVector2D MousePosition;
	if (GetWorld()->GetFirstPlayerController()->GetMousePosition(MousePosition.X, MousePosition.Y))
	{
		FVector WorldPosition, WorldDirection;
		if (GetWorld()->GetFirstPlayerController()->DeprojectScreenPositionToWorld(MousePosition.X, MousePosition.Y, WorldPosition, WorldDirection))
		{
			FHitResult HitResult;
			FCollisionQueryParams Params;

			FVector EndPosition = WorldPosition + WorldDirection * 10000.0f;
			bool bHit = GetWorld()->LineTraceSingleByChannel(HitResult, WorldPosition, EndPosition, ECC_Visibility, Params);

			if (bHit && HitResult.Component == PlaneComponent)
			{
				FVector LocalHit = PlaneComponent->GetComponentTransform().InverseTransformPosition(HitResult.Location);
				FVector2D GridPosition;
				GridPosition.X = (LocalHit.X + (PlaneComponent->GetStaticMesh()->GetBounds().BoxExtent.X)) / (PlaneComponent->GetStaticMesh()->GetBounds().BoxExtent.X * 2.0f) * Size;
				GridPosition.Y = (LocalHit.Y + (PlaneComponent->GetStaticMesh()->GetBounds().BoxExtent.Y)) / (PlaneComponent->GetStaticMesh()->GetBounds().BoxExtent.Y * 2.0f) * Size;

				int32 ClampedGridX = FMath::Clamp(static_cast<int32>(GridPosition.X), 1, Size - 2);
				int32 ClampedGridY = FMath::Clamp(static_cast<int32>(GridPosition.Y), 1, Size - 2);

				// Adjust the area affected by the hit point
				for (int32 offsetX = 0; offsetX <= 0; ++offsetX)
				{
					for (int32 offsetY = 0; offsetY <= 0; ++offsetY)
					{
						int32 x = FMath::Clamp(ClampedGridX + offsetX, 1, Size - 2);
						int32 y = FMath::Clamp(ClampedGridY + offsetY, 1, Size - 2);
						AddDensity(x, y, 10.0f);
						AddVelocity(x, y, 1.0f, 1.0f); // Example velocities
					}
				}

				StepSimulation();
				UpdateRenderTarget();
			}
		}
	}
}

void AFluidGrid::UpdateRenderTarget()
{
	FTextureRenderTargetResource* RenderTargetResource = RenderTarget->GameThread_GetRenderTargetResource();
	TArray<FColor> ColorData;
	ColorData.SetNum(Size * Size);

	for (int32 y = 0; y < Size; y++)
	{
		for (int32 x = 0; x < Size; x++)
		{
			float Value = Density[IX(x, y)];
			float Intensity = FMath::Clamp(Value, 0.0f, 1.0f);

			// Use the green gradient color function
			FColor Color = GetGradientColor(Intensity);

			ColorData[IX(x, y)] = Color;
		}
	}

	int32 LocalSize = Size;
	ENQUEUE_RENDER_COMMAND(UpdateRenderTarget)(
		[RenderTargetResource, ColorData, LocalSize](FRHICommandListImmediate& RHICmdList)
		{
			FUpdateTextureRegion2D UpdateRegion(0, 0, 0, 0, LocalSize, LocalSize);
			int32 Pitch = LocalSize * sizeof(FColor);
			RHICmdList.UpdateTexture2D(
				RenderTargetResource->GetRenderTargetTexture(), 0, UpdateRegion, Pitch, (uint8*)ColorData.GetData()
			);
		}
		);
}

FColor AFluidGrid::GetGradientColor(float Intensity)
{
	Intensity = FMath::Clamp(Intensity, 0.0f, 1.0f);

	// Define the green gradient colors (example: light green to dark green)
	FVector StartColor(0, 255, 128);  // Light green
	FVector MidColor(0, 204, 102);    // Medium green
	FVector EndColor(0, 153, 76);     // Dark green

	FVector Color;

	if (Intensity < 0.5f)
	{
		Color = FMath::Lerp(StartColor, MidColor, Intensity * 2.0f);
	}
	else
	{
		Color = FMath::Lerp(MidColor, EndColor, (Intensity - 0.5f) * 2.0f);
	}

	return FColor(Color.X, Color.Y, Color.Z, 255);
}



void AFluidGrid::AddDensity(int32 x, int32 y, float amount)
{
	int32 Index = IX(x, y);
	Density[Index] += amount;
}

void AFluidGrid::AddVelocity(int32 x, int32 y, float amountX, float amountY)
{
	int32 Index = IX(x, y);
	Vx[Index] += amountX;
	Vy[Index] += amountY;
}

void AFluidGrid::StepSimulation()
{
	TArray<float> Vx0 = Vx;
	TArray<float> Vy0 = Vy;
	TArray<float> Density0 = Density;

	// Adjusted parameters for improved fluid behavior
	float AdjustedViscosity = Viscosity * 2.0f;
	float AdjustedDiffusion = Diffusion * 2.0f;
	float AdjustedDt = Dt * 2.0f;

	// Diffuse velocities
	Diffuse(1, Vx, Vx0, AdjustedViscosity, AdjustedDt);
	Diffuse(2, Vy, Vy0, AdjustedViscosity, AdjustedDt);

	// Project velocities
	Project(Vx, Vy, Vx0, Vy0);

	// Advect velocities
	Advect(1, Vx, Vx0, Vx0, Vy0, AdjustedDt);
	Advect(2, Vy, Vy0, Vx0, Vy0, AdjustedDt);

	// Project velocities again
	Project(Vx, Vy, Vx0, Vy0);

	// Diffuse densities
	Diffuse(0, Density, Density0, AdjustedDiffusion, AdjustedDt);

	// Advect densities
	Advect(0, Density, Density0, Vx, Vy, AdjustedDt);

	// Apply boundary conditions
	SetBoundary(0, Density);
	SetBoundary(1, Vx);
	SetBoundary(2, Vy);
}


void AFluidGrid::Diffuse(int32 b, TArray<float>& x, TArray<float>& x0, float diff, float dt)
{
	float a = dt * diff * (Size - 2) * (Size - 2);
	LinearSolve(b, x, x0, a, 1 + 4 * a);
}


void AFluidGrid::Advect(int32 b, TArray<float>& d, TArray<float>& d0, TArray<float>& velocX, TArray<float>& velocY, float dt)
{
	float dtx = dt * (Size - 2);
	float dty = dt * (Size - 2);
	float Nfloat = Size - 2;
	int32 i, j;

	for (j = 1; j < Size - 1; j++)
	{
		for (i = 1; i < Size - 1; i++)
		{
			float x = i - dtx * velocX[IX(i, j)];
			float y = j - dty * velocY[IX(i, j)];

			x = FMath::Clamp(x, 0.5f, Nfloat + 0.5f);
			y = FMath::Clamp(y, 0.5f, Nfloat + 0.5f);

			int32 i0 = FMath::FloorToInt(x);
			int32 i1 = i0 + 1;
			int32 j0 = FMath::FloorToInt(y);
			int32 j1 = j0 + 1;

			float s1 = x - i0;
			float s0 = 1.0f - s1;
			float t1 = y - j0;
			float t0 = 1.0f - t1;

			d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) + s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
		}
	}

	SetBoundary(b, d);
}


void AFluidGrid::Project(TArray<float>& velocX, TArray<float>& velocY, TArray<float>& p, TArray<float>& div)
{
	for (int32 j = 1; j < Size - 1; j++)
	{
		for (int32 i = 1; i < Size - 1; i++)
		{
			div[IX(i, j)] = (-0.5f * (velocX[IX(i + 1, j)] - velocX[IX(i - 1, j)] + velocY[IX(i, j + 1)] - velocY[IX(i, j - 1)])) / Size;
			p[IX(i, j)] = 0;
		}
	}

	SetBoundary(0, div);
	SetBoundary(0, p);
	LinearSolve(0, p, div, 1, 6);

	for (int32 j = 1; j < Size - 1; j++)
	{
		for (int32 i = 1; i < Size - 1; i++)
		{
			velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * Size;
			velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * Size;
		}
	}

	SetBoundary(1, velocX);
	SetBoundary(2, velocY);
}

void AFluidGrid::LinearSolve(int32 b, TArray<float>& x, TArray<float>& x0, float a, float c)
{
	float cRecip = 1.0f / c;
	for (int32 t = 0; t < 20; t++)
	{
		for (int32 j = 1; j < Size - 1; j++)
		{
			for (int32 i = 1; i < Size - 1; i++)
			{
				x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i + 1, j)] + x[IX(i - 1, j)] + x[IX(i, j + 1)] + x[IX(i, j - 1)])) * cRecip;
			}
		}
		SetBoundary(b, x);
	}
}

void AFluidGrid::SetBoundary(int32 b, TArray<float>& x)
{
	for (int32 i = 1; i < Size - 1; i++)
	{
		x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, Size - 1)] = b == 2 ? -x[IX(i, Size - 2)] : x[IX(i, Size - 2)];
	}
	for (int32 j = 1; j < Size - 1; j++)
	{
		x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
		x[IX(Size - 1, j)] = b == 1 ? -x[IX(Size - 2, j)] : x[IX(Size - 2, j)];
	}

	x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
	x[IX(0, Size - 1)] = 0.5f * (x[IX(1, Size - 1)] + x[IX(0, Size - 2)]);
	x[IX(Size - 1, 0)] = 0.5f * (x[IX(Size - 2, 0)] + x[IX(Size - 1, 1)]);
	x[IX(Size - 1, Size - 1)] = 0.5f * (x[IX(Size - 2, Size - 1)] + x[IX(Size - 1, Size - 2)]);
}

int32 AFluidGrid::IX(int32 x, int32 y) const
{
	x = FMath::Clamp(x, 0, Size - 1);
	y = FMath::Clamp(y, 0, Size - 1);

	return x + (y * Size);
}
