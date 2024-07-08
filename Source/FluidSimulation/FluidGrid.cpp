#include "FluidGrid.h"
#include "Engine/World.h"
#include "TimerManager.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "UObject/ConstructorHelpers.h"
#include "Components/StaticMeshComponent.h"
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

	static ConstructorHelpers::FObjectFinder<UStaticMesh> PlaneMeshAsset(TEXT("/Game/Path/To/PlaneMesh"));
	if (PlaneMeshAsset.Succeeded())
	{
		PlaneComponent->SetStaticMesh(PlaneMeshAsset.Object);
		PlaneComponent->SetWorldScale3D(FVector(5.0f, 5.0f, 1.0f));
	}
}

void AFluidGrid::InitializeDynamicTexture()
{
	DynamicTexture = UTexture2D::CreateTransient(Size, Size, PF_B8G8R8A8);
	if (DynamicTexture)
	{
		DynamicTexture->CompressionSettings = TC_VectorDisplacementmap;
		DynamicTexture->SRGB = false;
		DynamicTexture->AddToRoot();
		DynamicTexture->UpdateResource();
		UE_LOG(LogTemp, Log, TEXT("DynamicTexture created successfully with size %d x %d."), Size, Size);
	}
	else
	{
		UE_LOG(LogTemp, Error, TEXT("Failed to create DynamicTexture."));
	}
}
void AFluidGrid::BeginPlay()
{
	Super::BeginPlay();

	if (!BaseMaterial)
	{
		UE_LOG(LogTemp, Error, TEXT("BaseMaterial is not assigned."));
		return;
	}

	InitializeDynamicTexture();

	DynamicMaterialInstance = UMaterialInstanceDynamic::Create(BaseMaterial, this);
	if (DynamicMaterialInstance)
	{
		DynamicMaterialInstance->SetTextureParameterValue(FName("DynamicTexture"), DynamicTexture);
		PlaneComponent->SetMaterial(0, DynamicMaterialInstance);
	}
	else
	{
		UE_LOG(LogTemp, Error, TEXT("Failed to create dynamic material instance."));
	}

	GetWorld()->GetTimerManager().SetTimer(TimerHandle, this, &AFluidGrid::StepSimulation, Dt, true);
}

void AFluidGrid::UpdateTexture()
{
	if (!DynamicTexture)
	{
		UE_LOG(LogTemp, Error, TEXT("DynamicTexture is not initialized."));
		return;
	}

	FTexture2DMipMap& Mip = DynamicTexture->GetPlatformData()->Mips[0];
	void* Data = Mip.BulkData.Lock(LOCK_READ_WRITE);

	TArray<FColor> ColorData;
	ColorData.SetNum(Size * Size);

	for (int32 y = 0; y < Size; y++)
	{
		for (int32 x = 0; x < Size; x++)
		{
			float Value = Density[IX(x, y)];
			uint8 Intensity = FMath::Clamp(Value * 255.0f, 0.0f, 255.0f);
			ColorData[IX(x, y)] = FColor(Intensity, Intensity, Intensity, 255); // Ensure RGBA format
		}
	}

	FMemory::Memcpy(Data, ColorData.GetData(), ColorData.Num() * sizeof(FColor));
	Mip.BulkData.Unlock();
	DynamicTexture->UpdateResource();
	UE_LOG(LogTemp, Log, TEXT("DynamicTexture updated successfully."));
}

void AFluidGrid::AddDensity(int32 x, int32 y, float amount)
{
	Density[IX(x, y)] += amount;
}

void AFluidGrid::AddVelocity(int32 x, int32 y, float amountX, float amountY)
{
	int32 index = IX(x, y);
	Vx[index] += amountX;
	Vy[index] += amountY;
}

void AFluidGrid::StepSimulation()
{
	AddDensity(Size / 2, Size / 2, 100.0f);
	AddVelocity(Size / 2, Size / 2, 1.0f, 0.0f);

	TArray<float> Vx0 = Vx;
	TArray<float> Vy0 = Vy;
	TArray<float> Density0 = Density;

	Diffuse(1, Vx, Vx0, Viscosity, Dt);
	Diffuse(2, Vy, Vy0, Viscosity, Dt);

	Project(Vx, Vy, Vx0, Vy0);

	Advect(1, Vx, Vx0, Vx0, Vy0, Dt);
	Advect(2, Vy, Vy0, Vx0, Vy0, Dt);

	Project(Vx, Vy, Vx0, Vy0);

	Diffuse(0, Density, Density0, Diffusion, Dt);

	Advect(0, Density, Density0, Vx, Vy, Dt);

	UpdateTexture();
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
	for (int32 t = 0; t < Size; t++)
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