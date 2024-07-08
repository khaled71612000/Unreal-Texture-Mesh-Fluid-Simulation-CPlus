#include "FluidGrid.h"
#include "Engine/World.h"
#include "TimerManager.h"
#include "Math/UnrealMathUtility.h"
#include "UObject/ConstructorHelpers.h"
#include "../../../../../../../Plugins/Runtime/ProceduralMeshComponent/Source/ProceduralMeshComponent/Public/ProceduralMeshComponent.h"

AFluidGrid::AFluidGrid()
{
	PrimaryActorTick.bCanEverTick = true;

	Size = 10; // Adjust size as needed
	Dt = 0.1f;
	Diffusion = 1.f;
	Viscosity = 1.f;

	int32 TotalSize = Size * Size * Size;
	S.SetNumZeroed(TotalSize);
	Density.SetNumZeroed(TotalSize);
	Vx.SetNumZeroed(TotalSize);
	Vy.SetNumZeroed(TotalSize);
	Vz.SetNumZeroed(TotalSize);
	Vx0.SetNumZeroed(TotalSize);
	Vy0.SetNumZeroed(TotalSize);
	Vz0.SetNumZeroed(TotalSize);

	ProceduralMesh = CreateDefaultSubobject<UProceduralMeshComponent>(TEXT("ProceduralMesh"));
	RootComponent = ProceduralMesh;

	ProceduralMesh->SetVisibility(true);
	ProceduralMesh->SetCollisionEnabled(ECollisionEnabled::NoCollision);
}

void AFluidGrid::BeginPlay()
{
	Super::BeginPlay();
	GenerateGrid();

	// Schedule recurring updates for the simulation
	GenerateGrid();

	// Assign a basic material to the mesh
	static ConstructorHelpers::FObjectFinder<UMaterial> Material(TEXT("/Game/StarterContent/Materials/M_Basic_Wall.M_Basic_Wall"));
	if (Material.Succeeded())
	{
		ProceduralMesh->SetMaterial(0, Material.Object);
	}

	// Schedule recurring updates for the simulation
	GetWorld()->GetTimerManager().SetTimer(TimerHandle, this, &AFluidGrid::StepSimulation, Dt, true);
}
	// Initialize a small portion of the Density array
	//Density[0] = 0.5f;

	//if (FluidMesh)
	//{
	//	UMaterialInstance* MaterialInstance = Cast<UMaterialInstance>(StaticLoadObject(UMaterialInstance::StaticClass(), nullptr, TEXT("/Game/FluidSimulation/MI_Fluid.MI_Fluid")));
	//	if (MaterialInstance)
	//	{
	//		DynamicMaterial = UMaterialInstanceDynamic::Create(MaterialInstance, FluidMesh);
	//		if (DynamicMaterial)
	//		{
	//			DynamicMaterial->SetTextureParameterValue(FName("FluidTexture"), DynamicTexture);
	//			FluidMesh->SetMaterial(0, DynamicMaterial);
	//			UE_LOG(LogTemp, Warning, TEXT("DynamicMaterial successfully created."));
	//		}
	//		else
	//		{
	//			UE_LOG(LogTemp, Error, TEXT("Failed to create DynamicMaterial."));
	//		}
	//	}
	//	else
	//	{
	//		UE_LOG(LogTemp, Error, TEXT("Failed to load the material instance."));
	//	}
	//}
	//else
	//{
	//	UE_LOG(LogTemp, Error, TEXT("FluidMesh is null."));
	//}

	//UpdateMaterialParameters();

	//// Schedule recurring updates for the simulation
	//GetWorld()->GetTimerManager().SetTimer(TimerHandle, this, &AFluidGrid::StepSimulation, Dt, true);


void AFluidGrid::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);
}

void AFluidGrid::AddDensity(int32 x, int32 y, int32 z, float amount)
{
	int32 index = IX(x, y, z);
	if (!Density.IsValidIndex(index))
	{
		return;
	}
	Density[index] += amount;
}

void AFluidGrid::AddVelocity(int32 x, int32 y, int32 z, float amountX, float amountY, float amountZ)
{
	int32 index = IX(x, y, z);
	if (Vx.IsValidIndex(index) && Vy.IsValidIndex(index) && Vz.IsValidIndex(index))
	{
		Vx[index] += amountX;
		Vy[index] += amountY;
		Vz[index] += amountZ;
	}
}

int32 AFluidGrid::IX(int32 x, int32 y, int32 z)
{
	return x + (y * Size) + (z * Size * Size);
}

void AFluidGrid::diffuse(int32 b, TArray<float>& x, TArray<float>& x0, float diff, float dt)
{
	float a = dt * diff * (Size - 2) * (Size - 2);
	LinearSolve(b, x, x0, a, 1 + 6 * a);
}

void AFluidGrid::LinearSolve(int32 b, TArray<float>& xArray, const TArray<float>& x0, float a, float c)
{
	float cRecip = 1.0f / c;
	for (int32 k = 0; k < Iteration; k++) {
		for (int32 z = 1; z < Size - 1; z++) {
			for (int32 y = 1; y < Size - 1; y++) {
				for (int32 x = 1; x < Size - 1; x++) {
					xArray[IX(x, y, z)] =
						(x0[IX(x, y, z)]
							+ a * (xArray[IX(x + 1, y, z)]
								+ xArray[IX(x - 1, y, z)]
								+ xArray[IX(x, y + 1, z)]
								+ xArray[IX(x, y - 1, z)]
								+ xArray[IX(x, y, z + 1)]
								+ xArray[IX(x, y, z - 1)]
								)) * cRecip;
				}
			}
		}
		SetBoundary(b, xArray);
	}
}
void AFluidGrid::SetBoundary(int32 b, TArray<float>& x)
{
	for (int32 j = 1; j < Size - 1; j++) {
		for (int32 i = 1; i < Size - 1; i++) {
			x[IX(i, j, 0)] = b == 3 ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
			x[IX(i, j, Size - 1)] = b == 3 ? -x[IX(i, j, Size - 2)] : x[IX(i, j, Size - 2)];
		}
	}
	for (int32 k = 1; k < Size - 1; k++) {
		for (int32 i = 1; i < Size - 1; i++) {
			x[IX(i, 0, k)] = b == 2 ? -x[IX(i, 1, k)] : x[IX(i, 1, k)];
			x[IX(i, Size - 1, k)] = b == 2 ? -x[IX(i, Size - 2, k)] : x[IX(i, Size - 2, k)];
		}
	}
	for (int32 k = 1; k < Size - 1; k++) {
		for (int32 j = 1; j < Size - 1; j++) {
			x[IX(0, j, k)] = b == 1 ? -x[IX(1, j, k)] : x[IX(1, j, k)];
			x[IX(Size - 1, j, k)] = b == 1 ? -x[IX(Size - 2, j, k)] : x[IX(Size - 2, j, k)];
		}
	}

	// Corner cases handling
	x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)] + x[IX(0, 1, 0)] + x[IX(0, 0, 1)]);
	x[IX(0, Size - 1, 0)] = 0.33f * (x[IX(1, Size - 1, 0)] + x[IX(0, Size - 2, 0)] + x[IX(0, Size - 1, 1)]);
	x[IX(0, 0, Size - 1)] = 0.33f * (x[IX(1, 0, Size - 1)] + x[IX(0, 1, Size - 1)] + x[IX(0, 0, Size - 2)]);
	x[IX(0, Size - 1, Size - 1)] = 0.33f * (x[IX(1, Size - 1, Size - 1)] + x[IX(0, Size - 2, Size - 1)] + x[IX(0, Size - 1, Size - 2)]);
	x[IX(Size - 1, 0, 0)] = 0.33f * (x[IX(Size - 2, 0, 0)] + x[IX(Size - 1, 1, 0)] + x[IX(Size - 1, 0, 1)]);
	x[IX(Size - 1, Size - 1, 0)] = 0.33f * (x[IX(Size - 2, Size - 1, 0)] + x[IX(Size - 1, Size - 2, 0)] + x[IX(Size - 1, Size - 1, 1)]);
	x[IX(Size - 1, 0, Size - 1)] = 0.33f * (x[IX(Size - 2, 0, Size - 1)] + x[IX(Size - 1, 1, Size - 1)] + x[IX(Size - 1, 0, Size - 2)]);
	x[IX(Size - 1, Size - 1, Size - 1)] = 0.33f * (x[IX(Size - 2, Size - 1, Size - 1)] + x[IX(Size - 1, Size - 2, Size - 1)] + x[IX(Size - 1, Size - 1, Size - 2)]);
}

void AFluidGrid::project(TArray<float>& velocX, TArray<float>& velocY, TArray<float>& velocZ, TArray<float>& p, TArray<float>& div)
{
	for (int32 z = 1; z < Size - 1; z++) {
		for (int32 y = 1; y < Size - 1; y++) {
			for (int32 x = 1; x < Size - 1; x++) {
				div[IX(x, y, z)] = -0.5f * (
					velocX[IX(x + 1, y, z)] - velocX[IX(x - 1, y, z)] +
					velocY[IX(x, y + 1, z)] - velocY[IX(x, y - 1, z)] +
					velocZ[IX(x, y, z + 1)] - velocZ[IX(x, y, z - 1)]
					) / Size;
				p[IX(x, y, z)] = 0;
			}
		}
	}

	SetBoundary(0, div);
	SetBoundary(0, p);
	LinearSolve(0, p, div, 1, 6);

	for (int32 z = 1; z < Size - 1; z++) {
		for (int32 y = 1; y < Size - 1; y++) {
			for (int32 x = 1; x < Size - 1; x++) {
				velocX[IX(x, y, z)] -= 0.5f * (p[IX(x + 1, y, z)] - p[IX(x - 1, y, z)]) * Size;
				velocY[IX(x, y, z)] -= 0.5f * (p[IX(x, y + 1, z)] - p[IX(x, y - 1, z)]) * Size;
				velocZ[IX(x, y, z)] -= 0.5f * (p[IX(x, y, z + 1)] - p[IX(x, y, z - 1)]) * Size;
			}
		}
	}

	SetBoundary(1, velocX);
	SetBoundary(2, velocY);
	SetBoundary(3, velocZ);
}


void AFluidGrid::Advect(int32 b, TArray<float>& d, const TArray<float>& d0, const TArray<float>& velocX, const TArray<float>& velocY, const TArray<float>& velocZ, float dt)
{
	float i0, i1, j0, j1, k0, k1;
	float dtx = dt * (Size - 2);
	float dty = dt * (Size - 2);
	float dtz = dt * (Size - 2);
	float s0, s1, t0, t1, u0, u1;
	float tmp1, tmp2, tmp3, x, y, z;
	float Nfloat = Size;
	float ifloat, jfloat, kfloat;
	int32 i, j, k;

	for (k = 1, kfloat = 1; k < Size - 1; k++, kfloat++) {
		for (j = 1, jfloat = 1; j < Size - 1; j++, jfloat++) {
			for (i = 1, ifloat = 1; i < Size - 1; i++, ifloat++) {
				tmp1 = dtx * velocX[IX(i, j, k)];
				tmp2 = dty * velocY[IX(i, j, k)];
				tmp3 = dtz * velocZ[IX(i, j, k)];
				x = ifloat - tmp1;
				y = jfloat - tmp2;
				z = kfloat - tmp3;

				if (x < 0.5f) x = 0.5f;
				if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
				i0 = floorf(x);
				i1 = i0 + 1.0f;
				if (y < 0.5f) y = 0.5f;
				if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
				j0 = floorf(y);
				j1 = j0 + 1.0f;
				if (z < 0.5f) z = 0.5f;
				if (z > Nfloat + 0.5f) z = Nfloat + 0.5f;
				k0 = floorf(z);
				k1 = k0 + 1.0f;

				s1 = x - i0;
				s0 = 1.0f - s1;
				t1 = y - j0;
				t0 = 1.0f - t1;
				u1 = z - k0;
				u0 = 1.0f - u1;

				int32 i0i = static_cast<int32>(i0);
				int32 i1i = static_cast<int32>(i1);
				int32 j0i = static_cast<int32>(j0);
				int32 j1i = static_cast<int32>(j1);
				int32 k0i = static_cast<int32>(k0);
				int32 k1i = static_cast<int32>(k1);

				d[IX(i, j, k)] =
					s0 * (t0 * (u0 * d0[IX(i0i, j0i, k0i)]
						+ u1 * d0[IX(i0i, j0i, k1i)])
						+ (t1 * (u0 * d0[IX(i0i, j1i, k0i)]
							+ u1 * d0[IX(i0i, j1i, k1i)])))
					+ s1 * (t0 * (u0 * d0[IX(i1i, j0i, k0i)]
						+ u1 * d0[IX(i1i, j0i, k1i)])
						+ (t1 * (u0 * d0[IX(i1i, j1i, k0i)]
							+ u1 * d0[IX(i1i, j1i, k1i)])));
			}
		}
	}

	SetBoundary(b, d);
}
void AFluidGrid::GenerateGrid()
{
	TArray<FVector> Vertices;
	TArray<int32> Triangles;
	TArray<FVector> Normals;
	TArray<FVector2D> UV0;
	TArray<FLinearColor> VertexColors;
	TArray<FProcMeshTangent> Tangents;

	for (int32 y = 0; y < Size; y++)
	{
		for (int32 x = 0; x < Size; x++)
		{
			int32 Index = x + y * Size;

			Vertices.Add(FVector(x, y, 0));
			Vertices.Add(FVector(x + 1, y, 0));
			Vertices.Add(FVector(x + 1, y + 1, 0));
			Vertices.Add(FVector(x, y + 1, 0));

			Triangles.Add(Index * 4);
			Triangles.Add(Index * 4 + 1);
			Triangles.Add(Index * 4 + 2);

			Triangles.Add(Index * 4);
			Triangles.Add(Index * 4 + 2);
			Triangles.Add(Index * 4 + 3);

			VertexColors.Add(FLinearColor::White);
			VertexColors.Add(FLinearColor::White);
			VertexColors.Add(FLinearColor::White);
			VertexColors.Add(FLinearColor::White);
		}
	}

	ProceduralMesh->CreateMeshSection_LinearColor(0, Vertices, Triangles, Normals, UV0, VertexColors, Tangents, true);
	UE_LOG(LogTemp, Warning, TEXT("Grid Generated: %d vertices, %d triangles"), Vertices.Num(), Triangles.Num());
}


void AFluidGrid::UpdateMesh()
{
	TArray<FVector> Vertices;
	for (int32 y = 0; y < Size; y++)
	{
		for (int32 x = 0; x < Size; x++)
		{
			int32 Index = x + y * Size;
			float DensityValue = Density[Index];

			Vertices.Add(FVector(x, y, DensityValue * 10.0f));
			Vertices.Add(FVector(x + 1, y, DensityValue * 10.0f));
			Vertices.Add(FVector(x + 1, y + 1, DensityValue * 10.0f));
			Vertices.Add(FVector(x, y + 1, DensityValue * 10.0f));
		}
	}

	ProceduralMesh->UpdateMeshSection_LinearColor(0, Vertices, TArray<FVector>(), TArray<FVector2D>(), TArray<FLinearColor>(), TArray<FProcMeshTangent>());

	// Logging vertex positions for debugging
	for (int32 i = 0; i < Vertices.Num(); i++)
	{
		UE_LOG(LogTemp, Warning, TEXT("Vertex %d: %s"), i, *Vertices[i].ToString());
	}
}



void AFluidGrid::StepSimulation()
{
	AddDensity(1, 1, 0, 1.0f);
	AddVelocity(1, 1, 0, 1.0f, 0.0f, 0.0f);

	diffuse(1, Vx, Vx0, Viscosity, Dt);
	diffuse(2, Vy, Vy0, Viscosity, Dt);
	diffuse(3, Vz, Vz0, Viscosity, Dt);

	project(Vx, Vy, Vz, Vx0, Vy0);

	Advect(1, Vx, Vx0, Vx0, Vy0, Vz0, Dt);
	Advect(2, Vy, Vy0, Vx0, Vy0, Vz0, Dt);
	Advect(3, Vz, Vz0, Vx0, Vy0, Vz0, Dt);

	project(Vx, Vy, Vz, Vx0, Vy0);

	UpdateMesh();
}


//
//void AFluidGrid::StepSimulation()
//{
//	// Apply forces like vorticity confinement and noise
//	AddVorticityConfinement(Dt);
//	AddNoise(Dt);
//
//	// Existing simulation steps
//	Density[0] = 0.5f + 0.5f * FMath::Sin(GetWorld()->GetTimeSeconds());
//	UpdateMaterialParameters();
//	UpdateTexture();
//
//	UE_LOG(LogTemp, Warning, TEXT("Density[0]: %f"), Density[0]);
//}
//
//void AFluidGrid::ComputeVorticity()
//{
//	for (int32 z = 1; z < Size - 1; z++)
//	{
//		for (int32 y = 1; y < Size - 1; y++)
//		{
//			for (int32 x = 1; x < Size - 1; x++)
//			{
//				int32 index = IX(x, y, z);
//				float dx = Vz[IX(x, y + 1, z)] - Vz[IX(x, y - 1, z)];
//				float dy = Vx[IX(x, y, z + 1)] - Vx[IX(x, y, z - 1)];
//				float dz = Vy[IX(x + 1, y, z)] - Vy[IX(x - 1, y, z)];
//				Vorticity[index] = FVector(dx, dy, dz);
//			}
//		}
//	}
//}
//
//void AFluidGrid::AddVorticityConfinement(float dt)
//{
//	ComputeVorticity();
//	for (int32 z = 1; z < Size - 1; z++)
//	{
//		for (int32 y = 1; y < Size - 1; y++)
//		{
//			for (int32 x = 1; x < Size - 1; x++)
//			{
//				int32 index = IX(x, y, z);
//				FVector vorticityForce = Vorticity[index].GetClampedToMaxSize(MaxVorticity);
//				Vx[index] += vorticityForce.X * dt;
//				Vy[index] += vorticityForce.Y * dt;
//				Vz[index] += vorticityForce.Z * dt;
//			}
//		}
//	}
//}
//
//void AFluidGrid::AddNoise(float dt)
//{
//	for (int32 z = 0; z < Size; z++)
//	{
//		for (int32 y = 0; y < Size; y++)
//		{
//			for (int32 x = 0; x < Size; x++)
//			{
//				int32 index = IX(x, y, z);
//				float noiseX = FMath::PerlinNoise3D(FVector(x, y, z) * NoiseScale);
//				float noiseY = FMath::PerlinNoise3D(FVector(x + 100, y, z) * NoiseScale); // Offset to get different noise
//				float noiseZ = FMath::PerlinNoise3D(FVector(x, y + 100, z) * NoiseScale); // Offset to get different noise
//				Vx[index] += noiseX * dt;
//				Vy[index] += noiseY * dt;
//				Vz[index] += noiseZ * dt;
//			}
//		}
//	}
//}
//
//void AFluidGrid::UpdateMaterialParameters()
//{
//	if (DynamicMaterial)
//	{
//		float DensityValue = FMath::Clamp(Density[0], 0.0f, 1.0f);
//		DynamicMaterial->SetScalarParameterValue(FName("FluidDensity"), DensityValue);
//
//		FLinearColor ColorValue = FLinearColor(DensityValue, 0.0f, 1.0f - DensityValue);
//		DynamicMaterial->SetVectorParameterValue(FName("FluidColor"), ColorValue);
//
//		UE_LOG(LogTemp, Warning, TEXT("Material - FluidDensity: %f, FluidColor: %s"), DensityValue, *ColorValue.ToString());
//	}
//	else
//	{
//		UE_LOG(LogTemp, Error, TEXT("DynamicMaterial is null."));
//	}
//}
//
//void AFluidGrid::UpdateTexture()
//{
//	// Ensure the texture is valid and has a mipmap
//	if (!DynamicTexture || DynamicTexture->GetPlatformData()->Mips.Num() == 0)
//	{
//		UE_LOG(LogTemp, Error, TEXT("DynamicTexture is invalid or has no mips."));
//		return;
//	}
//
//	// Get the first mipmap level
//	FTexture2DMipMap& Mip = DynamicTexture->GetPlatformData()->Mips[0];
//
//	// Lock the texture so we can update its data
//	void* Data = Mip.BulkData.Lock(LOCK_READ_WRITE);
//
//	FColor* FormattedImageData = static_cast<FColor*>(Data);
//
//	for (int32 y = 0; y < Size; y++)
//	{
//		for (int32 x = 0; x < Size; x++)
//		{
//			float DensityValue = FMath::Clamp(Density[IX(x, y, 0)], 0.0f, 1.0f);
//			FColor Color = FLinearColor::LerpUsingHSV(FLinearColor::Blue, FLinearColor::Red, DensityValue).ToFColor(true);
//			FormattedImageData[x + y * Size] = Color;
//		}
//	}
//
//	Mip.BulkData.Unlock();
//	DynamicTexture->UpdateResource();
//}