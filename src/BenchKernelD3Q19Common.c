// --------------------------------------------------------------------------
//
// Copyright
//   Markus Wittmann, 2016-2017
//   RRZE, University of Erlangen-Nuremberg, Germany
//   markus.wittmann -at- fau.de or hpc -at- rrze.fau.de
//
//   Viktor Haag, 2016
//   LSS, University of Erlangen-Nuremberg, Germany
//
//  This file is part of the Lattice Boltzmann Benchmark Kernels (LbmBenchKernels).
//
//  LbmBenchKernels is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  LbmBenchKernels is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with LbmBenchKernels.  If not, see <http://www.gnu.org/licenses/>.
//
// --------------------------------------------------------------------------
#include "BenchKernelD3Q19Common.h"

#include "Memory.h"
#include "Vtk.h"

#include <inttypes.h>
#include <math.h>


// Forward definition.
void FNAME(D3Q19Kernel)(LatticeDesc * ld, struct KernelData_ * kd, CaseData * cd);

void FNAME(D3Q19BlkKernel)(LatticeDesc * ld, struct KernelData_ * kd, CaseData * cd);



static void FNAME(BcGetPdf)(KernelData * kd, int x, int y, int z, int dir, PdfT * pdf)
{
	Assert(kd != NULL);
	Assert(kd->PdfsActive != NULL);
	Assert(kd->PdfsActive == kd->Pdfs[0] || kd->PdfsActive == kd->Pdfs[1]);
	Assert(pdf != NULL);

	Assert(x >= 0);
	Assert(y >= 0);
	Assert(z >= 0);
	Assert(x < kd->Dims[0]);
	Assert(y < kd->Dims[1]);
	Assert(z < kd->Dims[2]);
	Assert(dir >= 0);
	Assert(dir < N_D3Q19);

	int oX = kd->Offsets[0];
	int oY = kd->Offsets[1];
	int oZ = kd->Offsets[2];

#ifdef PROP_MODEL_PUSH
	int nx = x;
	int ny = y;
	int nz = z;
#elif PROP_MODEL_PULL
	int nx = x - D3Q19_X[dir];
	int ny = y - D3Q19_Y[dir];
	int nz = z - D3Q19_Z[dir];
#endif

	#define I(x, y, z, dir) P_INDEX_5(kd->GlobalDims, (x), (y), (z), (dir))
	*pdf = kd->PdfsActive[I(nx + oX, ny + oY, nz + oZ, dir)];
	#undef I

	return;
}

static void FNAME(BcSetPdf)(KernelData * kd, int x, int y, int z, int dir, PdfT pdf)
{
	Assert(kd != NULL);
	Assert(kd->PdfsActive != NULL);
	Assert(kd->PdfsActive == kd->Pdfs[0] || kd->PdfsActive == kd->Pdfs[1]);
	Assert(x >= 0);
	Assert(y >= 0);
	Assert(z >= 0);
	Assert(x < kd->Dims[0]);
	Assert(y < kd->Dims[1]);
	Assert(z < kd->Dims[2]);
	Assert(dir >= 0);
	Assert(dir < N_D3Q19);

	int oX = kd->Offsets[0];
	int oY = kd->Offsets[1];
	int oZ = kd->Offsets[2];

#ifdef PROP_MODEL_PUSH
	int nx = x;
	int ny = y;
	int nz = z;
#elif PROP_MODEL_PULL
	int nx = x - D3Q19_X[dir];
	int ny = y - D3Q19_Y[dir];
	int nz = z - D3Q19_Z[dir];
#endif

	#define I(x, y, z, dir) P_INDEX_5(kd->GlobalDims, (x), (y), (z), (dir))
	kd->PdfsActive[I(nx + oX, ny + oY, nz + oZ, dir)] = pdf;
	#undef I


	return;
}


static void FNAME(GetNode)(KernelData * kd, int x, int y, int z, PdfT * pdfs)
{
	Assert(kd != NULL);
	Assert(kd->PdfsActive != NULL);
	Assert(kd->PdfsActive == kd->Pdfs[0] || kd->PdfsActive == kd->Pdfs[1]);
	Assert(pdfs != NULL);
	Assert(x >= 0);
	Assert(y >= 0);
	Assert(z >= 0);
	Assert(x < kd->Dims[0]);
	Assert(y < kd->Dims[1]);
	Assert(z < kd->Dims[2]);

	int oX = kd->Offsets[0];
	int oY = kd->Offsets[1];
	int oZ = kd->Offsets[2];


	#define I(x, y, z, dir) P_INDEX_5(kd->GlobalDims, (x), (y), (z), (dir))
#ifdef PROP_MODEL_PUSH
	#define X(name, idx, idxinv, _x, _y, _z)	pdfs[idx] = kd->PdfsActive[I(x + oX, y + oY, z + oZ, idx)];
#elif PROP_MODEL_PULL
	#define X(name, idx, idxinv, _x, _y, _z)	pdfs[idx] = kd->PdfsActive[I(x + oX - (_x), y + oY - (_y), z + oZ - (_z), idx)];
#endif
	D3Q19_LIST
	#undef X
	#undef I

#if 0		// DETECT NANs

	for (int d = 0; d < 19; ++d) {
		if (isnan(pdfs[d])) {
			printf("%d %d %d %d nan! get node\n", x, y, z, d);

			for (int d2 = 0; d2 < 19; ++d2) {
				printf("%d: %e\n", d2, pdfs[d2]);
			}

			exit(1);
		}
	}

#endif

	return;
}


static void FNAME(SetNode)(KernelData * kd, int x, int y, int z, PdfT * pdfs)
{
	Assert(kd != NULL);
	Assert(kd->PdfsActive != NULL);
	Assert(kd->PdfsActive == kd->Pdfs[0] || kd->PdfsActive == kd->Pdfs[1]);
	Assert(pdfs != NULL);

	Assert(x >= 0);
	Assert(y >= 0);
	Assert(z >= 0);
	Assert(x < kd->Dims[0]);
	Assert(y < kd->Dims[1]);
	Assert(z < kd->Dims[2]);

	int oX = kd->Offsets[0];
	int oY = kd->Offsets[1];
	int oZ = kd->Offsets[2];

	#define I(x, y, z, dir) P_INDEX_5(kd->GlobalDims, (x), (y), (z), (dir))
#ifdef PROP_MODEL_PUSH
	#define X(name, idx, idxinv, _x, _y, _z)	kd->PdfsActive[I(x + oX, y + oY, z + oZ, idx)] = pdfs[idx];
#elif PROP_MODEL_PULL
	#define X(name, idx, idxinv, _x, _y, _z)	kd->PdfsActive[I(x + oX - (_x), y + oY - (_y), z + oZ - (_z), idx)] = pdfs[idx];
#endif
	D3Q19_LIST
	#undef X
	#undef I

	return;
}


static void ParameterUsage()
{
	printf("Kernel parameters:\n");
	printf("  [-blk <n>] [-blk-[xyz] <n>]\n");

	return;
}

static void ParseParameters(Parameters * params, int * blk)
{
	Assert(blk != NULL);

	blk[0] = 0; blk[1] = 0; blk[2] = 0;

	#define ARG_IS(param)                   (!strcmp(params->KernelArgs[i], param))
	#define NEXT_ARG_PRESENT() \
		do { \
			if (i + 1 >= params->nKernelArgs) { \
				printf("ERROR: argument %s requires a parameter.\n", params->KernelArgs[i]); \
				exit(1); \
			} \
		} while (0)


	for (int i = 0; i < params->nKernelArgs; ++i) {
		if (ARG_IS("-blk") || ARG_IS("--blk")) {
			NEXT_ARG_PRESENT();

			int tmp = strtol(params->KernelArgs[++i], NULL, 0);

			if (tmp <= 0) {
				printf("ERROR: blocking parameter must be > 0.\n");
				exit(1);
			}

			blk[0] = blk[1] = blk[2] = tmp;
		}
		else if (ARG_IS("-blk-x") || ARG_IS("--blk-x")) {
			NEXT_ARG_PRESENT();

			int tmp = strtol(params->KernelArgs[++i], NULL, 0);

			if (tmp <= 0) {
				printf("ERROR: blocking parameter must be > 0.\n");
				exit(1);
			}

			blk[0] = tmp;
		}
		else if (ARG_IS("-blk-y") || ARG_IS("--blk-y")) {
			NEXT_ARG_PRESENT();

			int tmp = strtol(params->KernelArgs[++i], NULL, 0);

			if (tmp <= 0) {
				printf("ERROR: blocking parameter must be > 0.\n");
				exit(1);
			}

			blk[1] = tmp;
		}
		else if (ARG_IS("-blk-z") || ARG_IS("--blk-z")) {
			NEXT_ARG_PRESENT();

			int tmp = strtol(params->KernelArgs[++i], NULL, 0);

			if (tmp <= 0) {
				printf("ERROR: blocking parameter must be > 0.\n");
				exit(1);
			}

			blk[2] = tmp;
		}
		else if (ARG_IS("-h") || ARG_IS("-help") || ARG_IS("--help")) {
			ParameterUsage();
			exit(1);
		}
		else {
			printf("ERROR: unknown kernel parameter.\n");
			ParameterUsage();
			exit(1);
		}
	}

	#undef ARG_IS
	#undef NEXT_ARG_PRESENT

	return;
}


void FNAME(D3Q19BlkInit)(LatticeDesc * ld, KernelData ** kernelData, Parameters * params)
{
	KernelDataEx * kdex = NULL;
	MemAlloc((void **)&kdex, sizeof(KernelDataEx));

	kdex->Blk[0] = 0; kdex->Blk[1] = 0; kdex->Blk[2] = 0;

	KernelData * kd = &kdex->kd;
	*kernelData = kd;

	kd->nObstIndices = ld->nObst;

	// Ajust the dimensions according to padding, if used.
	kd->Dims[0] = ld->Dims[0];
	kd->Dims[1] = ld->Dims[1];
	kd->Dims[2] = ld->Dims[2];


	int * lDims = ld->Dims;
	int * gDims = kd->GlobalDims;

	gDims[0] = lDims[0] + 2;
	gDims[1] = lDims[1] + 2;
	gDims[2] = lDims[2] + 2;

	kd->Offsets[0] = 1;
	kd->Offsets[1] = 1;
	kd->Offsets[2] = 1;

	int lX = lDims[0];
	int lY = lDims[1];
	int lZ = lDims[2];

	int gX = gDims[0];
	int gY = gDims[1];
	int gZ = gDims[2];

	int oX = kd->Offsets[0];
	int oY = kd->Offsets[1];
	int oZ = kd->Offsets[2];

	int blk[3] = { 0 };

	int nCells = gX * gY * gZ;

	PdfT * pdfs[2];

	ParseParameters(params, blk);

	if (blk[0] == 0) blk[0] = gX;
	if (blk[1] == 0) blk[1] = gY;
	if (blk[2] == 0) blk[2] = gZ;

	printf("# blocking x: %3d y: %3d z: %3d\n", blk[0], blk[1], blk[2]);


	kdex->Blk[0] = blk[0]; kdex->Blk[1] = blk[1]; kdex->Blk[2] = blk[2];


	printf("# allocating data for %d LB nodes with padding (%lu bytes = %f MiB for both lattices)\n",
	       nCells, 2 * sizeof(PdfT) * nCells * N_D3Q19,
	       2 * sizeof(PdfT) * nCells * N_D3Q19 / 1024.0 / 1024.0);

	MemAlloc((void **)&pdfs[0], sizeof(PdfT) * nCells * N_D3Q19);
	MemAlloc((void **)&pdfs[1], sizeof(PdfT) * nCells * N_D3Q19);

	kd->Pdfs[0] = pdfs[0];
	kd->Pdfs[1] = pdfs[1];

	// Initialize PDFs with some (arbitrary) data for correct NUMA placement.
	// This depends on the chosen data layout.
	// The structure of the loop should resemble the same "execution layout"
	// as in the kernel!
#ifdef _OPENMP
	#pragma omp parallel for collapse(3)
#endif

	for (int bZ = 0; bZ < gZ; bZ += blk[2]) {
		for (int bY = 0; bY < gY; bY += blk[1]) {
			for (int bX = 0; bX < gX; bX += blk[0]) {

				// Must do everything here, else it would break collapse.
				int eZ = MIN(bZ + blk[2], gZ);
				int eY = MIN(bY + blk[1], gY);
				int eX = MIN(bX + blk[0], gX);

				for (int z = bZ; z < eZ; ++z) {
					for (int y = bY; y < eY; ++y) {
						for (int x = bX; x < eX; ++x) {

							for (int d = 0; d < N_D3Q19; ++d) {
								pdfs[0][P_INDEX_5(gDims, x, y, z, d)] = 1.0;
								pdfs[1][P_INDEX_5(gDims, x, y, z, d)] = 1.0;
							}

						}
					}
				}
			}
		}
	}

	// Initialize all PDFs to some standard value.
	for (int z = 0; z < gZ; ++z) {
		for (int y = 0; y < gY; ++y) {
			for (int x = 0; x < gX; ++x) {
				for (int d = 0; d < N_D3Q19; ++d) {
					pdfs[0][P_INDEX_5(gDims, x, y, z, d)] = 0.0;
					pdfs[1][P_INDEX_5(gDims, x, y, z, d)] = 0.0;
				}
			}
		}
	}


	// Count how many *PDFs* need bounce back treatment.

	uint64_t nPdfs = ((uint64_t)19) * gX * gY * gZ;

	if (nPdfs > ((2LU << 31) - 1)) {
		printf("ERROR: number of PDFs exceed 2^31.\n");
		exit(1);
	}

	// Compiler bug? Incorrect computation of nBounceBackPdfs when using icc 15.0.2.
	// Works when declaring nBounceBackPdfs as int64_t or using volatile.
	volatile int nBounceBackPdfs = 0;
	// int64_t nBounceBackPdfs = 0;
	int nx, ny, nz, px, py, pz;

	// TODO: apply blocking?

	for (int z = 0; z < lZ; ++z) {
		for (int y = 0; y < lY; ++y) {
			for (int x = 0; x < lX; ++x) {

				if (ld->Lattice[L_INDEX_4(ld->Dims, x, y, z)] != LAT_CELL_OBSTACLE) {
					for (int d = 0; d < N_D3Q19; ++d) {
#ifdef PROP_MODEL_PUSH
						nx = x + D3Q19_X[d];
						ny = y + D3Q19_Y[d];
						nz = z + D3Q19_Z[d];
#elif PROP_MODEL_PULL
						nx = x - D3Q19_X[d];
						ny = y - D3Q19_Y[d];
						nz = z - D3Q19_Z[d];
#else
	#error PROP_MODEL_NAME unknown.
#endif
						// Check if neighbor is inside the lattice.
						// if(nx < 0 || ny < 0 || nz < 0 || nx >= lX || ny >= lY || nz >= lZ) {
						// 	continue;
						// }
						if ((nx < 0 || nx >= lX) && ld->PeriodicX) {
							++nBounceBackPdfs; // Compiler bug --> see above
						}
						else if ((ny < 0 || ny >= lY) && ld->PeriodicY) {
							++nBounceBackPdfs; // Compiler bug --> see above
						}
						else if ((nz < 0 || nz >= lZ) && ld->PeriodicZ) {
							++nBounceBackPdfs; // Compiler bug --> see above
						}
						else if (nx < 0 || ny < 0 || nz < 0 || nx >= lX || ny >= lY || nz >= lZ) {
							continue;
						}
						else if (ld->Lattice[L_INDEX_4(lDims, nx, ny, nz)] == LAT_CELL_OBSTACLE) {
							++nBounceBackPdfs; // Compiler bug --> see above
						}
					}
				}
			}
		}
	}


	printf("# allocating %d indices for bounce back pdfs (%s for source and destination array)\n", nBounceBackPdfs, ByteToHuman(sizeof(int) * nBounceBackPdfs * 2));

	MemAlloc((void **) & (kd->BounceBackPdfsSrc), sizeof(int) * nBounceBackPdfs + 100);
	MemAlloc((void **) & (kd->BounceBackPdfsDst), sizeof(int) * nBounceBackPdfs + 100);

	kd->nBounceBackPdfs = nBounceBackPdfs;
	nBounceBackPdfs = 0;

	int srcIndex;
	int dstIndex;

	for (int z = 0; z < lZ; ++z) {
		for (int y = 0; y < lY; ++y) {
			for (int x = 0; x < lX; ++x) {

				if (ld->Lattice[L_INDEX_4(ld->Dims, x, y, z)] != LAT_CELL_OBSTACLE) {
					for (int d = 0; d < N_D3Q19; ++d) {
#ifdef PROP_MODEL_PUSH
						nx = x + D3Q19_X[d];
						ny = y + D3Q19_Y[d];
						nz = z + D3Q19_Z[d];
#elif PROP_MODEL_PULL
						nx = x - D3Q19_X[d];
						ny = y - D3Q19_Y[d];
						nz = z - D3Q19_Z[d];
#else
	#error PROP_MODEL_NAME unknown.
#endif

						if ( 	((nx < 0 || nx >= lX) && ld->PeriodicX) ||
								((ny < 0 || ny >= lY) && ld->PeriodicY) ||
								((nz < 0 || nz >= lZ) && ld->PeriodicZ)
						){
							// Implement periodic boundary in X direction.

							// If the target node reached through propagation is outside the lattice
							// the kernel stores it in some buffer around the domain.
							// From this position the PDF must be transported to the other side of the
							// geometry.

							// Take PDF from outside the domain.

							// x periodic
							if (nx < 0) {
								px = lX - 1;
							}
							else if (nx >= lX) {
								px = 0;
							} else {
								px = nx;
							}

							// y periodic
							if (ny < 0) {
								py = lY - 1;
							}
							else if (ny >= lY) {
								py = 0;
							} else {
								py = ny;
							}

							// z periodic
							if (nz < 0) {
								pz = lZ - 1;
							}
							else if (nz >= lZ) {
								pz = 0;
							} else {
								pz = nz;
							}

							if (ld->Lattice[L_INDEX_4(lDims, px, py, pz)] == LAT_CELL_OBSTACLE) {
#ifdef PROP_MODEL_PUSH
								srcIndex = P_INDEX_5(gDims, nx + oX, ny + oY, nz + oZ, d);
								dstIndex = P_INDEX_5(gDims,  x + oX,  y + oY,  z + oZ, D3Q19_INV[d]);
#elif PROP_MODEL_PULL
								srcIndex = P_INDEX_5(gDims,  x + oX,  y + oY,  z + oZ, D3Q19_INV[d]);
								dstIndex = P_INDEX_5(gDims, nx + oX, ny + oY, nz + oZ, d);
#endif
							}
							else {

#ifdef PROP_MODEL_PUSH
								srcIndex = P_INDEX_5(gDims, nx + oX, ny + oY, nz + oZ, d);
								// Put it on the other side back into the domain.
								dstIndex = P_INDEX_5(gDims, px + oX, py + oY, pz + oZ, d);
#elif PROP_MODEL_PULL
								srcIndex = P_INDEX_5(gDims, px + oX, py + oY, pz + oZ, d);
								// Put it on the other side back into the ghost layer.
								dstIndex = P_INDEX_5(gDims, nx + oX, ny + oY, nz + oZ, d);
#endif

								VerifyMsg(nBounceBackPdfs < kd->nBounceBackPdfs, "nBBPdfs %d < kd->nBBPdfs %d  xyz: %d %d %d d: %d\n", nBounceBackPdfs, kd->nBounceBackPdfs, x, y, z, d);

							}

							kd->BounceBackPdfsSrc[nBounceBackPdfs] = srcIndex;
							kd->BounceBackPdfsDst[nBounceBackPdfs] = dstIndex;

							++nBounceBackPdfs;

						}
						else if (nx < 0 || ny < 0 || nz < 0 || nx >= lX || ny >= lY || nz >= lZ) {
							continue;
						}
						else if (ld->Lattice[L_INDEX_4(lDims, nx, ny, nz)] == LAT_CELL_OBSTACLE) {
#ifdef PROP_MODEL_PUSH
							srcIndex = P_INDEX_5(gDims, nx + oX, ny + oY, nz + oZ, d);
							dstIndex = P_INDEX_5(gDims,  x + oX,  y + oY,  z + oZ, D3Q19_INV[d]);
#elif PROP_MODEL_PULL
							srcIndex = P_INDEX_5(gDims,  x + oX,  y + oY,  z + oZ, D3Q19_INV[d]);
							dstIndex = P_INDEX_5(gDims, nx + oX, ny + oY, nz + oZ, d);
							// srcIndex = P_INDEX_5(gDims,  x + oX,  y + oY,  z + oZ, d);
							// dstIndex = P_INDEX_5(gDims, nx + oX, ny + oY, nz + oZ, D3Q19_INV[d]);
#endif

							VerifyMsg(nBounceBackPdfs < kd->nBounceBackPdfs, "nBBPdfs %d < kd->nBBPdfs %d  xyz: %d %d %d d: %d\n", nBounceBackPdfs, kd->nBounceBackPdfs, x, y, z, d);

							kd->BounceBackPdfsSrc[nBounceBackPdfs] = srcIndex;
							kd->BounceBackPdfsDst[nBounceBackPdfs] = dstIndex;

							++nBounceBackPdfs;
						}
					}
				}
			}
		}
	}


	// Fill remaining KernelData structures
	kd->GetNode = FNAME(GetNode);
	kd->SetNode = FNAME(SetNode);

	kd->BoundaryConditionsGetPdf = FNAME(BcGetPdf);
	kd->BoundaryConditionsSetPdf = FNAME(BcSetPdf);

	kd->Kernel = FNAME(D3Q19BlkKernel);

	kd->DstPdfs = NULL;
	kd->PdfsActive = kd->Pdfs[0];

	return;
}

void FNAME(D3Q19BlkDeinit)(LatticeDesc * ld, KernelData ** kernelData)
{
	MemFree((void **) & ((*kernelData)->Pdfs[0]));
	MemFree((void **) & ((*kernelData)->Pdfs[1]));

	MemFree((void **) & ((*kernelData)->BounceBackPdfsSrc));
	MemFree((void **) & ((*kernelData)->BounceBackPdfsDst));

	MemFree((void **)kernelData);

	return;
}

// Kernels without blocking perform the same initialization/deinitialization as with
// blocking, except that a different kernel is called. Hence, no arguments are allowed.

void FNAME(D3Q19Init)(LatticeDesc * ld, KernelData ** kernelData, Parameters * params)
{
	Parameters p;

	if (params->nKernelArgs != 0) {
		printf("ERROR: unknown kernel parameter.\n");
		printf("This kernels accepts no parameters.\n");
		exit(1);
	}

	// Setup an empty parameters structure.
	p.nArgs        = params->nArgs;
	p.Args         = params->Args;
	p.nKernelArgs  = 0;
	p.KernelArgs   = NULL;

	// Call init routine for blocking kernel and override the
	// kernel function to be called later on.
	FNAME(D3Q19BlkInit)(ld, kernelData, &p);

	(*kernelData)->Kernel = FNAME(D3Q19Kernel);

	return;

}

void FNAME(D3Q19Deinit)(LatticeDesc * ld, KernelData ** kernelData)
{
	FNAME(D3Q19BlkDeinit)(ld, kernelData);
	return;
}
