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
#include "BenchKernelD3Q19AaVecCommon.h"

#include "Memory.h"
#include "Vtk.h"
#include "Vector.h"

#include <inttypes.h>
#include <math.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

// Forward definition.
void FNAME(D3Q19AaVecKernel)(LatticeDesc * ld, struct KernelData_ * kd, CaseData * cd);


static void FNAME(BcGetPdf)(KernelData * kd, int x, int y, int z, int dir, PdfT * pdf)
{
	Assert(kd != NULL);
	Assert(kd->PdfsActive != NULL);
	Assert(kd->PdfsActive == kd->Pdfs[0] || kd->PdfsActive == kd->Pdfs[1]);
	Assert(pdf != NULL);

	Assert(x >= 0); Assert(y >= 0); Assert(z >= 0);
	Assert(x < kd->Dims[0]); Assert(y < kd->Dims[1]); Assert(z < kd->Dims[2]);
	Assert(dir >= 0); Assert(dir < N_D3Q19);

	KernelDataAa * kda = KDA(kd);

	int oX = kd->Offsets[0];
	int oY = kd->Offsets[1];
	int oZ = kd->Offsets[2];

	if (kda->Iteration % 2 == 0) {
		// Pdfs are stored inverse, local PDFs are located in remote nodes
		int nx = x - D3Q19_X[dir];
		int ny = y - D3Q19_Y[dir];
		int nz = z - D3Q19_Z[dir];

		#define I(x, y, z, dir) P_INDEX_5(kd->GlobalDims, (x), (y), (z), (dir))
		*pdf = kd->PdfsActive[I(nx + oX, ny + oY, nz + oZ, D3Q19_INV[dir])];
		#undef I
	}
	else {
		int nx = x;
		int ny = y;
		int nz = z;

		#define I(x, y, z, dir) P_INDEX_5(kd->GlobalDims, (x), (y), (z), (dir))
		*pdf = kd->PdfsActive[I(nx + oX, ny + oY, nz + oZ, dir)];
		#undef I
	}


	return;
}

static void FNAME(BcSetPdf)(KernelData * kd, int x, int y, int z, int dir, PdfT pdf)
{
	Assert(kd != NULL);
	Assert(kd->PdfsActive != NULL);
	Assert(kd->PdfsActive == kd->Pdfs[0] || kd->PdfsActive == kd->Pdfs[1]);

	Assert(x >= 0); Assert(y >= 0); Assert(z >= 0);
	Assert(x < kd->Dims[0]); Assert(y < kd->Dims[1]); Assert(z < kd->Dims[2]);
	Assert(dir >= 0); Assert(dir < N_D3Q19);

	KernelDataAa * kda = KDA(kd);

	int oX = kd->Offsets[0];
	int oY = kd->Offsets[1];
	int oZ = kd->Offsets[2];

	if (kda->Iteration % 2 == 0) {
		// Pdfs are stored inverse, local PDFs are located in remote nodes
		int nx = x - D3Q19_X[dir];
		int ny = y - D3Q19_Y[dir];
		int nz = z - D3Q19_Z[dir];

		#define I(x, y, z, dir) P_INDEX_5(kd->GlobalDims, (x), (y), (z), (dir))
		pdf = kd->PdfsActive[I(nx + oX, ny + oY, nz + oZ, D3Q19_INV[dir])] = pdf;
		#undef I
	}
	else {
		int nx = x;
		int ny = y;
		int nz = z;

		#define I(x, y, z, dir) P_INDEX_5(kd->GlobalDims, (x), (y), (z), (dir))
		kd->PdfsActive[I(nx + oX, ny + oY, nz + oZ, dir)] = pdf;
		#undef I
	}

	return;
}


static void FNAME(GetNode)(KernelData * kd, int x, int y, int z, PdfT * pdfs)
{
	Assert(kd != NULL);
	Assert(kd->PdfsActive != NULL);
	Assert(kd->PdfsActive == kd->Pdfs[0] || kd->PdfsActive == kd->Pdfs[1]);
	Assert(pdfs != NULL);

	Assert(x >= 0); Assert(y >= 0); Assert(z >= 0);
	Assert(x < kd->Dims[0]); Assert(y < kd->Dims[1]); Assert(z < kd->Dims[2]);

	KernelDataAa * kda = KDA(kd);

	int oX = kd->Offsets[0];
	int oY = kd->Offsets[1];
	int oZ = kd->Offsets[2];


	if (kda->Iteration % 2 == 0) {
		// Pdfs are stored inverse, local PDFs are located in remote nodes

		#define I(x, y, z, dir) P_INDEX_5(kd->GlobalDims, (x), (y), (z), (dir))
		#define X(name, idx, idxinv, _x, _y, _z)	pdfs[idx] = kd->PdfsActive[I(x + oX - _x, y + oY - _y, z + oZ - _z, D3Q19_INV[idx])];
		D3Q19_LIST
		#undef X
		#undef I
	}
	else {
		#define I(x, y, z, dir) P_INDEX_5(kd->GlobalDims, (x), (y), (z), (dir))
		#define X(name, idx, idxinv, _x, _y, _z)	pdfs[idx] = kd->PdfsActive[I(x + oX, y + oY, z + oZ, idx)];
		D3Q19_LIST
		#undef X
		#undef I

	}

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

	Assert(x >= 0); Assert(y >= 0); Assert(z >= 0);
	Assert(x < kd->Dims[0]); Assert(y < kd->Dims[1]); Assert(z < kd->Dims[2]);

	KernelDataAa * kda = KDA(kd);

	int oX = kd->Offsets[0];
	int oY = kd->Offsets[1];
	int oZ = kd->Offsets[2];

	if (kda->Iteration % 2 == 0) {
		// Pdfs are stored inverse, local PDFs are located in remote nodes

		#define I(x, y, z, dir) P_INDEX_5(kd->GlobalDims, (x), (y), (z), (dir))
		#define X(name, idx, idxinv, _x, _y, _z)	kd->PdfsActive[I(x + oX - _x, y + oY - _y, z + oZ - _z, D3Q19_INV[idx])] = pdfs[idx];
		D3Q19_LIST
		#undef X
		#undef I
	}
	else {
		#define I(x, y, z, dir) P_INDEX_5(kd->GlobalDims, (x), (y), (z), (dir))
		#define X(name, idx, idxinv, _x, _y, _z)	kd->PdfsActive[I(x + oX, y + oY, z + oZ, idx)] = pdfs[idx];
		D3Q19_LIST
		#undef X
		#undef I
	}
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

			if (tmp < 0) {
				printf("ERROR: blocking parameter must be >= 0.\n");
				exit(1);
			}

			blk[0] = blk[1] = blk[2] = tmp;
		}
		else if (ARG_IS("-blk-x") || ARG_IS("--blk-x")) {
			NEXT_ARG_PRESENT();

			int tmp = strtol(params->KernelArgs[++i], NULL, 0);

			if (tmp < 0) {
				printf("ERROR: blocking parameter must be >= 0.\n");
				exit(1);
			}

			blk[0] = tmp;
		}
		else if (ARG_IS("-blk-y") || ARG_IS("--blk-y")) {
			NEXT_ARG_PRESENT();

			int tmp = strtol(params->KernelArgs[++i], NULL, 0);

			if (tmp < 0) {
				printf("ERROR: blocking parameter must be >= 0.\n");
				exit(1);
			}

			blk[1] = tmp;
		}
		else if (ARG_IS("-blk-z") || ARG_IS("--blk-z")) {
			NEXT_ARG_PRESENT();

			int tmp = strtol(params->KernelArgs[++i], NULL, 0);

			if (tmp < 0) {
				printf("ERROR: blocking parameter must be >= 0.\n");
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


void FNAME(D3Q19AaVecInit)(LatticeDesc * ld, KernelData ** kernelData, Parameters * params)
{
	KernelDataAa * kda = NULL;
	MemAlloc((void **)&kda, sizeof(KernelDataAa));

	kda->Blk[0] = 0; kda->Blk[1] = 0; kda->Blk[2] = 0;
	kda->Iteration = -1;

	KernelData * kd = &kda->kd;
	*kernelData = kd;

	kd->nObstIndices = ld->nObst;

	// Ajust the dimensions according to padding, if used.
	kd->Dims[0] = ld->Dims[0];
	kd->Dims[1] = ld->Dims[1];
	kd->Dims[2] = ld->Dims[2];


	int * lDims = ld->Dims;
	int * gDims = kd->GlobalDims;

	Assert(VSIZE <= 4);

	// TODO: only add enough ghost cells so we can compute everything vectorized.
	gDims[0] = lDims[0] + 2;
	gDims[1] = lDims[1] + 2;
	// TODO: fix this for aa-vec2-soa
	gDims[2] = lDims[2] + 2 + VSIZE - 2; // one ghost cell in front, one in the back, plus at most two at the back for VSIZE = 4

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

	PdfT * pdfs[2] = { NULL, NULL };

	ParseParameters(params, blk);

	if (blk[2] % VSIZE != 0) {
		blk[2] -= blk[2] % VSIZE;
		printf("WARNING: blocking factor for z direction must be a multiple of VSIZE = %d, adjusting it to %d.\n", VSIZE, blk[2]);
	}

	if (blk[0] == 0) blk[0] = gX;
	if (blk[1] == 0) blk[1] = gY;
	if (blk[2] == 0) blk[2] = gZ;

	printf("# blocking x: %3d y: %3d z: %3d\n", blk[0], blk[1], blk[2]);

	kda->Blk[0] = blk[0]; kda->Blk[1] = blk[1]; kda->Blk[2] = blk[2];


	printf("# allocating data for %d LB nodes with padding (%lu bytes = %f MiB for the single lattice)\n",
	       nCells,
		   sizeof(PdfT) * nCells * N_D3Q19,
	       sizeof(PdfT) * nCells * N_D3Q19 / 1024.0 / 1024.0);

#define PAGE_4K		4096

	MemAllocAligned((void **)&pdfs[0], sizeof(PdfT) * nCells * N_D3Q19, PAGE_4K);

	kd->Pdfs[0] = pdfs[0];
	kd->Pdfs[1] = NULL;


	// Initialize PDFs with some (arbitrary) data for correct NUMA placement.
	// This depends on the chosen data layout.
	// The structure of the loop should resemble the same "execution layout"
	// as in the kernel!

	int nThreads;
#ifdef _OPENMP
	nThreads = omp_get_max_threads();
#endif

#ifdef _OPENMP
	#pragma omp parallel for \
				shared(gDims, pdfs, \
				oX, oY, oZ, lX, lY, lZ, blk, nThreads, ld)
#endif
	for (int i = 0; i < nThreads; ++i) {

		int threadStartX = lX / nThreads * i;
		int threadEndX   = lX / nThreads * (i + 1);

		if (lX % nThreads > 0) {
			if (lX % nThreads > i) {
				threadStartX += i;
				threadEndX   += i + 1;
			}
			else {
				threadStartX += lX % nThreads;
				threadEndX   += lX % nThreads;
			}
		}

		for (int bX = oX + threadStartX; bX < threadEndX + oX; bX += blk[0]) {
		for (int bY = oY; bY < lY + oY; bY += blk[1]) {
		for (int bZ = oZ; bZ < lZ + oZ; bZ += blk[2]) {

			int eX = MIN(bX + blk[0], threadEndX + oX);
			int eY = MIN(bY + blk[1], lY + oY);
			int eZ = MIN(bZ + blk[2], lZ + oZ);

			// printf("%d: %d-%d  %d-%d  %d-%d  %d - %d\n", omp_get_thread_num(), bZ, eZ, bY, eY, bX, eX, threadStartX, threadEndX);

			for (int x = bX; x < eX; ++x) {
			for (int y = bY; y < eY; ++y) {
			for (int z = bZ; z < eZ; ++z) {

				for (int d = 0; d < N_D3Q19; ++d) {
					pdfs[0][P_INDEX_5(gDims, x, y, z, d)] = -50.0;
				}

			} } } // x, y, z
		} } } // bx, by, bz
	}


	// Initialize all PDFs to some standard value.
	for (int x = oX; x < lX + oX; ++x) {
	for (int y = oY; y < lY + oY; ++y) {
	for (int z = oZ; z < lZ + oZ; ++z) {
		for (int d = 0; d < N_D3Q19; ++d) {
			pdfs[0][P_INDEX_5(gDims, x, y, z, d)] = 0.0;
		}
	} } } // x, y, z


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


	for (int x = 0; x < lX; ++x) {
		for (int y = 0; y < lY; ++y) {
			for (int z = 0; z < lZ; ++z) {

				if (ld->Lattice[L_INDEX_4(ld->Dims, x, y, z)] != LAT_CELL_OBSTACLE) {
					for (int d = 0; d < N_D3Q19; ++d) {
						nx = x - D3Q19_X[d];
						ny = y - D3Q19_Y[d];
						nz = z - D3Q19_Z[d];

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

	// TODO: currently this is not NUMA-aware
    //       - maybe use the same blocking as for lattice initialization?
	//       - do place the bounce back index vector parallel?

	for (int x = 0; x < lX; ++x) {
		for (int y = 0; y < lY; ++y) {
			for (int z = 0; z < lZ; ++z) {

				if (ld->Lattice[L_INDEX_4(ld->Dims, x, y, z)] != LAT_CELL_OBSTACLE) {
					for (int d = 0; d < N_D3Q19; ++d) {
						nx = x + D3Q19_X[d];
						ny = y + D3Q19_Y[d];
						nz = z + D3Q19_Z[d];

						if ( 	((nx < 0 || nx >= lX) && ld->PeriodicX) ||
								((ny < 0 || ny >= lY) && ld->PeriodicY) ||
								((nz < 0 || nz >= lZ) && ld->PeriodicZ)
						){
							// For periodicity:

							// We assume we have finished odd time step (accessing neighbor PDFs) and are
							// before executing the even time step (accessing local PDFs only).

							// Assuming we are at the most east position of the lattice. Through the odd
							// time step propagation has put a PDF in the east slot of the ghost cell east
							// of us, i.e. nx, ny, nz.  We copy it to the east slot of the most west node.

							// In case of transition from even to odd time step , src and dest must be
							// exchanged.


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
								// See description of bounce back handling below.
								srcIndex = P_INDEX_5(gDims, nx + oX, ny + oY, nz + oZ, d);
								dstIndex = P_INDEX_5(gDims,  x + oX,  y + oY,  z + oZ, D3Q19_INV[d]);
							}
							else {

								srcIndex = P_INDEX_5(gDims, nx + oX, ny + oY, nz + oZ, d);
								// Put it on the other side back into the domain.
								dstIndex = P_INDEX_5(gDims, px + oX, py + oY, pz + oZ, d);

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
							// Depending on the time step we are in we have to exchange src and dst index.

							// We build the list for the case, when we have finished odd time step
							// (accessing neighbor PDFs) and before we start with the even time step
							// (accessing local PDFs only).

							// Assume our neighbor east of us, i.e. nx, ny, nz, is an obstacle cell.
							// Then we have to move the east PDF from the obstacle to our west position,
							// i.e. the inverse of east.

							// In case of transition from even to odd time step src and dest just
							// have to be exchanged.

							srcIndex = P_INDEX_5(gDims, nx + oX, ny + oY, nz + oZ, d);
							dstIndex = P_INDEX_5(gDims,  x + oX,  y + oY,  z + oZ, D3Q19_INV[d]);

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

	kd->Kernel = FNAME(D3Q19AaVecKernel);

	kd->DstPdfs = NULL;
	kd->PdfsActive = kd->Pdfs[0];

	return;
}

void FNAME(D3Q19AaVecDeinit)(LatticeDesc * ld, KernelData ** kernelData)
{
	MemFree((void **) & ((*kernelData)->Pdfs[0]));

	MemFree((void **) & ((*kernelData)->BounceBackPdfsSrc));
	MemFree((void **) & ((*kernelData)->BounceBackPdfsDst));

	MemFree((void **)kernelData);

	return;
}

