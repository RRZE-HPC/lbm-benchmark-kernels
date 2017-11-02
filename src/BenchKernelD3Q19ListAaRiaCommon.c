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
#include "BenchKernelD3Q19ListAaRiaCommon.h"

#include "Memory.h"
#include "Vtk.h"
#include "Padding.h"

#include <math.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

// Forward definition.
void FNAME(D3Q19ListAaRiaKernel)(LatticeDesc * ld, struct KernelData_ * kd, CaseData * cd);




// -----------------------------------------------------------------------
// Functions which are used as callback by the kernel to read or write
// PDFs and nodes.

static void FNAME(BCGetPdf)(KernelData * kd, int x, int y, int z, int dir, PdfT * pdf)
{
	Assert(kd != NULL);
	Assert(kd->PdfsActive != NULL);
	Assert(kd->PdfsActive == kd->Pdfs[0] || kd->PdfsActive == kd->Pdfs[1]);
	Assert(pdf != NULL);

	Assert(x >= 0); Assert(y >= 0); Assert(z >= 0);
	Assert(x < kd->Dims[0]); Assert(y < kd->Dims[1]); Assert(z < kd->Dims[2]);
	Assert(dir >= 0); Assert(dir < N_D3Q19);

	KernelDataList * kdl = (KernelDataList *)kd;

	if (kdl->Iteration % 2 == 0) {
		// Pdfs are stored inverse, local PDFs are located in remote nodes

		uint32_t nodeIndex = KDL(kd)->Grid[L_INDEX_4(kd->Dims, x, y, z)];

		if (dir != D3Q19_C) {
			uint32_t adjListIndex = nodeIndex * N_D3Q19_IDX;

			*pdf = kd->PdfsActive[KDL(kd)->AdjList[adjListIndex + D3Q19_INV[dir]]];
		}
		else {
			*pdf = kd->PdfsActive[P_INDEX_3(KDL(kd)->nCells, nodeIndex, dir)];
		}

	}
	else {
		*pdf = kd->PdfsActive[P_INDEX_5(KDL(kd), x, y, z, dir)];
	}


	return;
}

static void FNAME(BCSetPdf)(KernelData * kd, int x, int y, int z, int dir, PdfT pdf)
{
	Assert(kd != NULL);
	Assert(kd->PdfsActive != NULL);
	Assert(kd->PdfsActive == kd->Pdfs[0] || kd->PdfsActive == kd->Pdfs[1]);
	Assert(x >= 0); Assert(y >= 0); Assert(z >= 0);
	Assert(x < kd->Dims[0]); Assert(y < kd->Dims[1]); Assert(z < kd->Dims[2]);
	Assert(dir >= 0); Assert(dir < N_D3Q19);

#if 0
	if (isnan(pdf)) {
		printf("ERROR: setting nan %d %d %d %d %s\n", x, y, z, dir, D3Q19_NAMES[dir]);
		DEBUG_BREAK_POINT();
		exit(1);
	}
#endif

	KernelDataList * kdl = (KernelDataList *)kd;

	if (kdl->Iteration % 2 == 0) {
		// Pdfs are stored inverse, local PDFs are located in remote nodes

		uint32_t nodeIndex = KDL(kd)->Grid[L_INDEX_4(kd->Dims, x, y, z)];

		if (dir != D3Q19_C) {
			uint32_t adjListIndex = nodeIndex * N_D3Q19_IDX;

			kd->PdfsActive[KDL(kd)->AdjList[adjListIndex + D3Q19_INV[dir]]] = pdf;
		}
		else {
			kd->PdfsActive[P_INDEX_3(KDL(kd)->nCells, nodeIndex, dir)] = pdf;
		}

	}
	else {
		kd->PdfsActive[P_INDEX_5(KDL(kd), x, y, z, dir)] = pdf;
	}

	return;
}


static void GetNode(KernelData * kd, int x, int y, int z, PdfT * pdfs)
{
	Assert(kd != NULL);
	Assert(kd->PdfsActive != NULL);
	Assert(kd->PdfsActive == kd->Pdfs[0] || kd->PdfsActive == kd->Pdfs[1]);
	Assert(pdfs != NULL);
	Assert(x >= 0); Assert(y >= 0); Assert(z >= 0);
	Assert(x < kd->Dims[0]); Assert(y < kd->Dims[1]); Assert(z < kd->Dims[2]);

	KernelDataList * kdl = (KernelDataList *)kd;

	if(kdl->Iteration % 2 == 0){

		uint32_t nodeIndex = kdl->Grid[L_INDEX_4(kdl->kd.Dims, x, y, z)];
		uint32_t adjListIndex = nodeIndex * N_D3Q19_IDX;

		// Load PDFs of local cell: pdf_N = src[adjList[adjListIndex + D3Q19_S]]; ...
		pdfs[D3Q19_C] = kd->PdfsActive[P_INDEX_3(kdl->nCells, nodeIndex, D3Q19_C)];

		#define X(name, idx, idxinv, _x, _y, _z)	pdfs[idx] = kd->PdfsActive[kdl->AdjList[adjListIndex + idxinv]];
		D3Q19_LIST_WO_C
		#undef X

	} else {

		#define I(x, y, z, dir) P_INDEX_5(KDL(kd), (x), (y), (z), (dir))
		#define X(name, idx, idxinv, _x, _y, _z)	pdfs[idx] = kd->PdfsActive[I(x, y, z, idx)];
		D3Q19_LIST
		#undef X
		#undef I

	}

#if 0
	for (int d = 0; d < 19; ++d) {
		if(isnan(pdfs[d]) || isinf(pdfs[d])) {
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


static void SetNode(KernelData * kd, int x, int y, int z, PdfT * pdfs)
{
	Assert(kd != NULL);
	Assert(kd->PdfsActive != NULL);
	Assert(kd->PdfsActive == kd->Pdfs[0] || kd->PdfsActive == kd->Pdfs[1]);
	Assert(pdfs != NULL);

	Assert(x >= 0); Assert(y >= 0); Assert(z >= 0);
	Assert(x < kd->Dims[0]); Assert(y < kd->Dims[1]); Assert(z < kd->Dims[2]);

#if 0
	for (int d = 0; d < 19; ++d) {
		if(isnan(pdfs[d])) {
			printf("%d %d %d %d nan! get node\n", x, y, z, d);
						for (int d2 = 0; d2 < 19; ++d2) {
							printf("%d: %e\n", d2, pdfs[d2]);
						}
			exit(1);
		}
	}
#endif

	KernelDataList * kdl = (KernelDataList *)kd;

	if(kdl->Iteration % 2 == 0){

		uint32_t nodeIndex = kdl->Grid[L_INDEX_4(kdl->kd.Dims, x, y, z)];
		uint32_t adjListIndex = nodeIndex * N_D3Q19_IDX;

		// Load PDFs of local cell: pdf_N = src[adjList[adjListIndex + D3Q19_S]]; ...
		kd->PdfsActive[P_INDEX_3(kdl->nCells, nodeIndex, D3Q19_C)] = pdfs[D3Q19_C];

		#define X(name, idx, idxinv, _x, _y, _z)	kd->PdfsActive[kdl->AdjList[adjListIndex + idxinv]] = pdfs[idx];
		D3Q19_LIST_WO_C
		#undef X

	} else {

		#define I(x, y, z, dir) P_INDEX_5(KDL(kd), (x), (y), (z), (dir))
		#define X(name, idx, idxinv, _x, _y, _z)	kd->PdfsActive[I(x, y, z, idx)] = pdfs[idx];
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
#ifdef DATA_LAYOUT_SOA
	printf("  [-pad auto|modulus_1+offset_1(,modulus_n+offset_n)*]\n");
#endif

	return;
}

static void ParseParameters(Parameters * params, int * blk, PadInfo ** padInfo)
{
	Assert(blk != NULL);

	blk[0] = 0; blk[1] = 0; blk[2] = 0;
	*padInfo = NULL;

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
#ifdef DATA_LAYOUT_SOA
		else if (ARG_IS("-pad") || ARG_IS("--pad")) {
			NEXT_ARG_PRESENT();

			*padInfo = PadInfoFromStr(params->KernelArgs[++i]);
		}
#endif
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

static void SetupConsecNodes(LatticeDesc * ld, KernelDataListRia * kdlr, int nThreads)
{
	Assert(ld != NULL);
	Assert(kdlr != NULL);
	Assert(nThreads > 0);

	uint32_t * adjList = kdlr->kdl.AdjList;

	uint32_t nConsecNodes = 0;
	uint32_t consecIndex = 0;

	int nFluid = kdlr->kdl.nFluid;

	uint32_t * consecThreadIndices = (uint32_t *)malloc(sizeof(uint32_t) * (nThreads + 1));
	int * fluidNodeThreadIndices = (int *)malloc(sizeof(int) * (nThreads + 1));

	int nNodesPerThread = nFluid / nThreads;

	for (int i = 0; i < nThreads; ++i) {
		consecThreadIndices[i]     = i * nNodesPerThread + MinI(i, nFluid % nThreads);
		fluidNodeThreadIndices[i] = consecThreadIndices[i];
	}
	consecThreadIndices[nThreads]     = -1;
	fluidNodeThreadIndices[nThreads] = nFluid;

	int indexThread = 1;

	// We execute following code two times.
	// - The first time to get the count of how many entries we need for the
	//   consecNodes array.
	// - The second time to fill the array.

	// Loop over adjacency list of all nodes.
    // Compare if adjacent nodes share the same access pattern.
	for (int index = 1; index < nFluid; ++index) {

		int different = 0;

		// Loop over all directions except the center one.
		for(int d = 0; d < N_D3Q19 - 1; ++d) {
			Assert(d != D3Q19_C);

			if (adjList[index * N_D3Q19_IDX + d] != adjList[(index - 1) * N_D3Q19_IDX + d] + 1) {
				// Different access pattern.
				different = 1;
				break;
			}
		}

		if (consecThreadIndices[indexThread] == index) {
			// We are at a thread boundary. Starting from this index the fluids
			// belong to another thread. Force a break, if nodes are consecutive.
			++indexThread;
			different = 1;
		}

		if (different) {
			++consecIndex;
		}
	}

	if (nFluid > 0) {
		nConsecNodes = consecIndex + 1;
	}

	uint32_t * consecNodes;
	MemAlloc((void **)&consecNodes, sizeof(uint32_t) * nConsecNodes);

	consecIndex = 0;

	if (nFluid > 0) {
		consecNodes[consecIndex] = 1;
	}

	indexThread = 1;
	consecThreadIndices[0] = 0;

	// Loop over adjacency list of all nodes.
    // Compare if adjacent nodes share the same access pattern.
	for (int index = 1; index < nFluid; ++index) {

		int different = 0;

		// Loop over all directions except the center one.
		for(int d = 0; d < N_D3Q19 - 1; ++d) {
			Assert(d != D3Q19_C);

			if (adjList[index * N_D3Q19_IDX + d] != adjList[(index - 1) * N_D3Q19_IDX + d] + 1) {
				// Different access pattern.
				different = 1;
				break;
			}
		}

		if (consecThreadIndices[indexThread] == index) {
			// We are at a thread boundary. Starting from this index the fluids
			// belong to another thread. Force a break, if nodes are consecutive.
			consecThreadIndices[indexThread] = consecIndex + 1;
			++indexThread;
			different = 1;
		}

		if (different) {
			++consecIndex;
			Assert(consecIndex < nConsecNodes);
			consecNodes[consecIndex] = 1;
		}
		else {
			Assert(consecIndex < nConsecNodes);
			consecNodes[consecIndex] += 1;
		}
	}


	kdlr->ConsecNodes = consecNodes;
	kdlr->nConsecNodes = nConsecNodes;

	kdlr->ConsecThreadIndices  = consecThreadIndices;
	kdlr->nConsecThreadIndices = nThreads;

	kdlr->FluidNodeThreadIndices = fluidNodeThreadIndices;
	kdlr->nFluidNodeThreadIndices = nThreads;

	printf("# total fluid nodes: %d   consecutive blocks: %d\n", nFluid, nConsecNodes);

	return;
}

void FNAME(D3Q19ListAaRiaInit)(LatticeDesc * ld, KernelData ** kernelData, Parameters * params)
{
	KernelData * kd;
	KernelDataList * kdl;
	KernelDataListRia * kdlr;
	MemAlloc((void **)&kdlr, sizeof(KernelDataListRia));

	kd = (KernelData *)kdlr;
	kdl = KDL(kdlr);

	*kernelData = kd;

#ifdef DEBUG
	kd->Pdfs[0] = NULL;
	kd->Pdfs[1] = NULL;
	kd->PdfsActive = NULL;
	kd->DstPdfs = NULL;
	kd->SrcPdfs = NULL;
	kd->Dims[0] = -1;
	kd->Dims[1] = -1;
	kd->Dims[2] = -1;
	kd->GlobalDims[0] = -1;
	kd->GlobalDims[1] = -1;
	kd->GlobalDims[2] = -1;
	kd->Offsets[0] = -1;
	kd->Offsets[1] = -1;
	kd->Offsets[2] = -1;

	kd->ObstIndices = NULL;
	kd->nObstIndices = -1;
	kd->BounceBackPdfsSrc = NULL;
	kd->BounceBackPdfsDst = NULL;
	kd->nBounceBackPdfs = -1;

	kdl->AdjList = NULL;
	kdl->Coords = NULL;
	kdl->Grid = NULL;
	kdl->nCells = -1;
	kdl->nFluid = -1;

	kdlr->ConsecNodes = NULL;
	kdlr->nConsecNodes = 0;
	kdlr->ConsecThreadIndices = NULL;
	kdlr->nConsecThreadIndices = 0;
#endif

	int blk[3] = { 0 };
	PadInfo * padInfo = NULL;

	ParseParameters(params, blk, &padInfo);

	// Ajust the dimensions according to padding, if used.
	kd->Dims[0] = kd->GlobalDims[0] = ld->Dims[0];
	kd->Dims[1] = kd->GlobalDims[1] = ld->Dims[1];
	kd->Dims[2] = kd->GlobalDims[2] = ld->Dims[2];

	int * lDims = ld->Dims;

	int lX = lDims[0];
	int lY = lDims[1];
	int lZ = lDims[2];

	int nTotalCells = lX * lY * lZ;
	int nCells = ld->nFluid; // TODO: + padding
	int nFluid = ld->nFluid;

#ifdef DATA_LAYOUT_SOA
	{
		nCells = PadCellsAndReport(nCells, sizeof(PdfT), &padInfo);
		PadInfoFree(padInfo); padInfo = NULL;
	}
#endif

	kdl->nCells = nCells;
	kdl->nFluid = nFluid;

	PdfT * pdfs[2];


	if (blk[0] == 0) blk[0] = lX;
	if (blk[1] == 0) blk[1] = lY;
	if (blk[2] == 0) blk[2] = lZ;

	printf("# blocking               x: %3d y: %3d z: %3d\n", blk[0], blk[1], blk[2]);

	double latMiB      = nCells * sizeof(PdfT) * N_D3Q19 / 1024.0 / 1024.0;
	double latFluidMib = nFluid * sizeof(PdfT) * N_D3Q19 / 1024.0 / 1024.0;
	double latPadMib   = (nCells - nFluid) * sizeof(PdfT) * N_D3Q19 / 1024.0 / 1024.0;

	printf("# lattice size:          %e MiB\n", latMiB);
	printf("# fluid lattice size:    %e MiB\n", latFluidMib);
	printf("# lattice padding:       %e MiB\n", latPadMib);

#define PAGE_4K		4096

	printf("# aligning lattices to:  %d b\n", PAGE_4K);

	MemAllocAligned((void **)&pdfs[0], sizeof(PdfT) * nCells * N_D3Q19, PAGE_4K);

	kd->Pdfs[0] = pdfs[0];

	// Initialize PDFs with some (arbitrary) data for correct NUMA placement.
	// Here we touch only the fluid nodes as this loop is OpenMP parallel and
	// we want the same scheduling as in the kernel.
	#ifdef _OPENMP
		#pragma omp parallel for
	#endif
	for (int i = 0; i < nFluid; ++i) { for(int d = 0; d < N_D3Q19; ++d) {
		pdfs[0][P_INDEX_3(nCells, i, d)] = 1.0;
	} }

	// Initialize all PDFs to some standard value.
	for (int i = 0; i < nFluid; ++i) { for(int d = 0; d < N_D3Q19; ++d) {
		pdfs[0][P_INDEX_3(nCells, i, d)] = 0.0;
	} }

	// ----------------------------------------------------------------------
	// create grid which will hold the index numbers of the fluid nodes

	uint32_t * grid;

	if (MemAlloc((void **)&grid, nTotalCells * sizeof(uint32_t))) {
		printf("ERROR: allocating grid for numbering failed: %lu bytes.\n", nTotalCells * sizeof(uint32_t));
		exit(1);
	}
	kdl->Grid = grid;

	int latticeIndex;

#ifdef DEBUG
	for(int z = 0; z < lZ; ++z) {
		for(int y = 0; y < lY; ++y) {
			for(int x = 0; x < lX; ++x) {

				latticeIndex = L_INDEX_4(ld->Dims, x, y, z);

				grid[latticeIndex] = ~0;
			}
		}
	}
#endif

	// ----------------------------------------------------------------------
	// generate numbering over grid

	uint32_t * coords;

	if (MemAlloc((void **)&coords, nFluid * sizeof(uint32_t) * 3)) {
		printf("ERROR: allocating coords array failed: %lu bytes.\n", nFluid * sizeof(uint32_t) * 3);
		exit(1);
	}

	kdl->Coords = coords;

	// Index for the PDF nodes can start at 0 as we distinguish solid and fluid nodes
	// through the ld->Lattice array.
	int counter = 0;

	// Blocking is implemented via setup of the adjacency list. The kernel later will
	// walk through the lattice blocked automatically.
	for (int bX = 0; bX < lX; bX += blk[0]) {
	for (int bY = 0; bY < lY; bY += blk[1]) {
	for (int bZ = 0; bZ < lZ; bZ += blk[2]) {

		int eX = MIN(bX + blk[0], lX);
		int eY = MIN(bY + blk[1], lY);
		int eZ = MIN(bZ + blk[2], lZ);

		for (int x = bX; x < eX; ++x) {
		for (int y = bY; y < eY; ++y) {
		for (int z = bZ; z < eZ; ++z) {

			latticeIndex = L_INDEX_4(lDims, x, y, z);

			if (ld->Lattice[latticeIndex] != LAT_CELL_OBSTACLE) {
				grid[latticeIndex] = counter;

				coords[C_INDEX_X(counter)] = x;
				coords[C_INDEX_Y(counter)] = y;
				coords[C_INDEX_Z(counter)] = z;

				++counter;
			}
		} } }
	} } }

	Verify(counter == nFluid);

	uint32_t * adjList;

	double indexMib = nFluid * sizeof(uint32_t) * N_D3Q19_IDX / 1024.0 / 1024.0;

	printf("# index size:            %e MiB\n", indexMib);

	// AdjList only requires 18 instead of 19 entries per node, as
	// the center PDF needs no addressing.
	if (MemAlloc((void **)&adjList, nFluid * sizeof(uint32_t) * N_D3Q19_IDX)) {
		printf("ERROR: allocating adjList array failed: %lu bytes.\n", nFluid * sizeof(uint32_t) * N_D3Q19_IDX);
		exit(1);
	}

	kdl->AdjList = adjList;

	int x, y, z;

	uint32_t neighborIndex;
	uint32_t dstIndex;

	int nx, ny, nz, px, py, pz;

	// Loop over all fluid nodes and compute the indices to the neighboring
	// PDFs for configured data layout (AoS/SoA).
	#ifdef _OPENMP
		#pragma omp parallel for
	#endif
	for (int index = 0; index < nFluid; ++index) {
		for (int d = 0; d < N_D3Q19_IDX; ++d) {
			adjList[index * N_D3Q19_IDX + d] = -1;
		}
	}

	// #ifdef _OPENMP  --> add line continuation
	// 	#pragma omp parallel for default(none)
	// 		shared(nFluid, nCells, coords, D3Q19_INV, D3Q19_X, D3Q19_Y, D3Q19_Z,
	// 				stderr,
	// 				lDims, grid, ld, lX, lY, lZ, adjList)
	// 		private(x, y, z, nx, ny, nz, neighborIndex, dstIndex)
	// #endif
	for (int index = 0; index < nFluid; ++index) {
		x = coords[C_INDEX_X(index)];
		y = coords[C_INDEX_Y(index)];
		z = coords[C_INDEX_Z(index)];

		Assert(x >= 0 && x < lX);
		Assert(y >= 0 && y < lY);
		Assert(z >= 0 && z < lZ);

		Assert(ld->Lattice[L_INDEX_4(lDims, x, y, z)] != LAT_CELL_OBSTACLE);

		// Loop over all directions except the center one.
		for(int d = 0; d < N_D3Q19 - 1; ++d) {
			Assert(d != D3Q19_C);

#ifdef PROP_MODEL_PUSH
			nx = x + D3Q19_X[d];
			ny = y + D3Q19_Y[d];
			nz = z + D3Q19_Z[d];

#elif PROP_MODEL_PULL
			nx = x - D3Q19_X[d];
			ny = y - D3Q19_Y[d];
			nz = z - D3Q19_Z[d];
#else
			#error No implementation for this PROP_MODEL_NAME.
#endif
			// If the neighbor is outside the latcie in X direction and we have a
			// periodic boundary then we need to wrap around.
			if ( 	((nx < 0 || nx >= lX) && ld->PeriodicX) ||
					((ny < 0 || ny >= lY) && ld->PeriodicY) ||
					((nz < 0 || nz >= lZ) && ld->PeriodicZ)
																){
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
					dstIndex = P_INDEX_3(nCells, index, D3Q19_INV[d]);
				}
				else {
					neighborIndex = grid[L_INDEX_4(lDims, px, py, pz)];

					AssertMsg(neighborIndex != ~0, "Neighbor has no Index. (%d %d %d) direction %s (%d)\n", px, py, pz, D3Q19_NAMES[d], d);

					dstIndex = P_INDEX_3(nCells, neighborIndex, d);
				}
			}
			else if (nx < 0 || ny < 0 || nz < 0 || nx >= lX || ny >= lY || nz >= lZ) {
				dstIndex = P_INDEX_3(nCells, index, D3Q19_INV[d]);
			}
			else if (ld->Lattice[L_INDEX_4(lDims, nx, ny, nz)] == LAT_CELL_OBSTACLE) {
				dstIndex = P_INDEX_3(nCells, index, D3Q19_INV[d]);
			}
			else {
				neighborIndex = grid[L_INDEX_4(lDims, nx, ny, nz)];

				Assert(neighborIndex != ~0);

				dstIndex = P_INDEX_3(nCells, neighborIndex, d);
			}

			Assert(dstIndex >= 0);
			Assert(dstIndex < nCells * N_D3Q19);

			adjList[index * N_D3Q19_IDX + d] = dstIndex;
		}
	}

	int nThreads = 1;

#ifdef _OPENMP
	nThreads = omp_get_max_threads();
#endif

	SetupConsecNodes(ld, KDLR(kd), nThreads);

	double loopBalanceEven = 2.0 * 19 * sizeof(PdfT);
	double loopBalanceOdd  = 2.0 * 19 * sizeof(PdfT) + (double)kdlr->nConsecNodes / nFluid * (18 * 4.0 + 4.0);
	double loopBalance     = (loopBalanceEven + loopBalanceOdd) / 2.0;

	printf("# loop balance:          %.2f B/FLUP  even: %.2f B/FLUP  odd %.2f B/FLUP\n",
			loopBalance, loopBalanceEven, loopBalanceOdd);

	// Fill remaining KernelData structures
	kd->GetNode = GetNode;
	kd->SetNode = SetNode;

	kd->BoundaryConditionsGetPdf = FNAME(BCGetPdf);
	kd->BoundaryConditionsSetPdf = FNAME(BCSetPdf);

	kd->Kernel = FNAME(D3Q19ListAaRiaKernel);

	kd->DstPdfs = NULL;
	kd->PdfsActive = kd->Pdfs[0];

	return;
}

void FNAME(D3Q19ListAaRiaDeinit)(LatticeDesc * ld, KernelData ** kernelData)
{
	KernelDataListRia ** kdlr = (KernelDataListRia **)kernelData;

	MemFree((void **)&((*kdlr)->ConsecNodes));

	if ((*kdlr)->ConsecThreadIndices != NULL) {
		MemFree((void **)&((*kdlr)->ConsecThreadIndices));
	}

	if ((*kdlr)->FluidNodeThreadIndices != NULL) {
		MemFree((void **)&((*kdlr)->FluidNodeThreadIndices));
	}

	KernelDataList ** kdl = (KernelDataList **)kernelData;

	MemFree((void **)&((*kdl)->AdjList));
	MemFree((void **)&((*kdl)->Coords));
	MemFree((void **)&((*kdl)->Grid));

	MemFree((void **)&((*kernelData)->Pdfs[0]));

	MemFree((void **)kernelData);
	return;
}

