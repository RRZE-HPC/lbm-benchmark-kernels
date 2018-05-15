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
//   Michael Hussnaetter, 2017-2018
//   University of Erlangen-Nuremberg, Germany
//   michael.hussnaetter -at- fau.de
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
#include "BenchKernelD3Q19ListAaPvGatherCommon.h"

#include "Memory.h"
#include "Vector.h"
#include "Vtk.h"

#include <math.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#define PAGE_4K		4096

#if ALLOC_ADJ_LIST_IN_HBM == 1
	#define ADJ_LIST_ALLOCATOR HbwAllocAligned
	#define ADJ_LIST_FREE HbwFree
#else
	#define ADJ_LIST_ALLOCATOR MemAllocAligned
	#define ADJ_LIST_FREE MemFree
#endif

#if ALLOC_PDF_IN_HBM == 1
	#define PDF_ALLOCATOR HbwAllocAligned
	#define PDF_FREE HbwFree
#else
	#define PDF_ALLOCATOR MemAllocAligned
	#define PDF_FREE MemFree
#endif

// Forward definition.
void FNAME(D3Q19ListAaPvGatherKernel)(LatticeDesc * ld, struct KernelData_ * kd, CaseData * cd);


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

	KernelDataList * kdl = KDL(kd);
	uint32_t * adjList = kdl->AdjList;

	if (kdl->Iteration % 2 == 0) {
		// Pdfs are stored inverse, local PDFs are located in remote nodes

		// getting node index
		uint32_t index = kdl->Grid[L_INDEX_4(kd->Dims, x, y, z)];

		if (dir != D3Q19_C) {
			#define ADJ_LIST(dir) adjList[(index - (index % VSIZE)) * N_D3Q19_IDX + (dir * VSIZE) + (index % VSIZE)]
			*pdf = kd->PdfsActive[ADJ_LIST(D3Q19_INV[dir])];
			#undef ADJ_LIST
		}
		else {
			*pdf = kd->PdfsActive[P_INDEX_3(kdl->nCells, index, dir)];
		}

	}
	else {
		*pdf = kd->PdfsActive[P_INDEX_5(kdl, x, y, z, dir)];
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

	if (isnan(pdf)) {
		printf("ERROR: setting nan %d %d %d %d %s\n", x, y, z, dir, D3Q19_NAMES[dir]);
		DEBUG_BREAK_POINT();
		exit(1);
	}

	KernelDataList * kdl = KDL(kd);
	uint32_t * adjList = kdl->AdjList;

	if (kdl->Iteration % 2 == 0) {
		// Pdfs are stored inverse, local PDFs are located in remote nodes

		// getting node index
		uint32_t index = kdl->Grid[L_INDEX_4(kd->Dims, x, y, z)];

		if (dir != D3Q19_C) {
			#define ADJ_LIST(dir) adjList[(index - (index % VSIZE)) * N_D3Q19_IDX + (dir * VSIZE) + (index % VSIZE)]
			kd->PdfsActive[ADJ_LIST(D3Q19_INV[dir])] = pdf;
			#undef ADJ_LIST
		} else {
			kd->PdfsActive[P_INDEX_3(kdl->nCells, index, dir)] = pdf;
		}

	} else {
		kd->PdfsActive[P_INDEX_5(kdl, x, y, z, dir)] = pdf;
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

	KernelDataList * kdl = KDL(kd);
	uint32_t * adjList = kdl->AdjList;

	if(kdl->Iteration % 2 == 0){

		uint32_t index = kdl->Grid[L_INDEX_4(kdl->kd.Dims, x, y, z)];

		// Load PDFs of local cell: pdf_N = src[adjList[adjListIndex + D3Q19_S]]; ...
		#define ADJ_LIST(dir) adjList[(index - (index % VSIZE)) * N_D3Q19_IDX + (dir * VSIZE) + (index % VSIZE)]
		#define X(name, idx, idxinv, _x, _y, _z)	pdfs[idx] = kd->PdfsActive[ADJ_LIST(idxinv)];
			D3Q19_LIST_WO_C
		#undef X
		#undef ADJ_LIST
		pdfs[D3Q19_C] = kd->PdfsActive[P_INDEX_3(kdl->nCells, index, D3Q19_C)];

	} else {

		#define I(x, y, z, dir) P_INDEX_5(KDL(kd), (x), (y), (z), (dir))
		#define X(name, idx, idxinv, _x, _y, _z)	pdfs[idx] = kd->PdfsActive[I(x, y, z, idx)];
			D3Q19_LIST
		#undef X
		#undef I

	}

	for (int d = 0; d < 19; ++d) {
		if(isnan(pdfs[d]) || isinf(pdfs[d])) {
			printf("%d %d %d %d nan! get node\n", x, y, z, d);
			for (int d2 = 0; d2 < 19; ++d2) {
				printf("%d: %e\n", d2, pdfs[d2]);
			}
			exit(1);
		}
	}

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

	for (int d = 0; d < 19; ++d) {
		if(isnan(pdfs[d])) {
			printf("%d %d %d %d nan! get node\n", x, y, z, d);
			for (int d2 = 0; d2 < 19; ++d2) {
				printf("%d: %e\n", d2, pdfs[d2]);
			}
			exit(1);
		}
	}

	KernelDataList * kdl = KDL(kd);
	uint32_t * adjList = kdl->AdjList;

	if(kdl->Iteration % 2 == 0){

		uint32_t index = kdl->Grid[L_INDEX_4(kdl->kd.Dims, x, y, z)];

		// Load PDFs of local cell: pdf_N = src[adjList[adjListIndex + D3Q19_S]]; ...
		kd->PdfsActive[P_INDEX_3(kdl->nCells, index, D3Q19_C)] = pdfs[D3Q19_C];

		#define ADJ_LIST(dir) adjList[(index - (index % VSIZE)) * N_D3Q19_IDX + (dir * VSIZE) + (index % VSIZE)]
		#define X(name, idx, idxinv, _x, _y, _z)	kd->PdfsActive[ADJ_LIST(idxinv)] = pdfs[idx];
			D3Q19_LIST_WO_C
		#undef X
		#undef ADJ_LIST

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

	int nNodesPerThread = nFluid / nThreads;

	for (int i = 0; i < nThreads; ++i) {
		consecThreadIndices[i] = i * nNodesPerThread + MinI(i, nFluid % nThreads);
	}
	consecThreadIndices[nThreads] = nFluid;

	int indexThread = 1;
	int similarPatterns = 0;
	int wasLastChunkThreadBoundary = 0;
	// We execute following code two times.
	// - The first time to get the count of how many entries we need for the
	//   consecNodes array.
	// - The second time to fill the array.

	// Loop over adjacency list of all nodes.
    // Compare if adjacent nodes share the same access pattern.
	for (int fluidBaseIndex = VSIZE; fluidBaseIndex < nFluid; fluidBaseIndex += VSIZE) {

		int hasSimilarAccessPattern = 1;

		// Loop over all directions except the center one.
		for(int d = 0; d < N_D3Q19 - 1; ++d) {
			Assert(d != D3Q19_C);

			// check if cache line itself has consecutive memory access pattern
			for(int inChunkIndex = 0; (inChunkIndex < VSIZE) && ((fluidBaseIndex + inChunkIndex) < nFluid); ++inChunkIndex){
				int index = fluidBaseIndex + inChunkIndex;

				Assert(index < nFluid);

				#define ADJ_LIST(idx, dir) adjList[((idx) - ((idx) % VSIZE)) * N_D3Q19_IDX + ((dir) * VSIZE) + ((idx) % VSIZE)]
				//if (adjList[index * N_D3Q19_IDX + d] != adjList[(index - 1) * N_D3Q19_IDX + d] + 1)
				if (ADJ_LIST(index, d) != ADJ_LIST(index-VSIZE, d) + VSIZE) {
					//printf("different @: ADJ_LST(%d,%d)=%d != %d=ADJ_LST(%d, %d) + VSIZE\n", index, d, ADJ_LIST(index,d), ADJ_LIST(index-VSIZE,d) + VSIZE, index-VSIZE, d);
					// Different access pattern.
					hasSimilarAccessPattern = 0;
					break;
				}
				#undef ADJ_LIST
			}

			if(!hasSimilarAccessPattern){
				break; //exit from nested loop
			}
		}

		long threadBoundaryIndex = consecThreadIndices[indexThread];
		if (fluidBaseIndex <= threadBoundaryIndex &&
				threadBoundaryIndex < fluidBaseIndex + VSIZE) {
			// Current chunk contains thread boundary.
			// These chunks are treated by scalar peel and reminder loops
			// in kernel of every thread to ensure VSIZE aligned access to
			// adjacency list

			// final cells of current thread
			++consecIndex;

			// first cells of next thread
			++indexThread;
			++consecIndex;

			wasLastChunkThreadBoundary = 1;
		}
		else {
			// We are not at a thread boundary
			if (hasSimilarAccessPattern && !wasLastChunkThreadBoundary){
				++similarPatterns;
			}
			else {
				++consecIndex;
			}

			wasLastChunkThreadBoundary = 0;

			/*
			if (!hasSimilarAccessPattern) {
				++consecIndex;
			}
			else {
				++similarPatterns;
			}
			*/
		}
	}

	if (nFluid > 0) {
		nConsecNodes = consecIndex + 1;
	}

	uint32_t * consecNodes;
	MemAlloc((void **)&consecNodes, sizeof(uint32_t) * nConsecNodes);

	unsigned long consecNodesByte = (nConsecNodes) * sizeof(uint32_t);

	printf("# Consec. Nodes Array Allocation:\n");
	printf("#   similar patterns\t\t%d\n", similarPatterns);
	printf("#   elements:      \t\t%d\n",   nConsecNodes);
	printf("#   size:          \t\t%e MiB\n", consecNodesByte / 1024.0 / 1024.0);
	printf("#   alignment:     \t\t%d b\n",   PAGE_4K);

	if (MemAllocAligned((void **)&consecNodesByte, consecNodesByte, PAGE_4K)) {
			printf("ERROR: allocating consecNodes array with MemAllocAligned failed: %lu bytes.\n", consecNodesByte);
			exit(1);
	}
	else {
		printf("#   allocator: \t\t\tMemAllocAligned()\n");
	}

	consecIndex = 0;

	if (nFluid > 0) {
		consecNodes[consecIndex] = VSIZE;
	}

	indexThread = 1;
	consecThreadIndices[0] = 0;

	//add first chunk manually to enable backward check for consecutive pattern
	consecNodes[consecIndex] = VSIZE;

	wasLastChunkThreadBoundary = 0;

	// Loop over adjacency list of all nodes.
    // Compare if access pattern does not change on chunk level
	// Since gather instructions are used, access pattern may not be consecutive
	for (int fluidBaseIndex = VSIZE; fluidBaseIndex < nFluid; fluidBaseIndex += VSIZE) {

		int hasSimilarAccessPattern = 1;

		// Loop over all directions except the center one.
		for(int d = 0; d < N_D3Q19 - 1; ++d) {
			Assert(d != D3Q19_C);

			// check if cache line itself has consecutive memory access pattern
			for(int inChunkIndex = 0; (inChunkIndex < VSIZE) && ((fluidBaseIndex + inChunkIndex) < nFluid); ++inChunkIndex){
				int index = fluidBaseIndex + inChunkIndex;

				Assert(index < nFluid);

				#define ADJ_LIST(idx, dir) adjList[((idx) - ((idx) % VSIZE)) * N_D3Q19_IDX + ((dir) * VSIZE) + ((idx) % VSIZE)]
				//if (adjList[index * N_D3Q19_IDX + d] != adjList[(index - 1) * N_D3Q19_IDX + d] + 1)
				if (ADJ_LIST(index, d) != ADJ_LIST(index-VSIZE, d) + VSIZE) {
					// Different access pattern.
					hasSimilarAccessPattern = 0;
					break;
				}
				#undef ADJ_LIST
			}

			if(!hasSimilarAccessPattern){
				break; //exit from nested loop
			}
		}

		long threadBoundaryIndex = consecThreadIndices[indexThread];
		if (fluidBaseIndex <= threadBoundaryIndex &&
				threadBoundaryIndex < fluidBaseIndex + VSIZE) {
			// Current chunk contains thread boundary.
			// These chunks are treated by scalar peel and reminder loops
			// in kernel of every thread to ensure VSIZE aligned access to
			// adjacency list

			// final cells of current thread
			++consecIndex;
			//consecThreadIndices[indexThread] = consecIndex;
			consecNodes[consecIndex] = threadBoundaryIndex - fluidBaseIndex;


			// first cells of next thread
			++consecIndex;
			consecThreadIndices[indexThread] = consecIndex;
			consecNodes[consecIndex] = (fluidBaseIndex + VSIZE) - threadBoundaryIndex;
			++indexThread;

			wasLastChunkThreadBoundary = 1;

		}
		else {
			// We are not at a thread boundary
			if (hasSimilarAccessPattern && !wasLastChunkThreadBoundary){
				Assert(consecIndex < nConsecNodes);
				consecNodes[consecIndex] += VSIZE;
			}
			else {
				++consecIndex;
				Assert(consecIndex < nConsecNodes);
				consecNodes[consecIndex] = VSIZE;
			}

			/*
			if (!hasSimilarAccessPattern) {
				++consecIndex;
				Assert(consecIndex < nConsecNodes);
				consecNodes[consecIndex] = VSIZE;
			}
			else {
				Assert(consecIndex < nConsecNodes);
				consecNodes[consecIndex] += VSIZE;
			}
			*/
			wasLastChunkThreadBoundary = 0;

		}
	}

	/*
	printf("consecNodes:\n");
	for(int i = 0; i < nConsecNodes + 5; ++i){
		printf("%d ", consecNodes[i]);
	}
	printf("\n");
	*/
	/*
	printf("consecThreadIndices:\n");
	for(int i = 0; i < nThreads + 5; ++i){
		printf("%d ", consecThreadIndices[i]);
	}
	printf("\n");
	*/

	kdlr->ConsecNodes = consecNodes;
	kdlr->nConsecNodes = nConsecNodes;

	kdlr->ConsecThreadIndices  = consecThreadIndices;
	kdlr->nConsecThreadIndices = nThreads;

	double loopBalanceEven = 2.0 * 19 * sizeof(PdfT);
	//N_D3Q19 - 1: C lookup not required, +1: transfer of consecValue
	double loopBalanceOdd  = 2.0 * 19 * sizeof(PdfT) + ((double)nConsecNodes *((N_D3Q19 - 1) * VSIZE + 1)) / nFluid * sizeof(int);
	double loopBalance     = (loopBalanceEven + loopBalanceOdd) / 2.0;

	kdlr->kdl.kd.LoopBalance = loopBalance;

	printf("# loop balance:\n");
	printf("#   even timestep:  \t\t%.2f B/FLUP\n", loopBalanceEven);
	printf("#   odd timestep:   \t\t%.2f B/FLUP\n", loopBalanceOdd);
	printf("#   average:        \t\t%.2f B/FLUP\n", loopBalance);

	return;
}

void FNAME(D3Q19ListAaPvGatherInit)(LatticeDesc * ld, KernelData ** kernelData, Parameters * params)
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

	// TODO: check nCells/nFluid do not exceed 2^31. This actually has to be
	// done during lattice setup.
	kdl->nCells = nCells;
	kdl->nFluid = nFluid;

	PdfT * pdfs[2];

	int blk[3] = { 0 };

	ParseParameters(params, blk);

	if (blk[0] == 0) blk[0] = lX;
	if (blk[1] == 0) blk[1] = lY;
	if (blk[2] == 0) blk[2] = lZ;

	printf("# blocking:            \t\tx: %3d y: %3d z: %3d\n", blk[0], blk[1], blk[2]);

	unsigned long latByte      = nCells * sizeof(PdfT) * N_D3Q19;
	unsigned long latFluidByte = nFluid * sizeof(PdfT) * N_D3Q19;
	unsigned long latPadByte   = (nCells - nFluid) * sizeof(PdfT) * N_D3Q19;

	printf("# Lattice Array Allocation:\n");
	printf("#   lattice size:      \t\t%e MiB\n", latByte      / 1024.0 / 1024.0);
	printf("#   fluid lattice size:\t\t%e MiB\n", latFluidByte / 1024.0 / 1024.0);
	printf("#   lattice padding:   \t\t%e MiB\n", latPadByte   / 1024.0 / 1024.0);


	printf("#   alignment:         \t\t%d b\n", PAGE_4K);

	if (PDF_ALLOCATOR((void **)&pdfs[0], latFluidByte, PAGE_4K)) {
		printf("ERROR: allocating PDF array with %s() failed: %lu bytes.\n", STRINGIFY(PDF_ALLOCATOR), latFluidByte);
		exit(1);
	}
	else {
		printf("#   allocator: \t\t\t%s()\n", STRINGIFY(PDF_ALLOCATOR));
	}

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
	for (int bZ = 0; bZ < lZ; bZ += blk[2]) {
	for (int bY = 0; bY < lY; bY += blk[1]) {
	for (int bX = 0; bX < lX; bX += blk[0]) {

		int eX = MIN(bX + blk[0], lX);
		int eY = MIN(bY + blk[1], lY);
		int eZ = MIN(bZ + blk[2], lZ);


		for (int z = bZ; z < eZ; ++z) {
		for (int y = bY; y < eY; ++y) {
		for (int x = bX; x < eX; ++x) {

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

	// AoSoA addressing for adjList needs padding for (nFluid % VSIZE) != 0
	unsigned long adjListBytes = nFluid * sizeof(int) * N_D3Q19_IDX;

	printf("# Adjacency List Allocation:\n");
	printf("#   size:              \t\t%e MiB\n", adjListBytes / 1024.0 / 1024.0);
	printf("#   alignment:         \t\t%d b\n", PAGE_4K);

	// AdjList only requires 18 instead of 19 entries per node, as
	// the center PDF needs no addressing.
	if (ADJ_LIST_ALLOCATOR((void **)&adjList, adjListBytes, PAGE_4K)) {
		printf("ERROR: allocating adjList array with %s() failed: %lu bytes.\n", STRINGIFY(ADJ_LIST_ALLOCATOR), adjListBytes);
		exit(1);
	}
	else {
		printf("#   allocator: \t\t\t%s()\n", STRINGIFY(ADJ_LIST_ALLOCATOR));
	}

	kdl->AdjList = adjList;

	int x, y, z;

	uint32_t neighborIndex;
	uint32_t dstIndex;

	int nx, ny, nz, px, py, pz;

	// Loop over all fluid nodes and compute the indices to the neighboring
	// PDFs for configured data layout (AoS/SoA).
	// Parallelized loop to ensure correct NUMA placement.
	// #ifdef _OPENMP  --> add line continuation
	// 	#pragma omp parallel for default(none)
	// 		shared(nFluid, nCells, coords, D3Q19_INV, D3Q19_X, D3Q19_Y, D3Q19_Z,
	// 				stderr,
	// 				lDims, grid, ld, lX, lY, lZ, adjList)
	// 		private(x, y, z, nx, ny, nz, neighborIndex, dstIndex)
	// #endif
	for (int fluidBaseIndex = 0; fluidBaseIndex < nFluid; fluidBaseIndex+=VSIZE) {


		// Loop over all directions except the center one.
		for(int d = 0; d < N_D3Q19 - 1; ++d) {
			Assert(d != D3Q19_C);

			for(int inChunkIndex = 0; (inChunkIndex < VSIZE) && ((fluidBaseIndex + inChunkIndex) < nFluid); ++inChunkIndex){
				int index = fluidBaseIndex + inChunkIndex;

				Assert(index < nFluid);

				x = coords[C_INDEX_X(index)];
				y = coords[C_INDEX_Y(index)];
				z = coords[C_INDEX_Z(index)];

				Assert(x >= 0 && x < lX);
				Assert(y >= 0 && y < lY);
				Assert(z >= 0 && z < lZ);

				Assert(ld->Lattice[L_INDEX_4(lDims, x, y, z)] != LAT_CELL_OBSTACLE);

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
				// If the neighbor is outside the latice in X direction and we have a
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

				adjList[(index - (index % VSIZE)) * N_D3Q19_IDX + (d * VSIZE) + (index % VSIZE)] = dstIndex;
			}
		}
	}

	/*
	printf("============\n");
	for(int baseIndex = 0; baseIndex < nFluid; baseIndex+=VSIZE){
		for(int i = 0; i < VSIZE; ++i){
			int index = baseIndex + i;

			printf("%d ", adjList[(index - (index % VSIZE)) * N_D3Q19_IDX + (0 * VSIZE) + (index % VSIZE)]);
		}
		printf("\n");
	}
	printf("============\n");
	*/

	int nThreads = 1;

#ifdef _OPENMP
	nThreads = omp_get_max_threads();
#endif

	SetupConsecNodes(ld, KDLR(kd), nThreads);

	// Fill remaining KernelData structures
	kd->GetNode = GetNode;
	kd->SetNode = SetNode;

	kd->BoundaryConditionsGetPdf = FNAME(BCGetPdf);
	kd->BoundaryConditionsSetPdf = FNAME(BCSetPdf);

	kd->Kernel = FNAME(D3Q19ListAaPvGatherKernel);

	kd->DstPdfs = NULL;
	kd->PdfsActive = kd->Pdfs[0];

	return;
}

void FNAME(D3Q19ListAaPvGatherDeinit)(LatticeDesc * ld, KernelData ** kernelData)
{
	KernelDataListRia ** kdlr = (KernelDataListRia **)kernelData;

	MemFree((void **)&((*kdlr)->ConsecNodes));

	if ((*kdlr)->ConsecThreadIndices != NULL) {
		MemFree((void **)&((*kdlr)->ConsecThreadIndices));
	}

	KernelDataList ** kdl = (KernelDataList **)kernelData;

	ADJ_LIST_FREE((void **)&((*kdl)->AdjList));

	MemFree((void **)&((*kdl)->Coords));
	MemFree((void **)&((*kdl)->Grid));

	PDF_FREE((void **)&((*kernelData)->Pdfs[0]));

	MemFree((void **)kernelData);
	return;
}

