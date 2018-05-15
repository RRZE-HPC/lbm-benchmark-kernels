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
#include "BenchKernelD3Q19ListAaPvGatherHybridCommon.h"

#include "Memory.h"
#include "Vtk.h"
#include "Vector.h"

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
void FNAME(D3Q19ListAaPvGatherHybridKernel)(LatticeDesc * ld, struct KernelData_ * kd, CaseData * cd);



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

static void SetuploopStartIndices(LatticeDesc * ld, KernelDataListRia * kdlr, int nThreads)
{
	//#define ADJ_LIST(dir) adjList[(index - (index % VSIZE)) * N_D3Q19_IDX + (dir * VSIZE) + (index % VSIZE)]
	Assert(ld != NULL);
	Assert(kdlr != NULL);
	Assert(nThreads > 0);

	//uint32_t * adjList = kdlr->kdl.AdjList;
	uint32_t * adjList = kdlr->kdl.AdjList;

	uint32_t nLoopStartIndices = 0;
	uint32_t loopStartIndex = 2;

	int nFluid = kdlr->kdl.nFluid;

	int * oddKernelThreadStartIndices = (int *)malloc(sizeof(int) * (nThreads + 1));

	int nNodesPerThread = nFluid / nThreads;
	//printf("nodesPerThread: %d\n", nNodesPerThread);

	for (int i = 0; i < nThreads; ++i) {
		oddKernelThreadStartIndices[i] = i * nNodesPerThread + MinI(i, nFluid % nThreads);
	}

	oddKernelThreadStartIndices[nThreads] = nFluid;

	/*
	   for (int i = 0; i <= nThreads; ++i) {
	   printf("oddKernelThreadStartIndices[%d] = %d\n", i, oddKernelThreadStartIndices[i]);
	   }
	   */

	int threadIndex = 1;

	// We execute following code two times.
	// - The first time to get the count of how many entries we need for the
	//   loopStartIndices array.
	// - The second time to fill the array.

	// Loop over adjacency list of all nodes.
	// Compare if adjacent nodes within one cache line share the same access pattern.
	// First vectorized access is assumed to be consecutive (-> may be loaded with regular load).

	int lastCacheLineConsecutive = 1;

	for (int fluidBaseIndex = 1; fluidBaseIndex < nFluid + 1; fluidBaseIndex += VSIZE) {

		int currentCacheLineConsecutive = 1;

		// Loop over all directions except the center one.
		for(int d = 0; d < N_D3Q19 - 1; ++d) {
			Assert(d != D3Q19_C);

			// check if cache line itself has consecutive memory access pattern
			for(int inChunkIndex = 0; (inChunkIndex < VSIZE - 1) && ((fluidBaseIndex + inChunkIndex) < nFluid); ++inChunkIndex) {
				int index = fluidBaseIndex + inChunkIndex;

				Assert(index < nFluid);

				#define ADJ_LIST(idx, dir) adjList[((idx) - ((idx) % VSIZE)) * N_D3Q19_IDX + ((dir) * VSIZE) + ((idx) % VSIZE)]
				//if (adjList[index * N_D3Q19_IDX + d] != adjList[(index - 1) * N_D3Q19_IDX + d] + 1)
				if (ADJ_LIST(index, d) != ADJ_LIST(index-1, d) + 1) {
					//printf("no match for index: %d\n", d);
					//printf("ADJ_LlST(%d,%d) = %d  !=  %d = ADJ_LlST(%d,%d) + 1\n", index, d, ADJ_LIST(index,d), ADJ_LIST(index-1,d), index-1, d);
					// Different access pattern.
					currentCacheLineConsecutive = 0;
					break;
				}
				#undef ADJ_LIST

			}

			if(!currentCacheLineConsecutive){
				break; //exit from nested loop
			}

		}

		int interCacheLineConsecutive = 1;

		if(currentCacheLineConsecutive && lastCacheLineConsecutive){
			// check if cache line has consecutive memory access pattern to last entry of previous cache line
			int lastIdxOfPreviousCacheLine = fluidBaseIndex - 2;
			if (lastIdxOfPreviousCacheLine > 0) {
				for(int d = 0; d < N_D3Q19 - 1; ++d) {
					Assert(d != D3Q19_C);
					#define ADJ_LIST(idx, dir) adjList[((idx) - ((idx) % VSIZE)) * N_D3Q19_IDX + ((dir) * VSIZE) + ((idx) % VSIZE)]
					if (ADJ_LIST(fluidBaseIndex-1, d) != ADJ_LIST(lastIdxOfPreviousCacheLine, d) + 1) {
						// Different access pattern.
						//printf("not interCacheConsecutive\n");
						interCacheLineConsecutive = 0;
						break;
					}
					#undef ADJ_LIST

				}
			}
		}
		int threadBoundaryIndex = oddKernelThreadStartIndices[threadIndex];
		if (fluidBaseIndex - 1 <= threadBoundaryIndex &&
				threadBoundaryIndex < fluidBaseIndex + VSIZE - 1) {
			// Current cache line contains thread boundary.
			// These cache lines are treated by scalar peel and
			// reminder loops in kernel of every thread.
			// TODO maybe replace these loops with masked gather / scatter?!
			if (loopStartIndex % 2 == 0) { // current index would be gather/scatter index
				++loopStartIndex; // reserving gather/scatter index
			}
			++loopStartIndex; // reserving space for start remainder loop of thread n

			if (threadIndex < nThreads){
				++loopStartIndex; // reserving space for starting peel loop of thread n+1

				if (fluidBaseIndex - 1 == threadBoundaryIndex){
					if(!currentCacheLineConsecutive){
						++loopStartIndex;
					}
				}
				else {
					currentCacheLineConsecutive = 1;
				}

				//++loopStartIndex; // reserving space for ending peel loop / starting load/store of thread n+1
				++loopStartIndex; // 1st load/store end == 1st gather/scatter start OR 1st gather/scatter end == 2nd load/start start
			}
			++threadIndex;
		}
		else {
			// We are not at a thread boundary.
			if (currentCacheLineConsecutive) {
				if(lastCacheLineConsecutive && !interCacheLineConsecutive){
					loopStartIndex+=2;
				}
				else if(!lastCacheLineConsecutive){
					++loopStartIndex;
				}
			}
			else {
				if(lastCacheLineConsecutive){
					++loopStartIndex;
				}
			}
		}

		// treating special case when last thread has no remainder loop
		if (oddKernelThreadStartIndices[nThreads] == fluidBaseIndex + VSIZE - 1) {
			//printf("--> special case 111. loopStartIndex: %d \n", loopStartIndex);
			if (loopStartIndex % 2 != 0) { // current index is gather/scatter end and load/store start index
				++loopStartIndex; //set load/store end (gather/scatter start) to same value as scalar remainder start => no more access to gather/scatter loop
			}

			++loopStartIndex; // gather/scatter end and scalar remainder start
			++loopStartIndex; // scalar remainder end and scalar peel start

		}

		lastCacheLineConsecutive = currentCacheLineConsecutive;
	}

	if (nFluid > 0) {
		nLoopStartIndices = loopStartIndex;
	}

	int * loopStartIndices;
	unsigned long loopStartIndicesByte = (nLoopStartIndices + 1) * sizeof(int);

	printf("# Loop Start Index Array Allocation:\n");
	printf("#   elements:      \t\t%d\n",   nLoopStartIndices + 1);
	printf("#   size:          \t\t%e MiB\n", loopStartIndicesByte / 1024.0 / 1024.0);
	printf("#   alignment:     \t\t%d b\n",   PAGE_4K);

	if (MemAllocAligned((void **)&loopStartIndices, loopStartIndicesByte, PAGE_4K)) {
		printf("ERROR: allocating loopStartIndices array with MemAllocAligned failed: %lu bytes.\n", loopStartIndicesByte);
		exit(1);
	}
	else {
		printf("#   allocator: \t\t\tMemAllocAligned()\n");
	}

	oddKernelThreadStartIndices[0] = 0;
	loopStartIndices[0] = 0; //first scalar loop would start with 0
	loopStartIndices[1] = 0; //no peel loop expected -> first load/store loop may start at index==0
	loopStartIndices[2] = 0; //may not be set in case first access is gather/scatter -> therefore its set here

	// resetting values to default
	threadIndex = 1;
	lastCacheLineConsecutive = 1;
	loopStartIndex = 2;

	// Loop over adjacency list of all nodes.
	// Compare if adjacent nodes share the same access pattern.

	int indexAccumulator = 0;

	// for statistical reasons:
	int gatherAccumulator = 0;
	int loadAccumulator = 0;
	int scalarLookups = 0;
	int loadLookups = 0;


	for (int fluidBaseIndex = 1; fluidBaseIndex < nFluid + 1; fluidBaseIndex += VSIZE) {
		int currentCacheLineConsecutive = 1;
		//printf("fluidbaseIndex: %d\n", fluidBaseIndex);
		// Loop over all directions except the center one.
		for(int d = 0; d < N_D3Q19 - 1; ++d) {
			Assert(d != D3Q19_C);

			// check if cache line itself has consecutive memory access pattern
			for(int inChunkIndex = 0; (inChunkIndex < VSIZE - 1) && ((fluidBaseIndex + inChunkIndex) < nFluid); ++inChunkIndex){
				int index = fluidBaseIndex + inChunkIndex;

				Assert(index < nFluid);

				#define ADJ_LIST(idx, dir) adjList[((idx) - ((idx) % VSIZE)) * N_D3Q19_IDX + ((dir) * VSIZE) + ((idx) % VSIZE)]
				//if (adjList[index * N_D3Q19_IDX + d] != adjList[(index - 1) * N_D3Q19_IDX + d] + 1)
				if (ADJ_LIST(index, d) != ADJ_LIST(index-1, d) + 1) {
					// Different access pattern.
					currentCacheLineConsecutive = 0;
					break;
				}
				#undef ADJ_LIST
			}

			if(!currentCacheLineConsecutive){
				break; //exit from nested loop
			}
		}

		int interCacheLineConsecutive = 1;

		if(currentCacheLineConsecutive && lastCacheLineConsecutive){
			// check if cache line has consecutive memory access pattern to last entry of previous cache line
			int lastIdxOfPreviousCacheLine = fluidBaseIndex - 2;
			if (lastIdxOfPreviousCacheLine > 0) {
				for(int d = 0; d < N_D3Q19 - 1; ++d) {
					Assert(d != D3Q19_C);
					#define ADJ_LIST(idx, dir) adjList[((idx) - ((idx) % VSIZE)) * N_D3Q19_IDX + ((dir) * VSIZE) + ((idx) % VSIZE)]
					if (ADJ_LIST(fluidBaseIndex-1, d) != ADJ_LIST(lastIdxOfPreviousCacheLine, d) + 1) {
						// Different access pattern.
						interCacheLineConsecutive = 0;
						break;
					}
					#undef ADJ_LIST

				}
			}
		}

		int threadBoundaryIndex = oddKernelThreadStartIndices[threadIndex];
		//if (fluidBaseIndex > 3500)
		//	printf("threadBoundaryIndex: %d  fluidBaseIndex-1: %d fluidBaseIndex + VSIZE - 1: %d\n", threadBoundaryIndex, fluidBaseIndex-1, fluidBaseIndex + VSIZE -1);

		if (fluidBaseIndex - 1 <= threadBoundaryIndex &&
				threadBoundaryIndex < fluidBaseIndex + VSIZE - 1) {
			// Current cache line contains thread boundary.
			// These cache lines are treated by scalar peel and
			// reminder loops in kernel of every thread.
			// TODO maybe replace these loops with masked gather / scatter?!
			if (loopStartIndex % 2 == 0) { // current index would be gather/scatter index
				//loopStartIndices[loopStartIndex] = fluidBaseIndex - 1; //same value as scalar remainder start => no more access to gather/scatter loop
				loopStartIndices[loopStartIndex] = indexAccumulator; //same value as scalar remainder start => no more access to gather/scatter loop
				++loopStartIndex;
			}

			//loopStartIndices[loopStartIndex] = fluidBaseIndex - 1; // gather/scatter end and scalar remainder start
			loopStartIndices[loopStartIndex] = indexAccumulator; // gather/scatter end and scalar remainder start
			++loopStartIndex;

			// starting indices of thread n+1
			loopStartIndices[loopStartIndex] = threadBoundaryIndex; // scalar remainder of thread n end and scalar peel of thread n+1 start
			oddKernelThreadStartIndices[threadIndex] = loopStartIndex; // thread start is where scalar peel starts

			if (threadIndex < nThreads){
				indexAccumulator = ((threadBoundaryIndex + VSIZE - 1) / VSIZE ) * VSIZE; // rounding towards next multiple of VSIZE
				++loopStartIndex;
				loopStartIndices[loopStartIndex] = indexAccumulator; // scalar peel end and 1st load/store start

				// treating special case when there is no peel / remainder loop
				if (fluidBaseIndex - 1 == threadBoundaryIndex){
					if(!currentCacheLineConsecutive){
						++loopStartIndex;
						loopStartIndices[loopStartIndex] = indexAccumulator; // 1st load/store end and 1st gather/scatter start
						gatherAccumulator += VSIZE;
					}
					else {
						loadLookups += VSIZE;
					}
					indexAccumulator += VSIZE;
				}
				else {
					scalarLookups += VSIZE;
					currentCacheLineConsecutive = 1;
				}

				++loopStartIndex; // 1st load/store end == 1st gather/scatter start OR 1st gather/scatter end == 2nd load/start start
				loopStartIndices[loopStartIndex] = indexAccumulator; // 1st load/store end == 1st gather/scatter start OR 1st gather/scatter end == 2nd load/start start
			}
			++threadIndex;

		}
		else {
			// We are not at a thread boundary.
			int print = 0;
			if (currentCacheLineConsecutive) {
				loadAccumulator += VSIZE;

				if(lastCacheLineConsecutive && !interCacheLineConsecutive){
					loadLookups += VSIZE;
					if (print)
						printf("#1 loopStartIndex: %d\n", loopStartIndex);
					// loopStartIndices[loopStartIndex] is not incremented since pointers need to be fetched again.
					// loopStartIndices[loopStartIndex + 1] (-> start Load/Store and end Gather/Scatter)
					// gets same value as loopStartIndices[loopStartindex] (-> start Gather/Scatter)
					// this ensures that no gather/scatter iteration is executed
					++loopStartIndex;
					loopStartIndices[loopStartIndex] = indexAccumulator;

					// loopStartIndices[loopStartIndex + 2] (-> start Gather/Scatter and end Load/Store)
					// gets set to have one Load/Store iteration
					++loopStartIndex;
					indexAccumulator+=VSIZE;
					loopStartIndices[loopStartIndex] = indexAccumulator;

				}
				else if(!lastCacheLineConsecutive){
					loadLookups += VSIZE;
					if (print)
						printf("#2 loopStartIndex: %d\n", loopStartIndex);
					++loopStartIndex;
					indexAccumulator+=VSIZE;
					loopStartIndices[loopStartIndex] = indexAccumulator;
				}
				else { // (lastCacheLineConsecutive && interCacheLineConsecutive)
					if (print)
						printf("#3 loopStartIndex: %d\n", loopStartIndex);
					indexAccumulator+=VSIZE;
					loopStartIndices[loopStartIndex] = indexAccumulator;
				}
			}
			else {
				gatherAccumulator += VSIZE;
				if(lastCacheLineConsecutive){
					if (print)
						printf("#4 loopStartIndex: %d\n", loopStartIndex);
					++loopStartIndex;
					indexAccumulator+=VSIZE;
					loopStartIndices[loopStartIndex] = indexAccumulator;
				}
				else { // lastCacheLine without not consecutive memory access pattern
					if (print)
						printf("#5 loopStartIndex: %d\n", loopStartIndex);
					indexAccumulator+=VSIZE;
					loopStartIndices[loopStartIndex] = indexAccumulator;
				}
			}
		}

		// treating special case when last thread has no remainder loop
		if (oddKernelThreadStartIndices[nThreads] == fluidBaseIndex + VSIZE - 1) {
			//printf("--> special case. indexAccumulator: %d\n", indexAccumulator);
			if (loopStartIndex % 2 != 0) { // current index is gather/scatter end and load/store start index
				++loopStartIndex;
				loopStartIndices[loopStartIndex] = indexAccumulator; //set load/store end (gather/scatter start) to same value as scalar remainder start => no more access to gather/scatter loop
			}

			++loopStartIndex;
			loopStartIndices[loopStartIndex] = indexAccumulator; // gather/scatter end and scalar remainder start
			++loopStartIndex;
			loopStartIndices[loopStartIndex] = indexAccumulator; // scalar remainder end and scalar peel start

			oddKernelThreadStartIndices[threadIndex] = loopStartIndex; // thread start is where scalar peel starts
		}

		lastCacheLineConsecutive = currentCacheLineConsecutive;

	}

	if (nLoopStartIndices != loopStartIndex){
		printf("ERROR: nLoopStartIndices unequal loopStartIndex!\n");
	}

	/*
	printf("loopStartIndices:\n");
	for(int i = 0; i <= nLoopStartIndices; ++i){
		printf("%d ", loopStartIndices[i]);
	}
	printf("\n");
	printf("oddKernelThreadStartIndices:\n");
	for(int i = 0; i <= nThreads; ++i){
		printf("%d ", oddKernelThreadStartIndices[i]);
	}
	printf("\n");
	*/

	kdlr->loopStartIndices = loopStartIndices;
	kdlr->nLoopStartIndices = nLoopStartIndices;

	kdlr->oddKernelThreadStartIndices  = oddKernelThreadStartIndices;
	kdlr->nOddKernelThreadStartIndices = nThreads;

	printf("#   vload/vstore nodes:  \t% 10d   \t(%3.4f %% of total fluid nodes)\n", loadAccumulator, ((double) loadAccumulator / (double) nFluid) * 100);
	printf("#   gather/scatter nodes:\t% 10d   \t(%3.4f %% of total fluid nodes)\n", gatherAccumulator, ((double) gatherAccumulator / (double) nFluid) * 100.0);
	printf("#   vload/vstore lookups:\t% 10d \n", loadLookups * (N_D3Q19 - 1));
	printf("#   gather/scatter lookups:\t% 10d \n", gatherAccumulator * (N_D3Q19 - 1));
	printf("#   scalar lookups:      \t% 10d \n", scalarLookups * (N_D3Q19 - 1));

	double loopBalanceEven = 2.0 * 19 * sizeof(PdfT);
	double loopBalanceOdd  = 2.0 * 19 * sizeof(PdfT) /* actual PDFs */
		+ (((double)(gatherAccumulator + loadLookups + scalarLookups)) / nFluid) * sizeof(int) * (N_D3Q19 - 1) /* AdjList */
		+ (nLoopStartIndices / nFluid) * sizeof(int); // one lookup to loopStartIndices

	double loopBalance     = (loopBalanceEven + loopBalanceOdd) / 2.0;

	kdlr->kdl.kd.LoopBalance = loopBalance;

	printf("# loop balance:\n");
	printf("#   even timestep:  \t\t%.2f B/FLUP\n", loopBalanceEven);
	printf("#   odd timestep:   \t\t%.2f B/FLUP\n", loopBalanceOdd);
	printf("#   average:        \t\t%.2f B/FLUP\n", loopBalance);

	return;
}

void FNAME(D3Q19ListAaPvGatherHybridInit)(LatticeDesc * ld, KernelData ** kernelData, Parameters * params)
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

	kdlr->loopStartIndices = NULL;
	kdlr->nLoopStartIndices = 0;
	kdlr->oddKernelThreadStartIndices = NULL;
	kdlr->nOddKernelThreadStartIndices = 0;
#endif

	kdl->Iteration = -1;

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
	int nFluid_padded = ((nFluid + VSIZE - 1) / VSIZE) * VSIZE;
	unsigned long adjListBytes = nFluid_padded * sizeof(int) * N_D3Q19_IDX;

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

	for (int i = 0; i < nFluid_padded; ++i){
		adjList[i] = -1;
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

// Sets unused adjList entries to some extreme value which triggers and SIGSEG, whenever these values are accidently accessed.
	for(int index = nFluid; index < nFluid_padded; ++index){
		for(int d = 0; d < N_D3Q19 - 1; ++d) {
			adjList[(index - (index % VSIZE)) * N_D3Q19_IDX + (d * VSIZE) + (index % VSIZE)] = -10*1000*1000;
		}
	}

/*
	printf("============\n");
	for (int i = 0; i < nFluid_padded * (N_D3Q19_IDX + 20);){
		for (int j = 0; j < VSIZE; ++j){
			printf("%d ",adjList[i]);
			++i;
		}
		printf("\n");
	}
	   for(int dir = 0; dir < N_D3Q19; ++dir){
	   printf("dir: %d\n",dir);
	   for(int baseIndex = 0; baseIndex < nFluid + VSIZE; baseIndex+=VSIZE){
	   for(int i = 0; i < VSIZE; ++i){
	   int index = baseIndex + i;

	   printf("%d ", adjList[(index - (index % VSIZE)) * N_D3Q19_IDX + (dir * VSIZE) + (index % VSIZE)]);
	   }
	   printf("\n");
	   }
	   printf("\n");
	   }
	printf("============\n");
*/


	int nThreads = 1;

#ifdef _OPENMP
	nThreads = omp_get_max_threads();
#endif

	SetuploopStartIndices(ld, KDLR(kd), nThreads);

	// Fill remaining KernelData structures
	kd->GetNode = GetNode;
	kd->SetNode = SetNode;

	kd->BoundaryConditionsGetPdf = FNAME(BCGetPdf);
	kd->BoundaryConditionsSetPdf = FNAME(BCSetPdf);

	kd->Kernel = FNAME(D3Q19ListAaPvGatherHybridKernel);

	kd->DstPdfs = NULL;
	kd->PdfsActive = kd->Pdfs[0];

	return;
}

void FNAME(D3Q19ListAaPvGatherHybridDeinit)(LatticeDesc * ld, KernelData ** kernelData)
{
	KernelDataListRia ** kdlr = (KernelDataListRia **)kernelData;

	MemFree((void **)&((*kdlr)->loopStartIndices));

	if ((*kdlr)->oddKernelThreadStartIndices != NULL) {
		MemFree((void **)&((*kdlr)->oddKernelThreadStartIndices));
	}

	KernelDataList ** kdl = (KernelDataList **)kernelData;

	ADJ_LIST_FREE((void **)&((*kdl)->AdjList));

	MemFree((void **)&((*kdl)->Coords));
	MemFree((void **)&((*kdl)->Grid));

	PDF_FREE((void **)&((*kernelData)->Pdfs[0]));

	MemFree((void **)kernelData);
	return;
}
#undef PAGE_4K
#undef ADJ_LIST
