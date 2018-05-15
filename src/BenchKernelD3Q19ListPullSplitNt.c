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
#include "BenchKernelD3Q19ListPullSplitNtCommon.h"

#include "Memory.h"
#include "Vtk.h"
#include "Vector.h"
#include "LikwidIf.h"

#include <inttypes.h>
#include <math.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#define TMP_UX 18
#define TMP_UY 19
#define TMP_UZ 20
#define TMP_W1 21
#define TMP_W2 22

#define N_TMP 23

#define TMP_INDEX(tmp_index, tmp_dir)	nTmpArray * (tmp_dir) + (tmp_index)

void FNAME(KernelPullSplitNt1S)(LatticeDesc * ld, KernelData * kernelData, CaseData * cd)
{

	Assert(ld != NULL);
	Assert(kernelData != NULL);
	Assert(cd != NULL);

	Assert(cd->Omega > F(0.0));
	Assert(cd->Omega < F(2.0));

	KernelData        * kd   = (KernelData *)kernelData;
	KernelDataList    * kdl  = KDL(kernelData);
	KernelDataListRia * kdlr = KDLR(kernelData);

	PdfT omega = cd->Omega;
	const PdfT omegaEven = omega;

	PdfT magicParam = F(1.0) / F(12.0);
	const PdfT omegaOdd = F(1.0) / (F(0.5) + magicParam / (F(1.0) / omega - F(0.5)));


	const PdfT w_0 = F(1.0) / F( 3.0);
	const PdfT w_1 = F(1.0) / F(18.0);
	const PdfT w_2 = F(1.0) / F(36.0);

	const PdfT w_1_x3 = w_1 * F(3.0);	const PdfT w_1_nine_half = w_1 * F(9.0) / F(2.0);
	const PdfT w_2_x3 = w_2 * F(3.0);	const PdfT w_2_nine_half = w_2 * F(9.0) / F(2.0);

	const VPDFT vw_1_x3 = VSET(w_1_x3);
	const VPDFT vw_2_x3 = VSET(w_2_x3);

	const VPDFT vw_1_nine_half = VSET(w_1_nine_half);
	const VPDFT vw_2_nine_half = VSET(w_2_nine_half);

	const VPDFT vomegaEven = VSET(omegaEven);
	const VPDFT vomegaOdd  = VSET(omegaOdd);

	const VPDFT voneHalf = VSET(F(0.5));

	// uint32_t nConsecNodes = kdlr->nConsecNodes;
	// uint32_t * consecNodes = kdlr->ConsecNodes;
	// uint32_t consecIndex = 0;
	// uint32_t consecValue = 0;

	PdfT * src = kd->Pdfs[0];
	PdfT * dst = kd->Pdfs[1];
	PdfT * tmp;

	int maxIterations  = cd->MaxIterations;

	int nFluid         = kdl->nFluid;
	int nCells         = kdl->nCells;

	int nTmpArray      = kdlr->nTmpArray;

	Assert(nTmpArray % VSIZE == 0);

	uint32_t * adjList = kdl->AdjList;

	#ifdef VTK_OUTPUT
		if (cd->VtkOutput) {
			kd->PdfsActive = src;
			VtkWrite(ld, kd, cd, -1);
		}
	#endif

	#ifdef STATISTICS
		kd->PdfsActive = src;
		KernelStatistics(kd, ld, cd, 0);
	#endif

	X_KERNEL_START(kernelData);

	X_LIKWID_START("list-pull-split-nt-1s");

	#ifdef _OPENMP
		#pragma omp parallel default(none) \
			shared(nFluid, nCells, kd, kdl, adjList, src, dst, \
			cd, maxIterations, ld, tmp, nTmpArray, \
			stderr )
	#endif
	{
		uint32_t adjListIndex;

		PdfT ux, uy, uz, ui;
		VPDFT vux, vuy, vuz, vui;

		#define X(name, idx, idxinv, x, y, z)	PdfT JOIN(pdf_,name);
		D3Q19_LIST
		#undef X
		VPDFT vpdf_a, vpdf_b;

		PdfT evenPart, oddPart, dir_indep_trm, dens;
		PdfT w_1_indep, w_2_indep;
		VPDFT vevenPart, voddPart;
		VPDFT vw_1_indep, vw_2_indep;

		int indexMax;

		PdfT * tmpArray;
		MemAllocAligned((void **)&tmpArray, sizeof(PdfT) * nTmpArray * N_TMP, VSIZE * sizeof(PdfT));

		int nThreads = 1;
		int threadId = 0;

#ifdef _OPENMP
		nThreads = omp_get_max_threads();
		threadId = omp_get_thread_num();
#endif

		int nCellsThread = nFluid / nThreads;
		int blIndexStart = threadId * nCellsThread;

		if (threadId < nFluid % nThreads) {
			blIndexStart += threadId;
			nCellsThread += 1;
		}
		else {
			blIndexStart += nFluid % nThreads;
		}

		int blIndexStop = blIndexStart + nCellsThread;

		// We have three loops:
		// 1. Peeling to ensure alignment for non-temporal stores in loop 2 is correct.
		// 2. Vectorized handling of nodes.
		// 3. Remaining nodes, less than vector size.

		unsigned long addrStart = (unsigned long)&(src[P_INDEX_3(nCells, blIndexStart, 0)]);
		int nCellsUnaligned = (VSIZE - (int)((addrStart / sizeof(PdfT)) % VSIZE)) % VSIZE;

		int nCellsVectorized = nCellsThread - nCellsUnaligned;
		nCellsVectorized = nCellsVectorized - (nCellsVectorized % VSIZE);

		int blIndexVec       = blIndexStart + nCellsUnaligned;
		int blIndexRemaining = blIndexStart + nCellsUnaligned + nCellsVectorized;

		// printf("%d [%d, %d, %d, %d[\n", threadId, blIndexStart, blIndexVec, blIndexRemaining, blIndexStop);

		for(int iter = 0; iter < maxIterations; ++iter) {


#if 1
			#define INDEX_START	blIndexStart
			#define INDEX_STOP  blIndexVec
			#include "BenchKernelD3Q19ListPullSplitNt1SScalar.h"

			#define INDEX_START blIndexVec
			#define INDEX_STOP  blIndexRemaining
			#include "BenchKernelD3Q19ListPullSplitNt1SIntrinsics.h"

			#define INDEX_START blIndexRemaining
			#define INDEX_STOP	blIndexStop
			#include "BenchKernelD3Q19ListPullSplitNt1SScalar.h"
#else
			#define INDEX_START blIndexStart
			#define INDEX_STOP	blIndexStop
			#include "BenchKernelD3Q19ListPullSplitNt1SScalar.h"
#endif


			#pragma omp barrier

			#pragma omp single
			{
				#ifdef VERIFICATION
					kd->PdfsActive = dst;
					KernelAddBodyForce(kd, ld, cd);
				#endif

				#ifdef VTK_OUTPUT
					if (cd->VtkOutput && (iter % cd->VtkModulus) == 0) {
						kd->PdfsActive = dst;
						VtkWrite(ld, kd, cd, iter);
					}
				#endif

				#ifdef STATISTICS
					kd->PdfsActive = dst;
					KernelStatistics(kd, ld, cd, iter);
				#endif

				// swap grids
				tmp = src;
				src = dst;
				dst = tmp;
			}

			#pragma omp barrier

		} // for (int iter = 0; ...

		MemFree((void **)&tmpArray);
	}


	X_LIKWID_STOP("list-pull-split-nt-1s");

	X_KERNEL_END(kernelData);

#ifdef VTK_OUTPUT
	if (cd->VtkOutput) {
		kd->PdfsActive = src;
		VtkWrite(ld, kd, cd, maxIterations);
	}
#endif

#ifdef STATISTICS
	kd->PdfsActive = src;
	KernelStatistics(kd, ld, cd, maxIterations);
#endif

	return;
}

void FNAME(KernelPullSplitNt2S)(LatticeDesc * ld, KernelData * kernelData, CaseData * cd)
{

	Assert(ld != NULL);
	Assert(kernelData != NULL);
	Assert(cd != NULL);

	Assert(cd->Omega > F(0.0));
	Assert(cd->Omega < F(2.0));

	KernelData        * kd   = (KernelData *)kernelData;
	KernelDataList    * kdl  = KDL(kernelData);
	KernelDataListRia * kdlr = KDLR(kernelData);

	PdfT omega = cd->Omega;
	const PdfT omegaEven = omega;

	PdfT magicParam = F(1.0) / F(12.0);
	const PdfT omegaOdd = F(1.0) / (F(0.5) + magicParam / (F(1.0) / omega - F(0.5)));

	const PdfT w_0 = F(1.0) / F( 3.0);
	const PdfT w_1 = F(1.0) / F(18.0);
	const PdfT w_2 = F(1.0) / F(36.0);

	const PdfT w_1_x3 = w_1 * F(3.0);	const PdfT w_1_nine_half = w_1 * F(9.0) / F(2.0);
	const PdfT w_2_x3 = w_2 * F(3.0);	const PdfT w_2_nine_half = w_2 * F(9.0) / F(2.0);

	const VPDFT vw_1_x3 = VSET(w_1_x3);
	const VPDFT vw_2_x3 = VSET(w_2_x3);

	const VPDFT vw_1_nine_half = VSET(w_1_nine_half);
	const VPDFT vw_2_nine_half = VSET(w_2_nine_half);

	const VPDFT vomegaEven = VSET(omegaEven);
	const VPDFT vomegaOdd  = VSET(omegaOdd);

	const VPDFT voneHalf = VSET(F(0.5));

	// uint32_t nConsecNodes = kdlr->nConsecNodes;
	// uint32_t * consecNodes = kdlr->ConsecNodes;
	// uint32_t consecIndex = 0;
	// uint32_t consecValue = 0;

	PdfT * src = kd->Pdfs[0];
	PdfT * dst = kd->Pdfs[1];
	PdfT * tmp;

	int maxIterations  = cd->MaxIterations;

	int nFluid         = kdl->nFluid;
	int nCells         = kdl->nCells;

	int nTmpArray      = kdlr->nTmpArray;

	Assert(nTmpArray % VSIZE == 0);

	uint32_t * adjList = kdl->AdjList;

	#ifdef VTK_OUTPUT
		if (cd->VtkOutput) {
			kd->PdfsActive = src;
			VtkWrite(ld, kd, cd, -1);
		}
	#endif

	#ifdef STATISTICS
		kd->PdfsActive = src;
		KernelStatistics(kd, ld, cd, 0);
	#endif


	X_KERNEL_START(kernelData);

	X_LIKWID_START("list-pull-split-nt-2s");


	#ifdef _OPENMP
		#pragma omp parallel default(none) \
			shared(nFluid, nCells, kd, kdl, adjList, src, dst, \
			cd, maxIterations, ld, tmp, nTmpArray, \
			stderr )
	#endif
	{
		uint32_t adjListIndex;

		PdfT ux, uy, uz, ui;
		VPDFT vux, vuy, vuz, vui;

		#define X(name, idx, idxinv, x, y, z)	PdfT JOIN(pdf_,name);
		D3Q19_LIST
		#undef X
		VPDFT vpdf_a, vpdf_b;

		PdfT evenPart, oddPart, dir_indep_trm, dens;
		PdfT w_1_indep, w_2_indep;
		VPDFT vevenPart, voddPart;
		VPDFT vw_1_indep, vw_2_indep;

		int indexMax;

		PdfT * tmpArray;
		MemAlloc((void **)&tmpArray, sizeof(PdfT) * nTmpArray * N_TMP);

		int nThreads = 1;
		int threadId = 0;

#ifdef _OPENMP
		nThreads = omp_get_max_threads();
		threadId = omp_get_thread_num();
#endif

		int nCellsThread = nFluid / nThreads;
		int blIndexStart = threadId * nCellsThread;

		if (threadId < nFluid % nThreads) {
			blIndexStart += threadId;
			nCellsThread += 1;
		}
		else {
			blIndexStart += nFluid % nThreads;
		}

		int blIndexStop = blIndexStart + nCellsThread;

		// We have three loops:
		// 1. Peeling to ensure alignment for non-temporal stores in loop 2 is correct.
		// 2. Vectorized handling of nodes.
		// 3. Remaining nodes, less than vector size.

		unsigned long addrStart = (unsigned long)&(src[P_INDEX_3(nCells, blIndexStart, 0)]);
		int nCellsUnaligned = (VSIZE - (int)((addrStart / sizeof(PdfT)) % VSIZE)) % VSIZE;

		int nCellsVectorized = nCellsThread - nCellsUnaligned;
		nCellsVectorized = nCellsVectorized - (nCellsVectorized % VSIZE);

		int blIndexVec       = blIndexStart + nCellsUnaligned;
		int blIndexRemaining = blIndexStart + nCellsUnaligned + nCellsVectorized;

		// printf("%d [%d, %d, %d, %d[\n", threadId, blIndexStart, blIndexVec, blIndexRemaining, blIndexStop);

		for(int iter = 0; iter < maxIterations; ++iter) {

#if 1
			#define INDEX_START	blIndexStart
			#define INDEX_STOP  blIndexVec
			#include "BenchKernelD3Q19ListPullSplitNt2SScalar.h"

			#define INDEX_START blIndexVec
			#define INDEX_STOP  blIndexRemaining
			#include "BenchKernelD3Q19ListPullSplitNt2SIntrinsics.h"

			#define INDEX_START blIndexRemaining
			#define INDEX_STOP	blIndexStop
			#include "BenchKernelD3Q19ListPullSplitNt2SScalar.h"
#else
			#define INDEX_START blIndexStart
			#define INDEX_STOP	blIndexStop
			#include "BenchKernelD3Q19ListPullSplitNt2SScalar.h"
#endif
			#pragma omp barrier


			#pragma omp single
			{
				#ifdef VERIFICATION
					kd->PdfsActive = dst;
					KernelAddBodyForce(kd, ld, cd);
				#endif

				#ifdef VTK_OUTPUT
					if (cd->VtkOutput && (iter % cd->VtkModulus) == 0) {
						kd->PdfsActive = dst;
						VtkWrite(ld, kd, cd, iter);
					}
				#endif

				#ifdef STATISTICS
					kd->PdfsActive = dst;
					KernelStatistics(kd, ld, cd, iter);
				#endif

				// swap grids
				tmp = src;
				src = dst;
				dst = tmp;
			}

			#pragma omp barrier

		} // for (int iter = 0; ...

		MemFree((void **)&tmpArray);
	}

	X_LIKWID_STOP("list-pull-split-nt-2s");

	X_KERNEL_END(kernelData);

#ifdef VTK_OUTPUT
	if (cd->VtkOutput) {
		kd->PdfsActive = src;
		VtkWrite(ld, kd, cd, maxIterations);
	}
#endif

#ifdef STATISTICS
	kd->PdfsActive = src;
	KernelStatistics(kd, ld, cd, maxIterations);
#endif

	return;
}

