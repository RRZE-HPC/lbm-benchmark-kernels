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
#include "LikwidIf.h"
#include "Vector.h"
#include "Vector.h"

#include <inttypes.h>
#include <math.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

static void KernelEven(LatticeDesc * ld, KernelData * kd, CaseData * cd);
static void KernelOddVecSl(LatticeDesc * ld, KernelData * kd, CaseData * cd);

#if 1 // {{{
void DumpPdfs(LatticeDesc * ld, KernelData * kd, int zStart, int zStop, int iter, const char * prefix, int dir)
{
	int * gDims = kd->GlobalDims;

	int nX = gDims[0];
	int nY = gDims[1];
	// int nZ = gDims[2];

	PdfT pdfs[N_D3Q19];

	int localZStart = zStart;
	int localZStop  = zStop;

	if (localZStart == -1) localZStart = 0;
	if (localZStop  == -1) localZStop  = gDims[2] - 1;

	printf("D iter: %d dir: %d %s\n", iter, dir, D3Q19_NAMES[dir]);

//	for (int dir = 0; dir < 19; ++dir) {
		for (int z = localZStop; z >= localZStart; --z) {
			printf("D [%2d][%2d][%s] plane % 2d\n", iter, dir, prefix, z);

			for(int y = 0; y < nY; ++y) {
			// for(int y = 2; y < nY - 2; ++y) {
				printf("D [%2d][%2d][%s] %2d  ", iter, dir, prefix, y);

				for(int x = 0; x < nX; ++x) {

					if (1) { // ld->Lattice[L_INDEX_4(ld->Dims, x, y, z)] != LAT_CELL_OBSTACLE) {

						#define I(x, y, z, dir)	P_INDEX_5(gDims, (x), (y), (z), (dir))
						pdfs[dir] = kd->PdfsActive[I(x, y, z, dir)];
						#undef I
					}
					else {
						pdfs[dir] = -1.0;
					}

					printf("%.16e ", pdfs[dir]);
					// printf("%08.0f ", pdfs[dir]);
				}

				printf("\n");
			}
		}
//	}
}
#endif  // }}}

void FNAME(D3Q19AaVecSlKernel)(LatticeDesc * ld, KernelData * kd, CaseData * cd)
{
	Assert(ld != NULL);
	Assert(kd != NULL);
	Assert(cd != NULL);

	Assert(cd->Omega > 0.0);
	Assert(cd->Omega < 2.0);

	KernelDataAa * kda = KDA(kd);

	PdfT * src = kd->PdfsActive;

	int maxIterations = cd->MaxIterations;

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

	Assert((maxIterations % 2) == 0);

	X_KERNEL_START(kd);

	#ifdef _OPENMP
		#pragma omp parallel default(none) shared(kda, kd, ld, cd, src, maxIterations)
	#endif
	{
		for (int iter = 0; iter < maxIterations; iter += 2) {

			// --------------------------------------------------------------------
			// even time step
			// --------------------------------------------------------------------

			X_LIKWID_START("aa-vec-even");

			KernelEven(ld, kd, cd);
			#ifdef _OPENMP
				#pragma omp barrier
			#endif

			X_LIKWID_STOP("aa-vec-even");

			// Fixup bounce back PDFs.
			#ifdef _OPENMP
				#pragma omp for
			#endif
			#ifdef INTEL_OPT_DIRECTIVES
				#pragma ivdep
			#endif
			for (int i = 0; i < kd->nBounceBackPdfs; ++i) {
				src[kd->BounceBackPdfsSrc[i]] = src[kd->BounceBackPdfsDst[i]];
			}

			#ifdef _OPENMP
				#pragma omp single
			#endif
			{
				// save current iteration
				kda->Iteration = iter;

				#ifdef VERIFICATION
					kd->PdfsActive = src;
					KernelAddBodyForce(kd, ld, cd);
				#endif

				#ifdef VTK_OUTPUT
					if (cd->VtkOutput && (iter % cd->VtkModulus) == 0) {
						kd->PdfsActive = src;
						VtkWrite(ld, kd, cd, iter);
					}
				#endif

				#ifdef STATISTICS
					kd->PdfsActive = src;
					KernelStatistics(kd, ld, cd, iter);
				#endif
			}
			#ifdef _OPENMP
				#pragma omp barrier
			#endif


			// --------------------------------------------------------------------
			// odd time step
			// --------------------------------------------------------------------

			X_LIKWID_START("aa-vec-odd");


			KernelOddVecSl(ld, kd, cd);
			#ifdef _OPENMP
				#pragma omp barrier
			#endif

			// Stop counters before bounce back. Else computing loop balance will
			// be incorrect.

			X_LIKWID_STOP("aa-vec-odd");

			// Fixup bounce back PDFs.
			#ifdef _OPENMP
				#pragma omp for
			#endif
			#ifdef INTEL_OPT_DIRECTIVES
				#pragma ivdep
			#endif
			for (int i = 0; i < kd->nBounceBackPdfs; ++i) {
				src[kd->BounceBackPdfsDst[i]] = src[kd->BounceBackPdfsSrc[i]];
			}

			#ifdef _OPENMP
				#pragma omp single
			#endif
			{
				// save current iteration
				kda->Iteration = iter + 1;

				#ifdef VERIFICATION
					kd->PdfsActive = src;
					KernelAddBodyForce(kd, ld, cd);
				#endif

				#ifdef VTK_OUTPUT
					if (cd->VtkOutput && ((iter + 1) % cd->VtkModulus) == 0) {
						kd->PdfsActive = src;
						VtkWrite(ld, kd, cd, iter + 1);
					}
				#endif

				#ifdef STATISTICS
					kd->PdfsActive = src;
					KernelStatistics(kd, ld, cd, iter + 1);
				#endif
			}
			#ifdef _OPENMP
				#pragma omp barrier
			#endif
		} // for (int iter = 0; ...
	} // omp parallel

	X_KERNEL_END(kd);

	#ifdef VTK_OUTPUT

	if (cd->VtkOutput) {
		kd->PdfsActive = src;
		VtkWrite(ld, kd, cd, maxIterations);
	}

	#endif

	return;
}

static void KernelEven(LatticeDesc * ld, KernelData * kd, CaseData * cd) // {{{
{
	Assert(ld != NULL);
	Assert(kd != NULL);
	Assert(cd != NULL);

	Assert(cd->Omega > F(0.0));
	Assert(cd->Omega < F(2.0));

	KernelDataAa * kda  = KDA(kd);

	int nX = ld->Dims[0];
	int nY = ld->Dims[1];
	int nZ = ld->Dims[2];

	int * gDims = kd->GlobalDims;

	int oX = kd->Offsets[0];
	int oY = kd->Offsets[1];
	int oZ = kd->Offsets[2];

	int blk[3];
	blk[0] = kda->Blk[0];
	blk[1] = kda->Blk[1];
	blk[2] = kda->Blk[2];

	PdfT omega = cd->Omega;
	PdfT omegaEven = omega;

	PdfT magicParam = F(1.0) / F(12.0);
	PdfT omegaOdd = F(1.0) / (F(0.5) + magicParam / (F(1.0) / omega - F(0.5)));

	const PdfT w_0 = F(1.0) / F( 3.0);
	const PdfT w_1 = F(1.0) / F(18.0);
	const PdfT w_2 = F(1.0) / F(36.0);

	const PdfT w_1_x3 = w_1 * F(3.0);	const PdfT w_1_nine_half = w_1 * F(9.0) / F(2.0);
	const PdfT w_2_x3 = w_2 * F(3.0);	const PdfT w_2_nine_half = w_2 * F(9.0) / F(2.0);


	VPDFT VONE_HALF   = VSET(F(0.5));
	VPDFT VTHREE_HALF = VSET(F(3.0) / F(2.0));

	VPDFT vw_1_indep, vw_2_indep;
	VPDFT vw_0 = VSET(w_0);
	VPDFT vw_1 = VSET(w_1);
	VPDFT vw_2 = VSET(w_2);

	VPDFT vw_1_x3 = VSET(w_1_x3);
	VPDFT vw_2_x3 = VSET(w_2_x3);
	VPDFT vw_1_nine_half = VSET(w_1_nine_half);
	VPDFT vw_2_nine_half = VSET(w_2_nine_half);

	VPDFT vui, vux, vuy, vuz, vdens;

	VPDFT vevenPart, voddPart, vdir_indep_trm;

	VPDFT vomegaEven = VSET(omegaEven);
	VPDFT vomegaOdd  = VSET(omegaOdd);

	VPDFT vpdf_a, vpdf_b;

	// Declare pdf_N, pdf_E, pdf_S, pdf_W, ...
	#define X(name, idx, idxinv, x, y, z)		VPDFT JOIN(vpdf_,name); PdfT * JOIN(ppdf_,name);
		D3Q19_LIST
	#undef X

	PdfT * src = kd->Pdfs[0];

	int nThreads = 1;
	int threadId = 0;

	#ifdef _OPENMP
		nThreads = omp_get_max_threads();
		threadId = omp_get_thread_num();
	#endif

	const int nodesPlane = gDims[1] * gDims[2];
	const int nodesCol   = gDims[2];

	#define I(x, y, z, dir)	P_INDEX_5(gDims, (x), (y), (z), (dir))

// TODO: make inline function out of macros.

	#define IMPLODE(_x, _y, _z)	(nodesPlane * (_x) + nodesCol * (_y) + (_z))
	#define EXPLODE(index, _x, _y, _z)	_x = index / (nodesPlane); _y = (index - nodesPlane * (_x)) / nodesCol; _z = index - nodesPlane * (_x) - nodesCol * (_y);

	int startX = oX;
	int startY = oY;
	int startZ = oZ;

	int indexStart  = IMPLODE(startX, startY, startZ);
	int indexEnd    = IMPLODE(startX + nX - 1, startY + nY - 1, startZ + nZ - 1);

	// How many cells as multiples of VSIZE do we have (rounded up)?
	int nVCells = (indexEnd - indexStart + 1 + VSIZE - 1) / VSIZE;

	int threadStart = nVCells / nThreads * threadId;
	int threadEnd   = nVCells / nThreads * (threadId + 1);

	if (nVCells % nThreads > threadId) {
		threadStart += threadId;
		threadEnd   += threadId + 1;
	}
	else {
		threadStart += nVCells % nThreads;
		threadEnd   += nVCells % nThreads;
	}

	threadStart *= VSIZE;
	threadEnd   *= VSIZE;

	// As threadStart/End is now in the granularity of cells we add the start offset.
	threadStart += indexStart;
	threadEnd   += indexStart;

	EXPLODE(threadStart, startX, startY, startZ);

	#undef EXPLODE
	#undef IMPLODE

	#define X(name, idx, idxinv, _x, _y, _z)	JOIN(ppdf_,name) = &src[I(startX, startY, startZ, idx)];
		D3Q19_LIST
	#undef X

	// printf("e thread %d   idx start: %d end: %d   thread start: %d end: %d\n",
	// 		threadId, indexStart, indexEnd, threadStart, threadEnd);


	for (int i = threadStart; i < threadEnd; i += VSIZE) { // LOOP aa-vec-sl-even

		// Load PDFs of local cell: pdf_N = src[I(x, y, z, D3Q19_N)]; ...
		// #define X(name, idx, idxinv, _x, _y, _z)	JOIN(vpdf_,name) = VLDU(&src[I(x, y, z, idx)]);
		#define X(name, idx, idxinv, _x, _y, _z)	JOIN(vpdf_,name) = VLDU(JOIN(ppdf_,name));
			D3Q19_LIST
		#undef X


		vux = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_E,VADD(vpdf_NE,vpdf_SE)),VADD(vpdf_TE,vpdf_BE)),vpdf_W),vpdf_NW),vpdf_SW),vpdf_TW),vpdf_BW);
		vuy = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_N,VADD(vpdf_NE,vpdf_NW)),VADD(vpdf_TN,vpdf_BN)),vpdf_S),vpdf_SE),vpdf_SW),vpdf_TS),vpdf_BS);
		vuz = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_T,VADD(vpdf_TE,vpdf_TW)),VADD(vpdf_TN,vpdf_TS)),vpdf_B),vpdf_BE),vpdf_BW),vpdf_BN),vpdf_BS);

		vdens = VADD(VADD(VADD(VADD(VADD(VADD(VADD(VADD(VADD(vpdf_C,VADD(vpdf_N,vpdf_E)),VADD(vpdf_S,vpdf_W)),VADD(vpdf_NE,vpdf_SE)),
				VADD(vpdf_SW,vpdf_NW)),VADD(vpdf_T,vpdf_TN)),VADD(vpdf_TE,vpdf_TS)),VADD(vpdf_TW,vpdf_B)),
				VADD(vpdf_BN,vpdf_BE)),VADD(vpdf_BS,vpdf_BW));

		vdir_indep_trm = VSUB(vdens,VMUL(VADD(VADD(VMUL(vux,vux),VMUL(vuy,vuy)),VMUL(vuz,vuz)),VTHREE_HALF));

		VSTU(ppdf_C, VSUB(vpdf_C,VMUL(vomegaEven,VSUB(vpdf_C,VMUL(vw_0,vdir_indep_trm)))));

		vw_1_indep = VMUL(vw_1,vdir_indep_trm);
		vw_2_indep = VMUL(vw_2,vdir_indep_trm);

#if defined(LOOP_1) || defined(LOOP_2)
	#error Loop macros are not allowed to be defined here.
#endif

		#define LOOP_1(_dir1, _dir2, _vel) \
				vui         = _vel; \
				vpdf_a      = JOIN(vpdf_,_dir1); \
				vpdf_b      = JOIN(vpdf_,_dir2); \
				\
				vevenPart = VMUL(vomegaEven, VSUB(VSUB(VMUL(VONE_HALF, VADD(vpdf_a, vpdf_b)), VMUL(vui, VMUL(vui, vw_1_nine_half))), vw_1_indep)); \
				voddPart  = VMUL(vomegaOdd,  VSUB(     VMUL(VONE_HALF, VSUB(vpdf_a, vpdf_b)), VMUL(vui, vw_1_x3))); \
				\
				VSTU(JOIN(ppdf_,_dir2), VSUB(VSUB(vpdf_a, vevenPart), voddPart)); \
				VSTU(JOIN(ppdf_,_dir1), VADD(VSUB(vpdf_b, vevenPart), voddPart));

		#define LOOP_2(_dir1, _dir2, _expr) \
				vui = _expr; \
				vpdf_a 		= JOIN(vpdf_,_dir1); \
				vpdf_b 		= JOIN(vpdf_,_dir2); \
				\
				vevenPart = VMUL(vomegaEven, VSUB(VSUB(VMUL(VONE_HALF, VADD(vpdf_a, vpdf_b)), VMUL(vui, VMUL(vui, vw_2_nine_half))), vw_2_indep)); \
				voddPart  = VMUL(vomegaOdd,  VSUB(     VMUL(VONE_HALF, VSUB(vpdf_a, vpdf_b)), VMUL(vui, vw_2_x3))); \
				\
				VSTU(JOIN(ppdf_,_dir2), VSUB(VSUB(vpdf_a, vevenPart), voddPart)); \
				VSTU(JOIN(ppdf_,_dir1), VADD(VSUB(vpdf_b, vevenPart), voddPart));

		LOOP_1(N, S, vuy);
		LOOP_1(E, W, vux);
		LOOP_1(T, B, vuz);

		LOOP_2(NW, SE, VSUB(vuy, vux));
		LOOP_2(NE, SW, VADD(vuy, vux));
		LOOP_2(TW, BE, VSUB(vuz, vux));
		LOOP_2(TE, BW, VADD(vuz, vux));
		LOOP_2(TS, BN, VSUB(vuz, vuy));
		LOOP_2(TN, BS, VADD(vuz, vuy));

		#undef LOOP_1
		#undef LOOP_2

		#define X(name, idx, idxinv, _x, _y, _z)	JOIN(ppdf_,name) += VSIZE;
			D3Q19_LIST
		#undef X
	}

	#undef I

	return;
} // }}}


static void KernelOddVecSl(LatticeDesc * ld, KernelData * kd, CaseData * cd)  // {{{
{
	Assert(ld != NULL);
	Assert(kd != NULL);
	Assert(cd != NULL);

	Assert(cd->Omega > 0.0);
	Assert(cd->Omega < F(2.0));

	KernelDataAa * kda  = KDA(kd);

	int nX = ld->Dims[0];
	int nY = ld->Dims[1];
	int nZ = ld->Dims[2];

	int * gDims = kd->GlobalDims;

	int oX = kd->Offsets[0];
	int oY = kd->Offsets[1];
	int oZ = kd->Offsets[2];

	int blk[3];
	blk[0] = kda->Blk[0];
	blk[1] = kda->Blk[1];
	blk[2] = kda->Blk[2];

	PdfT omega = cd->Omega;
	PdfT omegaEven = omega;

	PdfT magicParam = F(1.0) / F(12.0);
	PdfT omegaOdd = F(1.0) / (F(0.5) + magicParam / (F(1.0) / omega - F(0.5)));

	const PdfT w_0 = F(1.0) / F( 3.0);
	const PdfT w_1 = F(1.0) / F(18.0);
	const PdfT w_2 = F(1.0) / F(36.0);

	const PdfT w_1_x3 = w_1 * F(3.0);	const PdfT w_1_nine_half = w_1 * F(9.0) / F(2.0);
	const PdfT w_2_x3 = w_2 * F(3.0);	const PdfT w_2_nine_half = w_2 * F(9.0) / F(2.0);

	VPDFT VONE_HALF   = VSET(F(0.5));
	VPDFT VTHREE_HALF = VSET(F(3.0) / F(2.0));

	VPDFT vw_1_indep, vw_2_indep;
	VPDFT vw_0 = VSET(w_0);
	VPDFT vw_1 = VSET(w_1);
	VPDFT vw_2 = VSET(w_2);

	VPDFT vw_1_x3 = VSET(w_1_x3);
	VPDFT vw_2_x3 = VSET(w_2_x3);
	VPDFT vw_1_nine_half = VSET(w_1_nine_half);
	VPDFT vw_2_nine_half = VSET(w_2_nine_half);

	VPDFT vui, vux, vuy, vuz, vdens;

	VPDFT vevenPart, voddPart, vdir_indep_trm;

	VPDFT vomegaEven = VSET(omegaEven);
	VPDFT vomegaOdd  = VSET(omegaOdd);

	VPDFT vpdf_a, vpdf_b;

	// Declare pdf_N, pdf_E, pdf_S, pdf_W, ...
	#define X(name, idx, idxinv, x, y, z) 		VPDFT JOIN(vpdf_,name); PdfT * JOIN(ppdf_,idx);
		D3Q19_LIST
	#undef X

	PdfT * src = kd->Pdfs[0];

	int nThreads = 1;
	int threadId = 0;

	#ifdef _OPENMP
		nThreads = omp_get_max_threads();
		threadId = omp_get_thread_num();
	#endif

	const int nodesPlane = gDims[1] * gDims[2];
	const int nodesCol   = gDims[2];

	#define I(x, y, z, dir)	P_INDEX_5(gDims, (x), (y), (z), (dir))

// TODO: make inline function out of macros.

	#define IMPLODE(_x, _y, _z)	(nodesPlane * (_x) + nodesCol * (_y) + (_z))
	#define EXPLODE(index, _x, _y, _z)	_x = index / (nodesPlane); _y = (index - nodesPlane * (_x)) / nodesCol; _z = index - nodesPlane * (_x) - nodesCol * (_y);

	int startX = oX;
	int startY = oY;
	int startZ = oZ;

	int indexStart = IMPLODE(startX, startY, startZ);
	int indexEnd   = IMPLODE(startX + nX - 1, startY + nY - 1, startZ + nZ - 1);

	// How many multiples of VSIZE cells (rounded up) do we have?
	int nVCells = (indexEnd - indexStart + 1 + VSIZE - 1) / VSIZE;

	int threadStart = nVCells / nThreads * threadId;
	int threadEnd   = nVCells / nThreads * (threadId + 1);

	if (nVCells % nThreads > threadId) {
		threadStart += threadId;
		threadEnd   += threadId + 1;
	}
	else {
		threadStart += nVCells % nThreads;
		threadEnd   += nVCells % nThreads;
	}

	threadStart *= VSIZE;
	threadEnd   *= VSIZE;

	// As threadStart/End is now in the granularity of cells we add the start offset.
	threadStart += indexStart;
	threadEnd   += indexStart;

	EXPLODE(threadStart, startX, startY, startZ);

	#undef EXPLODE
	#undef IMPLODE

	// printf("o thread %d   idx start: %d end: %d   thread start: %d end: %d\n",
	// 		threadId, indexStart, indexEnd, threadStart, threadEnd);

	#define X(name, idx, idxinv, _x, _y, _z)	JOIN(ppdf_,idx) = &src[I(startX + _x, startY + _y, startZ + _z, idx)];
		D3Q19_LIST
	#undef X

#if DEBUG_EXTENDED

	#define X(name, idx, idxinv, x, y, z) 		PdfT * JOIN(ppdf_start_,idx), * JOIN(ppdf_end_,idx);
		D3Q19_LIST
	#undef X

	#define X(name, idx, idxinv, _x, _y, _z)	JOIN(ppdf_start_,idx) = &src[I(startX + _x, startY + _y, startZ + _z, idx)];
		D3Q19_LIST
	#undef X

	#define X(name, idx, idxinv, _x, _y, _z)	JOIN(ppdf_end_,idx) = &src[I(startX + nX - 1 + _x, startY + nY - 1 + _y, startZ  + nZ - 1 + _z, idx)];
		D3Q19_LIST
	#undef X

#if 0
	#define X(name, idx, idxinv, _x, _y, _z)	printf("%2s  ppdf_%d = %p (%d %d %d) (%d %d %d)\n", STRINGIFY(name), idx, JOIN(ppdf_,idx), \
startX , startY , startZ , startX + _x, startY + _y, startZ + _z);
		D3Q19_LIST
	#undef X
#endif

#endif // DEBUG_EXTENDED


	for (int i = threadStart; i < threadEnd; i += VSIZE) { // LOOP aa-vec-sl-odd

#if DEBUG_EXTENDED
		#define X(name, idx, idxinv, _x, _y, _z)	Assert((unsigned long)(JOIN(ppdf_,idx)) >= (unsigned long)(JOIN(ppdf_start_,idx))); Assert((unsigned long)(JOIN(ppdf_,idx)) <= (unsigned long)(JOIN(ppdf_end_,idx)));
			D3Q19_LIST
		#undef X
#endif

		#define X(name, idx, idxinv, _x, _y, _z)	JOIN(vpdf_,name) = VLDU(JOIN(ppdf_,idxinv));
			D3Q19_LIST
		#undef X

		vux = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_E,VADD(vpdf_NE,vpdf_SE)),VADD(vpdf_TE,vpdf_BE)),vpdf_W),vpdf_NW),vpdf_SW),vpdf_TW),vpdf_BW);
		vuy = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_N,VADD(vpdf_NE,vpdf_NW)),VADD(vpdf_TN,vpdf_BN)),vpdf_S),vpdf_SE),vpdf_SW),vpdf_TS),vpdf_BS);
		vuz = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_T,VADD(vpdf_TE,vpdf_TW)),VADD(vpdf_TN,vpdf_TS)),vpdf_B),vpdf_BE),vpdf_BW),vpdf_BN),vpdf_BS);

		vdens = VADD(VADD(VADD(VADD(VADD(VADD(VADD(VADD(VADD(vpdf_C,VADD(vpdf_N,vpdf_E)),VADD(vpdf_S,vpdf_W)),VADD(vpdf_NE,vpdf_SE)),
					 VADD(vpdf_SW,vpdf_NW)),VADD(vpdf_T,vpdf_TN)),VADD(vpdf_TE,vpdf_TS)),VADD(vpdf_TW,vpdf_B)),VADD(vpdf_BN,vpdf_BE)),VADD(vpdf_BS,vpdf_BW));

		vdir_indep_trm = VSUB(vdens,VMUL(VADD(VADD(VMUL(vux,vux),VMUL(vuy,vuy)),VMUL(vuz,vuz)),VTHREE_HALF));

		// ppdf_18 is the pointer to the center pdfs.
		VSTU(ppdf_18, VSUB(vpdf_C,VMUL(vomegaEven,VSUB(vpdf_C,VMUL(vw_0,vdir_indep_trm)))));

		vw_1_indep = VMUL(vw_1,vdir_indep_trm);
		vw_2_indep = VMUL(vw_2,vdir_indep_trm);

#if defined(LOOP_1) || defined(LOOP_2)
	#error Loop macros are not allowed to be defined here.
#endif

		#define LOOP_1(_dir1, _dir2, _idx1, _idx2, _vel) \
				vui         = _vel; \
				vpdf_a      = JOIN(vpdf_,_dir1); \
				vpdf_b      = JOIN(vpdf_,_dir2); \
				\
				vevenPart = VMUL(vomegaEven, VSUB(VSUB(VMUL(VONE_HALF, VADD(vpdf_a, vpdf_b)), VMUL(vui, VMUL(vui, vw_1_nine_half))), vw_1_indep)); \
				voddPart  = VMUL(vomegaOdd,  VSUB(     VMUL(VONE_HALF, VSUB(vpdf_a, vpdf_b)), VMUL(vui, vw_1_x3))); \
				\
				VSTU(JOIN(ppdf_,_idx1), VSUB(VSUB(vpdf_a, vevenPart), voddPart)); \
				VSTU(JOIN(ppdf_,_idx2), VADD(VSUB(vpdf_b, vevenPart), voddPart));

		#define LOOP_2(_dir1, _dir2, _idx1, _idx2, _expr) \
				vui = _expr; \
				vpdf_a 		= JOIN(vpdf_,_dir1); \
				vpdf_b 		= JOIN(vpdf_,_dir2); \
				\
				vevenPart = VMUL(vomegaEven, VSUB(VSUB(VMUL(VONE_HALF, VADD(vpdf_a, vpdf_b)), VMUL(vui, VMUL(vui, vw_2_nine_half))), vw_2_indep)); \
				voddPart  = VMUL(vomegaOdd,  VSUB(     VMUL(VONE_HALF, VSUB(vpdf_a, vpdf_b)), VMUL(vui, vw_2_x3))); \
				\
				VSTU(JOIN(ppdf_,_idx1), VSUB(VSUB(vpdf_a, vevenPart), voddPart)); \
				VSTU(JOIN(ppdf_,_idx2), VADD(VSUB(vpdf_b, vevenPart), voddPart));


		LOOP_1(N, S, D3Q19_N, D3Q19_S, vuy);
		LOOP_1(E, W, D3Q19_E, D3Q19_W, vux);
		LOOP_1(T, B, D3Q19_T, D3Q19_B, vuz);

		LOOP_2(NW, SE, D3Q19_NW, D3Q19_SE, VSUB(vuy, vux));
		LOOP_2(NE, SW, D3Q19_NE, D3Q19_SW, VADD(vuy, vux));
		LOOP_2(TW, BE, D3Q19_TW, D3Q19_BE, VSUB(vuz, vux));
		LOOP_2(TE, BW, D3Q19_TE, D3Q19_BW, VADD(vuz, vux));
		LOOP_2(TS, BN, D3Q19_TS, D3Q19_BN, VSUB(vuz, vuy));
		LOOP_2(TN, BS, D3Q19_TN, D3Q19_BS, VADD(vuz, vuy));

		#define X(name, idx, idxinv, _x, _y, _z)	JOIN(ppdf_,idx) += VSIZE;
			D3Q19_LIST
		#undef X
	}

	#undef I

	return;

}  // }}}
