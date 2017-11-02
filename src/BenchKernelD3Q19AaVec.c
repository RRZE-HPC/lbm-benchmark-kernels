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
static void KernelOdd( LatticeDesc * ld, KernelData * kd, CaseData * cd);

#if 0 // {{{
void DumpPdfs(LatticeDesc * ld, KernelData * kd, int zStart, int zStop, int iter, const char * prefix)
{
	return;

	int * lDims = ld->Dims;
	int * gDims = kd->GlobalDims;

	int nX = gDims[0];
	int nY = gDims[1];
	int nZ = gDims[2];

	PdfT pdfs[N_D3Q19];

	int localZStart = zStart;
	int localZStop  = zStop;

	if (localZStart == -1) localZStart = 0;
	if (localZStop  == -1) localZStop  = gDims[2] - 1;

	printf("D iter: %d\n", iter);

	for (int dir = 0; dir < 19; ++dir) {
		for (int z = localZStop; z >= localZStart; --z) {
			printf("D [%2d][%2d][%s] plane % 2d\n", iter, dir, prefix, z);

			for(int y = 0; y < nY; ++y) {
				printf("D [%2d][%2d][%s] %2d  ", iter, dir, prefix, y);

				for(int x = 0; x < nX; ++x) {

					if (1) { // ld->Lattice[L_INDEX_4(ld->Dims, x, y, z)] != LAT_CELL_OBSTACLE) {

						#define I(x, y, z, dir)	P_INDEX_5(gDims, (x), (y), (z), (dir))
						pdfs[dir] = kd->PdfsActive[I(x, y, z, dir)];
						#undef I
//						kd->GetNode(kd, x, y, z, pdfs);
					}
					else {
						pdfs[dir] = -1.0;
					}

					printf("%.16e ", pdfs[dir]);
				}

				printf("\n");
			}
		}
	}
}
#endif  // }}}

void FNAME(D3Q19AaVecKernel)(LatticeDesc * ld, KernelData * kd, CaseData * cd)
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

	for (int iter = 0; iter < maxIterations; iter += 2) {

		// --------------------------------------------------------------------
		// even time step
		// --------------------------------------------------------------------

		X_LIKWID_START("aa-vec-even");

		#pragma omp parallel
		{
			KernelEven(ld, kd, cd);
		}

		X_LIKWID_STOP("aa-vec-even");

		// Fixup bounce back PDFs.
		#ifdef _OPENMP
			#pragma omp parallel for default(none) \
					shared(kd, src)
		#endif
		for (int i = 0; i < kd->nBounceBackPdfs; ++i) {
			src[kd->BounceBackPdfsSrc[i]] = src[kd->BounceBackPdfsDst[i]];
		}

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

		// --------------------------------------------------------------------
		// odd time step
		// --------------------------------------------------------------------

		X_LIKWID_START("aa-vec-odd");

		#pragma omp parallel
		{
			KernelOdd(ld, kd, cd);
		}

		// Stop counters before bounce back. Else computing loop balance will
		// be incorrect.

		X_LIKWID_STOP("aa-vec-odd");

		// Fixup bounce back PDFs.
		#ifdef _OPENMP
			#pragma omp parallel for default(none) \
					shared(kd, src)
		#endif
		for (int i = 0; i < kd->nBounceBackPdfs; ++i) {
			src[kd->BounceBackPdfsDst[i]] = src[kd->BounceBackPdfsSrc[i]];
		}

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
		#endif // }}}


	} // for (int iter = 0; ...

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

	Assert(cd->Omega > 0.0);
	Assert(cd->Omega < 2.0);

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

	PdfT magicParam = 1.0 / 12.0;
	PdfT omegaOdd = 1.0 / (0.5 + magicParam / (1.0 / omega - 0.5));

	const PdfT w_0 = 1.0 /  3.0;
	const PdfT w_1 = 1.0 / 18.0;
	const PdfT w_2 = 1.0 / 36.0;

	const PdfT w_1_x3 = w_1 * 3.0;	const PdfT w_1_nine_half = w_1 * 9.0 / 2.0;
	const PdfT w_2_x3 = w_2 * 3.0;	const PdfT w_2_nine_half = w_2 * 9.0 / 2.0;


	VPDFT VONE_HALF   = VSET(0.5);
	VPDFT VTHREE_HALF = VSET(3.0 / 2.0);

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

	// Declare pdf_N, pdf_E, pdf_S, pdf_W, ...
	#define X(name, idx, idxinv, x, y, z)	VPDFT JOIN(vpdf_,name);
		D3Q19_LIST
	#undef X

	PdfT * src = kd->Pdfs[0];

	int nThreads = 1;
	int threadId = 0;

	#ifdef _OPENMP
		nThreads = omp_get_max_threads();
		threadId = omp_get_thread_num();
	#endif

	// TODO: Currently only a 1-D decomposition is applied. For achritectures
	//       with a lot of cores we want at least 2-D.

	int threadStartX = nX / nThreads * threadId;
	int threadEndX   = nX / nThreads * (threadId + 1);

	if (nX % nThreads > 0) {
		if (nX % nThreads > threadId) {
			threadStartX += threadId;
			threadEndX   += threadId + 1;
		}
		else {
			threadStartX += nX % nThreads;
			threadEndX   += nX % nThreads;
		}
	}

	AssertMsg((blk[2] % VSIZE == 0) || blk[2] >= nZ, "Blocking in z direction must be a multiple of VSIZE = %d or larger than z dimension.", VSIZE);

	for (int bX = oX + threadStartX; bX < threadEndX + oX; bX += blk[0]) {
	for (int bY = oY; bY < nY + oY; bY += blk[1]) {
	for (int bZ = oZ; bZ < nZ + oZ; bZ += blk[2]) {

		int eX = MIN(bX + blk[0], threadEndX + oX);
		int eY = MIN(bY + blk[1], nY + oY);
		int eZ = MIN(bZ + blk[2], nZ + oZ);

		for (int x = bX; x < eX; x += 1) {
		for (int y = bY; y < eY; y += 1) {
		for (int z = bZ; z < eZ; z += VSIZE) {

			#define I(x, y, z, dir)	P_INDEX_5(gDims, (x), (y), (z), (dir))

			// Load PDFs of local cell: pdf_N = src[I(x, y, z, D3Q19_N)]; ...
			#define X(name, idx, idxinv, _x, _y, _z)	JOIN(vpdf_,name) = VLDU(&src[I(x, y, z, idx)]);
				D3Q19_LIST
			#undef X


			vux = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_E,VADD(vpdf_NE,vpdf_SE)),VADD(vpdf_TE,vpdf_BE)),vpdf_W),vpdf_NW),vpdf_SW),vpdf_TW),vpdf_BW);
			vuy = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_N,VADD(vpdf_NE,vpdf_NW)),VADD(vpdf_TN,vpdf_BN)),vpdf_S),vpdf_SE),vpdf_SW),vpdf_TS),vpdf_BS);
			vuz = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_T,VADD(vpdf_TE,vpdf_TW)),VADD(vpdf_TN,vpdf_TS)),vpdf_B),vpdf_BE),vpdf_BW),vpdf_BN),vpdf_BS);

			vdens = VADD(VADD(VADD(VADD(VADD(VADD(VADD(VADD(VADD(vpdf_C,VADD(vpdf_N,vpdf_E)),VADD(vpdf_S,vpdf_W)),VADD(vpdf_NE,vpdf_SE)),
					VADD(vpdf_SW,vpdf_NW)),VADD(vpdf_T,vpdf_TN)),VADD(vpdf_TE,vpdf_TS)),VADD(vpdf_TW,vpdf_B)),
					VADD(vpdf_BN,vpdf_BE)),VADD(vpdf_BS,vpdf_BW));

			vdir_indep_trm = VSUB(vdens,VMUL(VADD(VADD(VMUL(vux,vux),VMUL(vuy,vuy)),VMUL(vuz,vuz)),VTHREE_HALF));

			VSTU(&src[I(x, y, z, D3Q19_C)],VSUB(vpdf_C,VMUL(vomegaEven,VSUB(vpdf_C,VMUL(vw_0,vdir_indep_trm)))));

			vw_1_indep = VMUL(vw_1,vdir_indep_trm);

			vui = vuy;
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_N,vpdf_S)),VMUL(vui,VMUL(vui,vw_1_nine_half))),vw_1_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_N,vpdf_S)),VMUL(vui,vw_1_x3)));
			VSTU(&src[I(x, y, z, D3Q19_S)],VSUB(VSUB(vpdf_N,vevenPart),voddPart));
			VSTU(&src[I(x, y, z, D3Q19_N)],VADD(VSUB(vpdf_S,vevenPart),voddPart));

			vui = vux;
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_E,vpdf_W)),VMUL(vui,VMUL(vui,vw_1_nine_half))),vw_1_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_E,vpdf_W)),VMUL(vui,vw_1_x3)));
			VSTU(&src[I(x, y, z, D3Q19_W)],VSUB(VSUB(vpdf_E,vevenPart),voddPart));
			VSTU(&src[I(x, y, z, D3Q19_E)],VADD(VSUB(vpdf_W,vevenPart),voddPart));

			vui = vuz;
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_T,vpdf_B)),VMUL(vui,VMUL(vui,vw_1_nine_half))),vw_1_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_T,vpdf_B)),VMUL(vui,vw_1_x3)));
			VSTU(&src[I(x, y, z, D3Q19_B)],VSUB(VSUB(vpdf_T,vevenPart),voddPart));
			VSTU(&src[I(x, y, z, D3Q19_T)],VADD(VSUB(vpdf_B,vevenPart),voddPart));

			vw_2_indep = VMUL(vw_2,vdir_indep_trm);

			vui = VSUB(vuy,vux);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_NW,vpdf_SE)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_NW,vpdf_SE)),VMUL(vui,vw_2_x3)));
			VSTU(&src[I(x, y, z, D3Q19_SE)],VSUB(VSUB(vpdf_NW,vevenPart),voddPart));
			VSTU(&src[I(x, y, z, D3Q19_NW)],VADD(VSUB(vpdf_SE,vevenPart),voddPart));

			vui = VADD(vux,vuy);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_NE,vpdf_SW)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_NE,vpdf_SW)),VMUL(vui,vw_2_x3)));
			VSTU(&src[I(x, y, z, D3Q19_SW)],VSUB(VSUB(vpdf_NE,vevenPart),voddPart));
			VSTU(&src[I(x, y, z, D3Q19_NE)],VADD(VSUB(vpdf_SW,vevenPart),voddPart));

			vui = VSUB(vuz,vux);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TW,vpdf_BE)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TW,vpdf_BE)),VMUL(vui,vw_2_x3)));
			VSTU(&src[I(x, y, z, D3Q19_BE)],VSUB(VSUB(vpdf_TW,vevenPart),voddPart));
			VSTU(&src[I(x, y, z, D3Q19_TW)],VADD(VSUB(vpdf_BE,vevenPart),voddPart));

			vui = VADD(vux,vuz);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TE,vpdf_BW)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TE,vpdf_BW)),VMUL(vui,vw_2_x3)));
			VSTU(&src[I(x, y, z, D3Q19_BW)],VSUB(VSUB(vpdf_TE,vevenPart),voddPart));
			VSTU(&src[I(x, y, z, D3Q19_TE)],VADD(VSUB(vpdf_BW,vevenPart),voddPart));

			vui = VSUB(vuz,vuy);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TS,vpdf_BN)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TS,vpdf_BN)),VMUL(vui,vw_2_x3)));
			VSTU(&src[I(x, y, z, D3Q19_BN)],VSUB(VSUB(vpdf_TS,vevenPart),voddPart));
			VSTU(&src[I(x, y, z, D3Q19_TS)],VADD(VSUB(vpdf_BN,vevenPart),voddPart));

			vui = VADD(vuy,vuz);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TN,vpdf_BS)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TN,vpdf_BS)),VMUL(vui,vw_2_x3)));
			VSTU(&src[I(x, y, z, D3Q19_BS)],VSUB(VSUB(vpdf_TN,vevenPart),voddPart));
			VSTU(&src[I(x, y, z, D3Q19_TN)],VADD(VSUB(vpdf_BS,vevenPart),voddPart));

			#undef I
		} } } // x, y, z
	} } } // blocked x, y, z



	return;
} // }}}


static void KernelOdd(LatticeDesc * ld, KernelData * kd, CaseData * cd)  // {{{
{
	Assert(ld != NULL);
	Assert(kd != NULL);
	Assert(cd != NULL);

	Assert(cd->Omega > 0.0);
	Assert(cd->Omega < 2.0);

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

	PdfT magicParam = 1.0 / 12.0;
	PdfT omegaOdd = 1.0 / (0.5 + magicParam / (1.0 / omega - 0.5));

	const PdfT w_0 = 1.0 /  3.0;
	const PdfT w_1 = 1.0 / 18.0;
	const PdfT w_2 = 1.0 / 36.0;

	const PdfT w_1_x3 = w_1 * 3.0;	const PdfT w_1_nine_half = w_1 * 9.0 / 2.0;
	const PdfT w_2_x3 = w_2 * 3.0;	const PdfT w_2_nine_half = w_2 * 9.0 / 2.0;

	VPDFT VONE_HALF   = VSET(0.5);
	VPDFT VTHREE_HALF = VSET(3.0 / 2.0);

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

	// Declare pdf_N, pdf_E, pdf_S, pdf_W, ...
	#define X(name, idx, idxinv, x, y, z) 		VPDFT JOIN(vpdf_,name);
		D3Q19_LIST
	#undef X

	PdfT * src = kd->Pdfs[0];

	int threadId = 0;
	int nThreads = 1;

	#ifdef _OPENMP
		threadId = omp_get_thread_num();
		nThreads = omp_get_max_threads();
	#endif

	// TODO: Currently only a 1-D decomposition is applied. For achritectures
	//       with a lot of cores we want at least 2-D.
	int threadStartX = nX / nThreads * threadId;
	int threadEndX   = nX / nThreads * (threadId + 1);

	if (nX % nThreads > 0) {
		if (nX % nThreads > threadId) {
			threadStartX += threadId;
			threadEndX   += threadId + 1;
		}
		else {
			threadStartX += nX % nThreads;
			threadEndX   += nX % nThreads;
		}
	}

	AssertMsg((blk[2] % VSIZE == 0) || blk[2] >= nZ, "Blocking in z direction must be a multiple of VSIZE = %d or larger than z dimension.", VSIZE);

	for (int bX = oX + threadStartX; bX < threadEndX + oX; bX += blk[0]) {
	for (int bY = oY; bY < nY + oY; bY += blk[1]) {
	for (int bZ = oZ; bZ < nZ + oZ; bZ += blk[2]) {

		int eX = MIN(bX + blk[0], threadEndX + oX);
		int eY = MIN(bY + blk[1], nY + oY);
		int eZ = MIN(bZ + blk[2], nZ + oZ);

		for (int x = bX; x < eX; ++x) {
		for (int y = bY; y < eY; ++y) {
		for (int z = bZ; z < eZ; z += VSIZE) {

			#define I(x, y, z, dir)	P_INDEX_5(gDims, (x), (y), (z), (dir))


			#define X(name, idx, idxinv, _x, _y, _z)	JOIN(vpdf_,name) = VLDU(&src[I(x - _x, y - _y, z - _z, idxinv)]);
			D3Q19_LIST
			#undef X

			vux = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_E,VADD(vpdf_NE,vpdf_SE)),VADD(vpdf_TE,vpdf_BE)),vpdf_W),vpdf_NW),vpdf_SW),vpdf_TW),vpdf_BW);
			vuy = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_N,VADD(vpdf_NE,vpdf_NW)),VADD(vpdf_TN,vpdf_BN)),vpdf_S),vpdf_SE),vpdf_SW),vpdf_TS),vpdf_BS);
			vuz = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_T,VADD(vpdf_TE,vpdf_TW)),VADD(vpdf_TN,vpdf_TS)),vpdf_B),vpdf_BE),vpdf_BW),vpdf_BN),vpdf_BS);

			vdens = VADD(VADD(VADD(VADD(VADD(VADD(VADD(VADD(VADD(vpdf_C,VADD(vpdf_N,vpdf_E)),VADD(vpdf_S,vpdf_W)),VADD(vpdf_NE,vpdf_SE)),
						 VADD(vpdf_SW,vpdf_NW)),VADD(vpdf_T,vpdf_TN)),VADD(vpdf_TE,vpdf_TS)),VADD(vpdf_TW,vpdf_B)),VADD(vpdf_BN,vpdf_BE)),VADD(vpdf_BS,vpdf_BW));

			vdir_indep_trm = VSUB(vdens,VMUL(VADD(VADD(VMUL(vux,vux),VMUL(vuy,vuy)),VMUL(vuz,vuz)),VTHREE_HALF));

			VSTU(&src[I(x, y, z, D3Q19_C)],VSUB(vpdf_C,VMUL(vomegaEven,VSUB(vpdf_C,VMUL(vw_0,vdir_indep_trm)))));

			vw_1_indep = VMUL(vw_1,vdir_indep_trm);

			vui = vuy;
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_N,vpdf_S)),VMUL(vui,VMUL(vui,vw_1_nine_half))),vw_1_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_N,vpdf_S)),VMUL(vui,vw_1_x3)));
			VSTU(&src[I(x, y + 1,  z, D3Q19_N)], VSUB(VSUB(vpdf_N,vevenPart),voddPart));
			VSTU(&src[I(x, y - 1,  z, D3Q19_S)], VADD(VSUB(vpdf_S,vevenPart),voddPart));

			vui = vux;
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_E,vpdf_W)),VMUL(vui,VMUL(vui,vw_1_nine_half))),vw_1_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_E,vpdf_W)),VMUL(vui,vw_1_x3)));
			VSTU(&src[I(x + 1, y, z, D3Q19_E)], VSUB(VSUB(vpdf_E,vevenPart),voddPart));
			VSTU(&src[I(x - 1, y, z, D3Q19_W)], VADD(VSUB(vpdf_W,vevenPart),voddPart));

			vui = vuz;
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_T,vpdf_B)),VMUL(vui,VMUL(vui,vw_1_nine_half))),vw_1_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_T,vpdf_B)),VMUL(vui,vw_1_x3)));
			VSTU(&src[I(x, y, z + 1, D3Q19_T)], VSUB(VSUB(vpdf_T,vevenPart),voddPart));
			VSTU(&src[I(x, y, z - 1, D3Q19_B)], VADD(VSUB(vpdf_B,vevenPart),voddPart));

			vw_2_indep = VMUL(vw_2,vdir_indep_trm);

			vui = VSUB(vuy,vux);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_NW,vpdf_SE)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_NW,vpdf_SE)),VMUL(vui,vw_2_x3)));
			VSTU(&src[I(x - 1, y + 1, z, D3Q19_NW)], VSUB(VSUB(vpdf_NW,vevenPart),voddPart));
			VSTU(&src[I(x + 1, y - 1, z, D3Q19_SE)], VADD(VSUB(vpdf_SE,vevenPart),voddPart));

			vui = VADD(vux,vuy);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_NE,vpdf_SW)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_NE,vpdf_SW)),VMUL(vui,vw_2_x3)));
			VSTU(&src[I(x + 1, y + 1, z, D3Q19_NE)], VSUB(VSUB(vpdf_NE,vevenPart),voddPart));
			VSTU(&src[I(x - 1, y - 1, z, D3Q19_SW)], VADD(VSUB(vpdf_SW,vevenPart),voddPart));

			vui = VSUB(vuz,vux);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TW,vpdf_BE)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TW,vpdf_BE)),VMUL(vui,vw_2_x3)));
			VSTU(&src[I(x - 1, y, z + 1, D3Q19_TW)], VSUB(VSUB(vpdf_TW,vevenPart),voddPart));
			VSTU(&src[I(x + 1, y, z - 1, D3Q19_BE)], VADD(VSUB(vpdf_BE,vevenPart),voddPart));

			vui = VADD(vux,vuz);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TE,vpdf_BW)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TE,vpdf_BW)),VMUL(vui,vw_2_x3)));
			VSTU(&src[I(x + 1, y, z + 1, D3Q19_TE)], VSUB(VSUB(vpdf_TE,vevenPart),voddPart));
			VSTU(&src[I(x - 1, y, z - 1, D3Q19_BW)], VADD(VSUB(vpdf_BW,vevenPart),voddPart));

			vui = VSUB(vuz,vuy);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TS,vpdf_BN)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TS,vpdf_BN)),VMUL(vui,vw_2_x3)));
			VSTU(&src[I(x, y - 1, z + 1, D3Q19_TS)], VSUB(VSUB(vpdf_TS,vevenPart),voddPart));
			VSTU(&src[I(x, y + 1, z - 1, D3Q19_BN)], VADD(VSUB(vpdf_BN,vevenPart),voddPart));

			vui = VADD(vuy,vuz);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TN,vpdf_BS)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TN,vpdf_BS)),VMUL(vui,vw_2_x3)));
			VSTU(&src[I(x, y + 1, z + 1, D3Q19_TN)], VSUB(VSUB(vpdf_TN,vevenPart),voddPart));
			VSTU(&src[I(x, y - 1, z - 1, D3Q19_BS)], VADD(VSUB(vpdf_BS,vevenPart),voddPart));

			#undef I
		} } } // x, y, z
	} } } // blocked x, y, z

	return;

}  // }}}
