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
#include "BenchKernelD3Q19ListAaPvGatherAoSoACommon.h"

#include "Base.h"
#include "Memory.h"
#include "Vtk.h"
#include "Vector.h"

#include <inttypes.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

//enable software prefetchting for vectorized gather/scatter loop in odd kernel
#ifndef SOFTWARE_PREFETCH_LOOKAHEAD_L2
#define SOFTWARE_PREFETCH_LOOKAHEAD_L2 (0) //prefetchting X SIMD widths ahead
#endif

#ifndef SOFTWARE_PREFETCH_LOOKAHEAD_L1
#define SOFTWARE_PREFETCH_LOOKAHEAD_L1 (0) //prefetchting X SIMD widths ahead
#endif

static void KernelEven(LatticeDesc * ld, KernelData * kernelData, CaseData * cd, int * threadIndices);
static void KernelOdd( LatticeDesc * ld, KernelData * kernelData, CaseData * cd, int * threadIndices);

void FNAME(D3Q19ListAaPvGatherAoSoAKernel)(LatticeDesc * ld, KernelData * kernelData, CaseData * cd)
{

	Assert(ld != NULL);
	Assert(kernelData != NULL);
	Assert(cd != NULL);

	Assert(cd->Omega > 0.0);
	Assert(cd->Omega < 2.0);

#if defined(VTK_OUTPUT) || defined(STATISTICS) || defined(VERIFICATION)
	KernelData * kd = (KernelData *)kernelData;
#endif
	KernelDataList * kdl = KDL(kernelData);

	int maxIterations = cd->MaxIterations;
	int nFluid = kdl->nFluid;

	printf("\n");
#if (SOFTWARE_PREFETCH_LOOKAHEAD_L2 > 0) || (SOFTWARE_PREFETCH_LOOKAHEAD_L1 > 0)
	printf("# Software prefetching enabled:\n");
	printf("#   Gather/Scatter prefetch lookahead L2: \t%d\n", SOFTWARE_PREFETCH_LOOKAHEAD_L2);
	printf("#   Gather/Scatter prefetch lookahead L1: \t%d\n", SOFTWARE_PREFETCH_LOOKAHEAD_L1);
#else
	printf("# Software prefetching disabled.\n");
#endif
	printf("\n");

	int nThreads = 1;
#ifdef _OPENMP
	nThreads = omp_get_max_threads();
#endif

	int * threadIndices = (int *)malloc(sizeof(int) * (nThreads + 1));
	for (int i = 0; i < nThreads; ++i) {
		threadIndices[i] = i * (nFluid / nThreads) + MinI(i, nFluid % nThreads);
	}

	threadIndices[nThreads] = nFluid;

#ifdef VTK_OUTPUT
	if (cd->VtkOutput) {
		kd->PdfsActive = kd->Pdfs[0];
		VtkWrite(ld, kd, cd, -1);
	}
#endif

#ifdef STATISTICS
	kd->PdfsActive = kd->Pdfs[0];
	KernelStatistics(kd, ld, cd, 0);
#endif

	LIKWID_MARKER_INIT;

	X_KERNEL_START(kernelData);

	// TODO: outer openmp parallel

	LIKWID_MARKER_START("OuterLoop");
	for(int iter = 0; iter < maxIterations; iter += 2) {

		// even time step

#ifdef _OPENMP
#pragma omp parallel default(none) shared(ld, kernelData, cd, threadIndices)
#endif
		{
			KernelEven(ld, kernelData, cd, threadIndices);
		}


#ifdef VERIFICATION
		kdl->Iteration = iter;
		kd->PdfsActive = kd->Pdfs[0];
		KernelAddBodyForce(kd, ld, cd);
#endif

		// odd time step

#ifdef _OPENMP
#pragma omp parallel default(none) shared(ld, kernelData, cd, threadIndices)
#endif
		{
			KernelOdd(ld, kernelData, cd, threadIndices);
		}


#ifdef VERIFICATION
		kdl->Iteration = iter + 1;
		kd->PdfsActive = kd->Pdfs[0];
		KernelAddBodyForce(kd, ld, cd);
#endif

#ifdef VTK_OUTPUT
		if (cd->VtkOutput && (iter % cd->VtkModulus) == 0) {
			kdl->Iteration = iter + 1;
			kd->PdfsActive = kd->Pdfs[0];
			VtkWrite(ld, kd, cd, iter);
		}
#endif

#ifdef STATISTICS
		kdl->Iteration = iter + 1;
		kd->PdfsActive = kd->Pdfs[0];
		KernelStatistics(kd, ld, cd, iter);
#endif

	} // for (int iter = 0; ...
	LIKWID_MARKER_STOP("OuterLoop");

	X_KERNEL_END(kernelData);

#ifdef VTK_OUTPUT
	if (cd->VtkOutput) {
		kd->PdfsActive = kd->Pdfs[0];
		VtkWrite(ld, kd, cd, maxIterations);
	}
#endif

#ifdef STATISTICS
	kd->PdfsActive = kd->Pdfs[0];
	KernelStatistics(kd, ld, cd, maxIterations);
#endif

	LIKWID_MARKER_CLOSE;
	free(threadIndices);

	return;
}

static void KernelEven(LatticeDesc * ld, KernelData * kernelData, CaseData * cd, int * threadIndices)
{
	Assert(ld != NULL);
	Assert(kernelData != NULL);
	Assert(cd != NULL);

	Assert(cd->Omega > 0.0);
	Assert(cd->Omega < 2.0);

	KernelData * kd = (KernelData *)kernelData;
	KernelDataList * kdl = KDL(kernelData);

	PdfT omega = cd->Omega;
	PdfT omegaEven = omega;

	PdfT magicParam = 1.0 / 12.0;
	PdfT omegaOdd = 1.0 / (0.5 + magicParam / (1.0 / omega - 0.5));

	PdfT evenPart = 0.0;
	PdfT oddPart = 0.0;
	PdfT dir_indep_trm = 0.0;

	const PdfT w_0 = 1.0 /  3.0;
	const PdfT w_1 = 1.0 / 18.0;
	const PdfT w_2 = 1.0 / 36.0;

	const PdfT w_1_x3 = w_1 * 3.0;	const PdfT w_1_nine_half = w_1 * 9.0 / 2.0;	PdfT w_1_indep = 0.0;
	const PdfT w_2_x3 = w_2 * 3.0;	const PdfT w_2_nine_half = w_2 * 9.0 / 2.0;	PdfT w_2_indep = 0.0;

	PdfT ux, uy, uz, ui;
	PdfT dens;

	VPDFT VONE_HALF = VSET(0.5);
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
	#define X(name, idx, idxinv, x, y, z) \
		PdfT JOIN(pdf_,name); \
		PdfT * JOIN(ppdf_,name); \
		VPDFT JOIN(vpdf_,name);
		D3Q19_LIST
	#undef X

	PdfT * src = kd->Pdfs[0];

	int nCells = kdl->nCells;

	int threadId = 0;
#ifdef _OPENMP
	threadId =  omp_get_thread_num();
#endif


	int indexStart    = threadIndices[threadId];
	int nFluidThread  = threadIndices[threadId + 1] - threadIndices[threadId];
	int indexStop     = indexStart + nFluidThread;

	int indexStartVec = ((indexStart + VSIZE - 1) / VSIZE) * VSIZE;
	int nFluidVec     = (indexStop / VSIZE) * VSIZE - indexStartVec;
	int indexStopVec  = indexStartVec + nFluidVec;

	#define I(index, dir)	P_INDEX_3((nCells), (index), (dir))

	#define X(name, idx, idxinv, _x, _y, _z)	JOIN(ppdf_,name) = &(src[I(indexStart, idx)]);
			D3Q19_LIST
	#undef X

	for (int index = indexStart; index < indexStartVec; ++index) {

		#define X(name, idx, idxinv, _x, _y, _z)	JOIN(pdf_,name) = *(JOIN(ppdf_,name));
			D3Q19_LIST
		#undef X

		ux = pdf_E + pdf_NE + pdf_SE + pdf_TE + pdf_BE -
			pdf_W - pdf_NW - pdf_SW - pdf_TW - pdf_BW;
		uy = pdf_N + pdf_NE + pdf_NW + pdf_TN + pdf_BN -
			pdf_S - pdf_SE - pdf_SW - pdf_TS - pdf_BS;
		uz = pdf_T + pdf_TE + pdf_TW + pdf_TN + pdf_TS -
			pdf_B - pdf_BE - pdf_BW - pdf_BN - pdf_BS;

		dens = pdf_C +
			pdf_N  + pdf_E  + pdf_S  + pdf_W  +
			pdf_NE + pdf_SE + pdf_SW + pdf_NW +
			pdf_T  + pdf_TN + pdf_TE + pdf_TS + pdf_TW +
			pdf_B  + pdf_BN + pdf_BE + pdf_BS + pdf_BW;

		dir_indep_trm = dens - (ux * ux + uy * uy + uz * uz)*3.0/2.0;

		// direction: w_0
		*ppdf_C  = pdf_C - omegaEven*(pdf_C - w_0*dir_indep_trm);

		// direction: w_1
		w_1_indep = w_1*dir_indep_trm;

		#define COLLIDE_AA_S(tmpUi, dir1, dir2) \
			ui = tmpUi; \
			evenPart = omegaEven * (0.5 * (JOIN(pdf_,dir1) + JOIN(pdf_,dir2)) - ui * ui * w_1_nine_half - w_1_indep); \
			oddPart  = omegaOdd  * (0.5 * (JOIN(pdf_,dir1) - JOIN(pdf_,dir2)) - ui * w_1_x3); \
			*(JOIN(ppdf_,dir2))  = JOIN(pdf_,dir1) - evenPart - oddPart; \
			*(JOIN(ppdf_,dir1))  = JOIN(pdf_,dir2) - evenPart + oddPart;

		COLLIDE_AA_S(uy, N, S)
		COLLIDE_AA_S(ux, E, W)
		COLLIDE_AA_S(uz, T, B)

		#undef COLLIDE_AA_S

		// direction: w_2
		w_2_indep = w_2*dir_indep_trm;

		#define COLLIDE_UA_S(tmpUi, dir1, dir2) \
			ui = tmpUi; \
			evenPart = omegaEven * (0.5 * (JOIN(pdf_,dir1) + JOIN(pdf_,dir2)) - ui * ui * w_2_nine_half - w_2_indep); \
			oddPart  = omegaOdd  * (0.5 * (JOIN(pdf_,dir1) - JOIN(pdf_,dir2)) - ui * w_2_x3); \
			*(JOIN(ppdf_,dir2)) = JOIN(pdf_,dir1) - evenPart - oddPart; \
			*(JOIN(ppdf_,dir1)) = JOIN(pdf_,dir2) - evenPart + oddPart;

		COLLIDE_UA_S((-ux + uy), NW, SE)
		COLLIDE_UA_S(( ux + uy), NE, SW)
		COLLIDE_UA_S((-ux + uz), TW, BE)
		COLLIDE_UA_S(( ux + uz), TE, BW)
		COLLIDE_UA_S((-uy + uz), TS, BN)
		COLLIDE_UA_S(( uy + uz), TN, BS)

		#undef COLLIDE_UA_S

		#define X(name, idx, idxinv, _x, _y, _z)	JOIN(ppdf_,name)++;
			D3Q19_LIST
		#undef X
	}

	#define X(name, idx, idxinv, _x, _y, _z)	JOIN(ppdf_,name) = &(src[I(indexStartVec, idx)]);
			D3Q19_LIST
	#undef X

	for (int index = indexStartVec; index < indexStopVec; index += VSIZE) {

		#if (SOFTWARE_PREFETCH_LOOKAHEAD_L2 > 0)
			#define X(name, idx, idxinv, _x, _y, _z) _mm_prefetch((char const *)(JOIN(ppdf_,name) + SOFTWARE_PREFETCH_LOOKAHEAD_L2 * VSIZE * N_D3Q19), _MM_HINT_T1);
				D3Q19_LIST
			#undef X
		#endif

		#if (SOFTWARE_PREFETCH_LOOKAHEAD_L1 > 0)
			#define X(name, idx, idxinv, _x, _y, _z) _mm_prefetch((char const *)(JOIN(ppdf_,name) + SOFTWARE_PREFETCH_LOOKAHEAD_L1 * VSIZE * N_D3Q19), _MM_HINT_T0);
				D3Q19_LIST
			#undef X
		#endif

		#define X(name, idx, idxinv, _x, _y, _z)	JOIN(vpdf_,name) = VLDU(JOIN(ppdf_,name));
			D3Q19_LIST
		#undef X

		//vux = vpdf_E + vpdf_NE + vpdf_SE + vpdf_TE + vpdf_BE -
		//	     vpdf_W - vpdf_NW - vpdf_SW - vpdf_TW - vpdf_BW;
		vux = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_E,VADD(vpdf_NE,vpdf_SE)),VADD(vpdf_TE,vpdf_BE)),vpdf_W),vpdf_NW),vpdf_SW),vpdf_TW),vpdf_BW);
		//vuy = vpdf_N + vpdf_NE + vpdf_NW + vpdf_TN + vpdf_BN -
		//	     vpdf_S - vpdf_SE - vpdf_SW - vpdf_TS - vpdf_BS;
		vuy = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_N,VADD(vpdf_NE,vpdf_NW)),VADD(vpdf_TN,vpdf_BN)),vpdf_S),vpdf_SE),vpdf_SW),vpdf_TS),vpdf_BS);
		//vuz = vpdf_T + vpdf_TE + vpdf_TW + vpdf_TN + vpdf_TS -
		//	     vpdf_B - vpdf_BE - vpdf_BW - vpdf_BN - vpdf_BS;
		vuz = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_T,VADD(vpdf_TE,vpdf_TW)),VADD(vpdf_TN,vpdf_TS)),vpdf_B),vpdf_BE),vpdf_BW),vpdf_BN),vpdf_BS);

		//vdens = vpdf_C +
		//	    vpdf_N  + vpdf_E  + vpdf_S  + vpdf_W  +
		//	    vpdf_NE + vpdf_SE + vpdf_SW + vpdf_NW +
		//	    vpdf_T  + vpdf_TN + vpdf_TE + vpdf_TS + vpdf_TW +
		//	    vpdf_B  + vpdf_BN + vpdf_BE + vpdf_BS + vpdf_BW;
		vdens = VADD(VADD(VADD(VADD(VADD(VADD(VADD(VADD(VADD(vpdf_C,VADD(vpdf_N,vpdf_E)),VADD(vpdf_S,vpdf_W)),VADD(vpdf_NE,vpdf_SE)),VADD(vpdf_SW,vpdf_NW)),VADD(vpdf_T,vpdf_TN)),VADD(vpdf_TE,vpdf_TS)),VADD(vpdf_TW,vpdf_B)),VADD(vpdf_BN,vpdf_BE)),VADD(vpdf_BS,vpdf_BW));

		//vdir_indep_trm = vdens - (vux * vux + vuy * vuy + vuz * vuz) * VTHREE_HALF;
		vdir_indep_trm = VSUB(vdens,VMUL(VADD(VADD(VMUL(vux,vux),VMUL(vuy,vuy)),VMUL(vuz,vuz)),VTHREE_HALF));

		//src[I(index, D3Q19_C)]  =[UA] vpdf_C - vomegaEven * (vpdf_C - vw_0 * vdir_indep_trm);
		VSTU(ppdf_C,VSUB(vpdf_C,VMUL(vomegaEven,VSUB(vpdf_C,VMUL(vw_0,vdir_indep_trm)))));

		//vw_1_indep = vw_1 * vdir_indep_trm;
		vw_1_indep = VMUL(vw_1,vdir_indep_trm);

		#define COLLIDE_AA_V(tmpVui, dir1, dir2) \
			vui = tmpVui; \
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(JOIN(vpdf_,dir1),JOIN(vpdf_,dir2))),VMUL(vui,VMUL(vui,vw_1_nine_half))),vw_1_indep));\
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(JOIN(vpdf_,dir1),JOIN(vpdf_,dir2))),VMUL(vui,vw_1_x3)));\
			VSTU(JOIN(ppdf_,dir2),VSUB(VSUB(JOIN(vpdf_,dir1),vevenPart),voddPart));\
			VSTU(JOIN(ppdf_,dir1),VADD(VSUB(JOIN(vpdf_,dir2),vevenPart),voddPart));

		COLLIDE_AA_V(vuy, N, S)
		COLLIDE_AA_V(vux, E, W)
		COLLIDE_AA_V(vuz, T, B)

		#undef COLLIDE_AA_V

		//vw_2_indep = vw_2 * vdir_indep_trm;
		vw_2_indep = VMUL(vw_2,vdir_indep_trm);

		// collide axis unaligned pdfs vectorized
		#define COLLIDE_UA_V(tmpVui, dir1, dir2) \
			vui = tmpVui; \
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(JOIN(vpdf_,dir1),JOIN(vpdf_,dir2))),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));\
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(JOIN(vpdf_,dir1),JOIN(vpdf_,dir2))),VMUL(vui,vw_2_x3)));\
			VSTU(JOIN(ppdf_,dir2),VSUB(VSUB(JOIN(vpdf_,dir1),vevenPart),voddPart)); \
			VSTU(JOIN(ppdf_,dir1),VADD(VSUB(JOIN(vpdf_,dir2),vevenPart),voddPart));

		COLLIDE_UA_V(VSUB(vuy,vux), NW, SE)
		COLLIDE_UA_V(VADD(vux,vuy), NE, SW)
		COLLIDE_UA_V(VSUB(vuz,vux), TW, BE)
		COLLIDE_UA_V(VADD(vux,vuz), TE, BW)
		COLLIDE_UA_V(VSUB(vuz,vuy), TS, BN)
		COLLIDE_UA_V(VADD(vuy,vuz), TN, BS)

		#undef COLLIDE_UA_V

		#define X(name, idx, idxinv, _x, _y, _z)	JOIN(ppdf_,name) +=(VSIZE * N_D3Q19);
			D3Q19_LIST
		#undef X
	} // loop over fluid nodes

	for (int index = indexStopVec; index < indexStop; ++index) {

		#define X(name, idx, idxinv, _x, _y, _z)	JOIN(pdf_,name) = *(JOIN(ppdf_,name));
			D3Q19_LIST
		#undef X

		ux = pdf_E + pdf_NE + pdf_SE + pdf_TE + pdf_BE -
			pdf_W - pdf_NW - pdf_SW - pdf_TW - pdf_BW;
		uy = pdf_N + pdf_NE + pdf_NW + pdf_TN + pdf_BN -
			pdf_S - pdf_SE - pdf_SW - pdf_TS - pdf_BS;
		uz = pdf_T + pdf_TE + pdf_TW + pdf_TN + pdf_TS -
			pdf_B - pdf_BE - pdf_BW - pdf_BN - pdf_BS;

		dens = pdf_C +
			pdf_N  + pdf_E  + pdf_S  + pdf_W  +
			pdf_NE + pdf_SE + pdf_SW + pdf_NW +
			pdf_T  + pdf_TN + pdf_TE + pdf_TS + pdf_TW +
			pdf_B  + pdf_BN + pdf_BE + pdf_BS + pdf_BW;

		dir_indep_trm = dens - (ux * ux + uy * uy + uz * uz)*3.0/2.0;

		// direction: w_0
		*ppdf_C  = pdf_C - omegaEven*(pdf_C - w_0*dir_indep_trm);

		// direction: w_1
		w_1_indep = w_1*dir_indep_trm;

		#define COLLIDE_AA_S(tmpUi, dir1, dir2) \
			ui = tmpUi; \
			evenPart = omegaEven * (0.5 * (JOIN(pdf_,dir1) + JOIN(pdf_,dir2)) - ui * ui * w_1_nine_half - w_1_indep); \
			oddPart  = omegaOdd  * (0.5 * (JOIN(pdf_,dir1) - JOIN(pdf_,dir2)) - ui * w_1_x3); \
			*(JOIN(ppdf_,dir2))  = JOIN(pdf_,dir1) - evenPart - oddPart; \
			*(JOIN(ppdf_,dir1))  = JOIN(pdf_,dir2) - evenPart + oddPart;

		COLLIDE_AA_S(uy, N, S)
		COLLIDE_AA_S(ux, E, W)
		COLLIDE_AA_S(uz, T, B)

		#undef COLLIDE_AA_S

		// direction: w_2
		w_2_indep = w_2*dir_indep_trm;

		#define COLLIDE_UA_S(tmpUi, dir1, dir2) \
			ui = tmpUi; \
			evenPart = omegaEven * (0.5 * (JOIN(pdf_,dir1) + JOIN(pdf_,dir2)) - ui * ui * w_2_nine_half - w_2_indep); \
			oddPart  = omegaOdd  * (0.5 * (JOIN(pdf_,dir1) - JOIN(pdf_,dir2)) - ui * w_2_x3); \
			*(JOIN(ppdf_,dir2)) = JOIN(pdf_,dir1) - evenPart - oddPart; \
			*(JOIN(ppdf_,dir1)) = JOIN(pdf_,dir2) - evenPart + oddPart;

		COLLIDE_UA_S((-ux + uy), NW, SE)
		COLLIDE_UA_S(( ux + uy), NE, SW)
		COLLIDE_UA_S((-ux + uz), TW, BE)
		COLLIDE_UA_S(( ux + uz), TE, BW)
		COLLIDE_UA_S((-uy + uz), TS, BN)
		COLLIDE_UA_S(( uy + uz), TN, BS)

		#undef COLLIDE_UA_S

		#define X(name, idx, idxinv, _x, _y, _z)	JOIN(ppdf_,name)++;
			D3Q19_LIST
		#undef X
	} // loop over fluid nodes

	#undef I

	return;
}

static void KernelOdd(LatticeDesc * ld, KernelData * kernelData, CaseData * cd, int * threadIndices)
{

	Assert(ld != NULL);
	Assert(kernelData != NULL);
	Assert(cd != NULL);

	Assert(cd->Omega > 0.0);
	Assert(cd->Omega < 2.0);

	KernelData * kd = (KernelData *)kernelData;
	KernelDataList * kdl = KDL(kernelData);
	KernelDataListRia * kdlr = KDLR(kernelData);
	PdfT omega = cd->Omega;
	PdfT omegaEven = omega;

	PdfT magicParam = 1.0 / 12.0;
	PdfT omegaOdd = 1.0 / (0.5 + magicParam / (1.0 / omega - 0.5));

	PdfT evenPart = 0.0;
	PdfT oddPart = 0.0;
	PdfT dir_indep_trm = 0.0;

	const PdfT w_0 = 1.0 /  3.0;
	const PdfT w_1 = 1.0 / 18.0;
	const PdfT w_2 = 1.0 / 36.0;

	const PdfT w_1_x3 = w_1 * 3.0;	const PdfT w_1_nine_half = w_1 * 9.0 / 2.0;	PdfT w_1_indep = 0.0;
	const PdfT w_2_x3 = w_2 * 3.0;	const PdfT w_2_nine_half = w_2 * 9.0 / 2.0;	PdfT w_2_indep = 0.0;

	PdfT ux, uy, uz, ui;
	PdfT dens;

	VPDFT VONE_HALF = VSET(0.5);
	VPDFT VTHREE_HALF = VSET(3.0 / 2.0);

	VPDFT vw_1_indep, vw_2_indep;
	VPDFT vw_0 = VSET(w_0);
	VPDFT vw_1 = VSET(w_1);
	VPDFT vw_2 = VSET(w_2);

	VPDFT vw_1_x3 = VSET(w_1_x3);
	VPDFT vw_2_x3 = VSET(w_2_x3);
	VPDFT vw_1_nine_half = VSET(w_1_nine_half);
	VPDFT vw_2_nine_half = VSET(w_2_nine_half);

	VPDFT vux, vuy, vuz, vui;
	VPDFT vdens;

	VPDFT vevenPart, voddPart, vdir_indep_trm;

	VPDFT vomegaEven = VSET(omegaEven);
	VPDFT vomegaOdd  = VSET(omegaOdd);

	// Declare pdf_N, pdf_E, pdf_S, pdf_W, ...
	#define X(name, idx, idxinv, x, y, z) \
		PdfT JOIN(pdf_,name) = 0; \
		PdfT * JOIN(ppdf_,name) = NULL; \
		VPDFT JOIN(vpdf_,name);
		D3Q19_LIST
	#undef X

	#define X(name, idx, idxinv, x, y, z) \
		__m256i JOIN(vgatheridx_,name) = _mm256_set1_epi32(0);
		D3Q19_LIST_WO_C
	#undef X

	__m256i vgatherinc = VSETI32(VSIZE * N_D3Q19);

	uint32_t * consecNodes = kdlr->ConsecNodes;
	uint32_t consecIndex = 0;
	uint32_t consecValue = 0;

	PdfT * src = kd->Pdfs[0];

	int nCells = kdl->nCells;

	int adjListIndex;
	uint32_t * adjList = kdl->AdjList;

	int threadId = 0;

#ifdef _OPENMP
	threadId = omp_get_thread_num();
#endif
	consecIndex = kdlr->ConsecThreadIndices[threadId];
	consecValue = 0;

	int nFluidThread = threadIndices[threadId + 1] - threadIndices[threadId];

	int indexStart = threadIndices[threadId];
	int indexStop  = threadIndices[threadId] + nFluidThread;

	#define I(index, dir)	P_INDEX_3((nCells), (index), (dir))
	#define ADJ_LIST(dir) adjList[adjListIndex + (dir * VSIZE)]

	int offset_ppdf_C = -1; //dummy init to detect errors

	for (int index = indexStart; index < indexStop; index += 1) {

		if (consecValue > 0) {
			--consecValue;
			// Increment all pdf pointers by an offset. If the previous iteration was
			// scalar, increment only by one. If the previous iteration was vectorized,
			// increment by the vector width. These offsets are set in the corresponding
			// if branches.

			//increment offsets

			#define X(name, idx, idxinv, _x, _y, _z) JOIN(vgatheridx_,name) = VADDI32(JOIN(vgatheridx_,name), vgatherinc);
				D3Q19_LIST_WO_C
			#undef X

			//printVector(vgatheridx_N);

			ppdf_C += offset_ppdf_C;

		}
		else {
			// Load new pointers to PDFs of local cell:
			Assert(consecIndex < nConsecNodes);

			consecValue = consecNodes[consecIndex] - 1;

			adjListIndex = (index - (index % VSIZE)) * N_D3Q19_IDX + (index % VSIZE);
			#define X(name, idx, idxinv, _x, _y, _z)	JOIN(vgatheridx_,name) = VLIU(&(ADJ_LIST(idxinv)));
				D3Q19_LIST_WO_C
			#undef X

			ppdf_C = &(src[P_INDEX_3(nCells, index, D3Q19_C)]);
			++consecIndex;
		}

		if (consecValue >= (VSIZE - 1)) {
			// Vectorized part.

			#if (SOFTWARE_PREFETCH_LOOKAHEAD_L2 > 0)
				int const indexPrefetchL2 = index + VSIZE * SOFTWARE_PREFETCH_LOOKAHEAD_L2;
				// make sure that adjList access is never out of bounds since it is an actual memory access and no prefetch
				if (indexPrefetchL2 < indexStop){
					// update pointers from adjacency list only if necessary
					if (consecValue >= (SOFTWARE_PREFETCH_LOOKAHEAD_L2 * VSIZE + VSIZE - 1)) {
						#define INCR_PTR(name)		(VADDI32(JOIN(vgatheridx_,name), VMULI32(vgatherinc, VSETI32(SOFTWARE_PREFETCH_LOOKAHEAD_L2))))
						#define X(name, idx, idxinv, _x, _y, _z) VPG32(INCR_PTR(name), (char const *) src, 8, _MM_HINT_T1);
							D3Q19_LIST_WO_C
						#undef X
						#undef INCR_PTR
					}
					else {
						adjListIndex = (indexPrefetchL2 - (indexPrefetchL2 % VSIZE)) * N_D3Q19_IDX + (indexPrefetchL2 % VSIZE);
						#define X(name, idx, idxinv, _x, _y, _z) VPG32(VLIU(&ADJ_LIST(idxinv)), (char const *) src, 8, _MM_HINT_T1);
							D3Q19_LIST_WO_C
						#undef X
					}

					_mm_prefetch((char const *) &(src[P_INDEX_3(nCells, indexPrefetchL2, D3Q19_C)]), _MM_HINT_T1);
				}
			#endif

			#if (SOFTWARE_PREFETCH_LOOKAHEAD_L1 > 0)
				int const indexPrefetchL1 = index + VSIZE * SOFTWARE_PREFETCH_LOOKAHEAD_L1;
				// make sure that adjList access is never out of bounds since it is an actual memory access and no prefetch
				if (indexPrefetchL1 < indexStop){
					// update pointers from adjacency list only if necessary
					if (consecValue > (SOFTWARE_PREFETCH_LOOKAHEAD_L1 * VSIZE + VSIZE - 1)) {
						#define INCR_PTR(name)		(VADDI32(JOIN(vgatheridx_,name), VMULI32(vgatherinc, VSETI32(SOFTWARE_PREFETCH_LOOKAHEAD_L1))))
						#define X(name, idx, idxinv, _x, _y, _z) VPG32(INCR_PTR(name), (char const *) src, 8, _MM_HINT_T0);
							D3Q19_LIST_WO_C
						#undef X
						#undef INCR_PTR
					}
					else {
						adjListIndex = (indexPrefetchL1 - (indexPrefetchL1 % VSIZE)) * N_D3Q19_IDX + (indexPrefetchL1 % VSIZE);
						#define X(name, idx, idxinv, _x, _y, _z) VPG32(VLIU(&ADJ_LIST(idxinv)), (char const *) src, 8, _MM_HINT_T0);
							D3Q19_LIST_WO_C
						#undef X
					}

					_mm_prefetch((char const *) &(src[P_INDEX_3(nCells, indexPrefetchL1, D3Q19_C)]), _MM_HINT_T0);
				}
			#endif

			#define X(name, idx, idxinv, _x, _y, _z)	JOIN(vpdf_,name) = VG32(JOIN(vgatheridx_,name), src, 8);
				D3Q19_LIST_WO_C
			#undef X

			vpdf_C = VLDU(ppdf_C);

			//vux = vpdf_E + vpdf_NE + vpdf_SE + vpdf_TE + vpdf_BE -
			//      vpdf_W - vpdf_NW - vpdf_SW - vpdf_TW - vpdf_BW;
			vux = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_E,VADD(vpdf_NE,vpdf_SE)),VADD(vpdf_TE,vpdf_BE)),vpdf_W),vpdf_NW),vpdf_SW),vpdf_TW),vpdf_BW);
			//vuy = vpdf_N + vpdf_NE + vpdf_NW + vpdf_TN + vpdf_BN -
			//      vpdf_S - vpdf_SE - vpdf_SW - vpdf_TS - vpdf_BS;
			vuy = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_N,VADD(vpdf_NE,vpdf_NW)),VADD(vpdf_TN,vpdf_BN)),vpdf_S),vpdf_SE),vpdf_SW),vpdf_TS),vpdf_BS);
			//vuz = vpdf_T + vpdf_TE + vpdf_TW + vpdf_TN + vpdf_TS -
			//      vpdf_B - vpdf_BE - vpdf_BW - vpdf_BN - vpdf_BS;
			vuz = VSUB(VSUB(VSUB(VSUB(VSUB(VADD(VADD(vpdf_T,VADD(vpdf_TE,vpdf_TW)),VADD(vpdf_TN,vpdf_TS)),vpdf_B),vpdf_BE),vpdf_BW),vpdf_BN),vpdf_BS);

			//vdens = vpdf_C +
			//        vpdf_N  + vpdf_E  + vpdf_S  + vpdf_W  +
			//        vpdf_NE + vpdf_SE + vpdf_SW + vpdf_NW +
			//        vpdf_T  + vpdf_TN + vpdf_TE + vpdf_TS + vpdf_TW +
			//        vpdf_B  + vpdf_BN + vpdf_BE + vpdf_BS + vpdf_BW;
			vdens = VADD(VADD(VADD(VADD(VADD(VADD(VADD(VADD(VADD(vpdf_C,VADD(vpdf_N,vpdf_E)),VADD(vpdf_S,vpdf_W)),VADD(vpdf_NE,vpdf_SE)),
										VADD(vpdf_SW,vpdf_NW)),VADD(vpdf_T,vpdf_TN)),VADD(vpdf_TE,vpdf_TS)),VADD(vpdf_TW,vpdf_B)),VADD(vpdf_BN,vpdf_BE)),VADD(vpdf_BS,vpdf_BW));

			//vdir_indep_trm = vdens - (vux * vux + vuy * vuy + vuz * vuz) * VTHREE_HALF;
			vdir_indep_trm = VSUB(vdens,VMUL(VADD(VADD(VMUL(vux,vux),VMUL(vuy,vuy)),VMUL(vuz,vuz)),VTHREE_HALF));

			//src[I(index, D3Q19_C)]  =[UA] vpdf_C - vomegaEven * (vpdf_C - vw_0 * vdir_indep_trm);
			VSTU(ppdf_C,VSUB(vpdf_C,VMUL(vomegaEven,VSUB(vpdf_C,VMUL(vw_0,vdir_indep_trm)))));

			// collide axis aligend pdfs vectorized
			#define SCAT(offsets, vsrc) VS32(src, offsets, vsrc, 8)

			//vw_1_indep = vw_1 * vdir_indep_trm;
			vw_1_indep = VMUL(vw_1,vdir_indep_trm);

			// collide axis aligend pdfs vectorized
			#define COLLIDE_AA_V(tmpVui, dir1, dir2) \
				vui = tmpVui; \
				vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(JOIN(vpdf_,dir1),JOIN(vpdf_,dir2))),VMUL(vui,VMUL(vui,vw_1_nine_half))),vw_1_indep));\
				voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(JOIN(vpdf_,dir1),JOIN(vpdf_,dir2))),VMUL(vui,vw_1_x3)));\
				SCAT(JOIN(vgatheridx_,dir2),VSUB(VSUB(JOIN(vpdf_,dir1),vevenPart),voddPart));\
				SCAT(JOIN(vgatheridx_,dir1),VADD(VSUB(JOIN(vpdf_,dir2),vevenPart),voddPart));

			COLLIDE_AA_V(vuy, N, S)
			COLLIDE_AA_V(vux, E, W)
			COLLIDE_AA_V(vuz, T, B)

			#undef COLLIDE_AA_V

			//vw_2_indep = vw_2 * vdir_indep_trm;
			vw_2_indep = VMUL(vw_2,vdir_indep_trm);

			// collide axis unaligned pdfs vectorized
			#define COLLIDE_UA_V(tmpVui, dir1, dir2) \
				vui = tmpVui; \
				vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(JOIN(vpdf_,dir1),JOIN(vpdf_,dir2))),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));\
				voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(JOIN(vpdf_,dir1),JOIN(vpdf_,dir2))),VMUL(vui,vw_2_x3)));\
				SCAT(JOIN(vgatheridx_,dir2),VSUB(VSUB(JOIN(vpdf_,dir1),vevenPart),voddPart)); \
				SCAT(JOIN(vgatheridx_,dir1),VADD(VSUB(JOIN(vpdf_,dir2),vevenPart),voddPart));

			COLLIDE_UA_V(VSUB(vuy,vux), NW, SE)
			COLLIDE_UA_V(VADD(vux,vuy), NE, SW)
			COLLIDE_UA_V(VSUB(vuz,vux), TW, BE)
			COLLIDE_UA_V(VADD(vux,vuz), TE, BW)
			COLLIDE_UA_V(VSUB(vuz,vuy), TS, BN)
			COLLIDE_UA_V(VADD(vuy,vuz), TN, BS)

			#undef COLLIDE_UA_V
			#undef SCAT

			consecValue   -= (VSIZE - 1);
			index         += (VSIZE - 1);
			offset_ppdf_C  = VSIZE * N_D3Q19;

		}
		else {
			// Scalar part.

			adjListIndex = (index - (index % VSIZE)) * N_D3Q19_IDX + (index % VSIZE);
			#define X(name, idx, idxinv, _x, _y, _z)	JOIN(ppdf_,name) = &(src[ADJ_LIST(idxinv)]);
				D3Q19_LIST_WO_C
			#undef X
			#define X(name, idx, idxinv, _x, _y, _z)	JOIN(pdf_,name) = *(JOIN(ppdf_,name));
				D3Q19_LIST_WO_C
			#undef X

			pdf_C = *ppdf_C;

			ux = pdf_E + pdf_NE + pdf_SE + pdf_TE + pdf_BE -
				pdf_W - pdf_NW - pdf_SW - pdf_TW - pdf_BW;
			uy = pdf_N + pdf_NE + pdf_NW + pdf_TN + pdf_BN -
				pdf_S - pdf_SE - pdf_SW - pdf_TS - pdf_BS;
			uz = pdf_T + pdf_TE + pdf_TW + pdf_TN + pdf_TS -
				pdf_B - pdf_BE - pdf_BW - pdf_BN - pdf_BS;

			dens = pdf_C +
				pdf_N  + pdf_E  + pdf_S  + pdf_W  +
				pdf_NE + pdf_SE + pdf_SW + pdf_NW +
				pdf_T  + pdf_TN + pdf_TE + pdf_TS + pdf_TW +
				pdf_B  + pdf_BN + pdf_BE + pdf_BS + pdf_BW;

			dir_indep_trm = dens - (ux * ux + uy * uy + uz * uz)*3.0/2.0;

			// direction: w_0
			*ppdf_C = pdf_C - omegaEven * (pdf_C - w_0 * dir_indep_trm);

			// direction: w_1
			w_1_indep = w_1 * dir_indep_trm;

			#define COLLIDE_AA_S(tmpUi, dir1, dir2) \
				ui = tmpUi; \
				evenPart = omegaEven * (0.5 * (JOIN(pdf_,dir1) + JOIN(pdf_,dir2)) - ui * ui * w_1_nine_half - w_1_indep); \
				oddPart  = omegaOdd  * (0.5 * (JOIN(pdf_,dir1) - JOIN(pdf_,dir2)) - ui * w_1_x3); \
				*(JOIN(ppdf_,dir2))  = JOIN(pdf_,dir1) - evenPart - oddPart; \
				*(JOIN(ppdf_,dir1))  = JOIN(pdf_,dir2) - evenPart + oddPart;

			COLLIDE_AA_S(uy, N, S)
			COLLIDE_AA_S(ux, E, W)
			COLLIDE_AA_S(uz, T, B)

			#undef COLLIDE_AA_S

			// direction: w_2
			w_2_indep = w_2 * dir_indep_trm;

			#define COLLIDE_UA_S(tmpUi, dir1, dir2) \
				ui = tmpUi; \
				evenPart = omegaEven * (0.5 * (JOIN(pdf_,dir1) + JOIN(pdf_,dir2)) - ui * ui * w_2_nine_half - w_2_indep); \
				oddPart  = omegaOdd  * (0.5 * (JOIN(pdf_,dir1) - JOIN(pdf_,dir2)) - ui * w_2_x3); \
				*(JOIN(ppdf_,dir2)) = JOIN(pdf_,dir1) - evenPart - oddPart; \
				*(JOIN(ppdf_,dir1)) = JOIN(pdf_,dir2) - evenPart + oddPart;

			COLLIDE_UA_S((-ux + uy), NW, SE)
			COLLIDE_UA_S(( ux + uy), NE, SW)
			COLLIDE_UA_S((-ux + uz), TW, BE)
			COLLIDE_UA_S(( ux + uz), TE, BW)
			COLLIDE_UA_S((-uy + uz), TS, BN)
			COLLIDE_UA_S(( uy + uz), TN, BS)

			#undef COLLIDE_UA_S

			offset_ppdf_C = 1;
		}

	} // loop over fluid nodes

#undef ADJ_LIST
#undef I
}
