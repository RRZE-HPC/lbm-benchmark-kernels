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
#include "BenchKernelD3Q19ListAaPvCommon.h"

#include "Memory.h"
#include "Vtk.h"
#include "Vector.h"
#include "LikwidIf.h"

#include <inttypes.h>
#include <math.h>

#ifdef _OPENMP
	#include <omp.h>
#endif


static void KernelEven(LatticeDesc * ld, KernelData * kernelData, CaseData * cd);
static void KernelOdd( LatticeDesc * ld, KernelData * kernelData, CaseData * cd);

void FNAME(D3Q19ListAaPvKernel)(LatticeDesc * ld, KernelData * kernelData, CaseData * cd)
{

	Assert(ld != NULL);
	Assert(kernelData != NULL);
	Assert(cd != NULL);

	Assert(cd->Omega > 0.0);
	Assert(cd->Omega < 2.0);

#if defined(VTK_OUTPUT) || defined(STATISTICS) || defined(VERIFICATION)
	KernelData     * kd  = (KernelData *)kernelData;
	KernelDataList * kdl = KDL(kernelData);
#endif

	int maxIterations = cd->MaxIterations;

	int nThreads = 1;
#ifdef _OPENMP
	nThreads = omp_get_max_threads();
#endif

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

	// TODO: outer openmp parallel

	for(int iter = 0; iter < maxIterations; iter += 2) {

		// ---------------------------------------------------
		// even time step
		// ---------------------------------------------------

		X_LIKWID_START("list-aa-pv-even");

		#ifdef _OPENMP
			#pragma omp parallel default(none) shared(ld, kernelData, cd)
		#endif
		{
			KernelEven(ld, kernelData, cd);
		}

		X_LIKWID_STOP("list-aa-pv-even");

		#ifdef VERIFICATION
			kdl->Iteration = iter;
			kd->PdfsActive = kd->Pdfs[0];
			KernelAddBodyForce(kd, ld, cd);
		#endif

		// ---------------------------------------------------
		// odd time step
		// ---------------------------------------------------

		X_LIKWID_START("list-aa-pv-odd");

		#ifdef _OPENMP
			#pragma omp parallel default(none) shared(ld, kernelData, cd)
		#endif
		{
			KernelOdd(ld, kernelData, cd);
		}

		X_LIKWID_STOP("list-aa-pv-odd");


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

	return;
}

static void KernelEven(LatticeDesc * ld, KernelData * kernelData, CaseData * cd)
{
	Assert(ld != NULL);
	Assert(kernelData != NULL);
	Assert(cd != NULL);

	Assert(cd->Omega > 0.0);
	Assert(cd->Omega < 2.0);

	KernelData        * kd   = (KernelData *)kernelData;
	KernelDataList    * kdl  = KDL(kernelData);
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

	PdfT ui;

	PdfT ux, uy, uz;
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
	#define X(name, idx, idxinv, x, y, z)	PdfT JOIN(pdf_,name); VPDFT JOIN(vpdf_,name);
	D3Q19_LIST
	#undef X

	PdfT * src = kd->Pdfs[0];

	int nCells = kdl->nCells;

	int threadId = 0;
#ifdef _OPENMP
	threadId =  omp_get_thread_num();
#endif

	int * threadIndices = kdlr->FluidNodeThreadIndices;

	int nFluidThread = threadIndices[threadId + 1] - threadIndices[threadId];
	int nFluidVec    = nFluidThread - (nFluidThread % VSIZE);

	int indexStartVec = threadIndices[threadId];
	int indexStopVec  = threadIndices[threadId] + nFluidVec;
	int indexStop     = threadIndices[threadId] + nFluidThread;

	#define I(index, dir)	P_INDEX_3((nCells), (index), (dir))

	for (int index = indexStartVec; index < indexStopVec; index += VSIZE) {


		#define X(name, idx, idxinv, _x, _y, _z)	JOIN(vpdf_,name) = VLDU(&src[I(index, idx)]);
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

		//src[I(index, D3Q19_C)             ]  =[UA] vpdf_C - vomegaEven * (vpdf_C - vw_0 * vdir_indep_trm);
		VSTU(&src[I(index, D3Q19_C)],VSUB(vpdf_C,VMUL(vomegaEven,VSUB(vpdf_C,VMUL(vw_0,vdir_indep_trm)))));

		//vw_1_indep = vw_1 * vdir_indep_trm;
		vw_1_indep = VMUL(vw_1,vdir_indep_trm);

		//vui = vuy;
		vui = vuy;
		//vevenPart = vomegaEven * (VONE_HALF * (vpdf_N + vpdf_S) - vui * vui * vw_1_nine_half - vw_1_indep);
		vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_N,vpdf_S)),VMUL(vui,VMUL(vui,vw_1_nine_half))),vw_1_indep));
		//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_N - vpdf_S) - vui * vw_1_x3);
		voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_N,vpdf_S)),VMUL(vui,vw_1_x3)));
		//src[I(index, D3Q19_S)]  =[UA] vpdf_N - vevenPart - voddPart;
		VSTU(&src[I(index, D3Q19_S)],VSUB(VSUB(vpdf_N,vevenPart),voddPart));
		//src[I(index, D3Q19_N)]  =[UA] vpdf_S - vevenPart + voddPart;
		VSTU(&src[I(index, D3Q19_N)],VADD(VSUB(vpdf_S,vevenPart),voddPart));

		//vui = vux;
		vui = vux;
		//vevenPart = vomegaEven * (VONE_HALF * (vpdf_E + vpdf_W) - vui * vui * vw_1_nine_half - vw_1_indep);
		vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_E,vpdf_W)),VMUL(vui,VMUL(vui,vw_1_nine_half))),vw_1_indep));
		//voddPart =  vomegaOdd  * (VONE_HALF * (vpdf_E - vpdf_W) - vui * vw_1_x3 );
		voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_E,vpdf_W)),VMUL(vui,vw_1_x3)));
		//src[I(index, D3Q19_W)]  =[UA] vpdf_E - vevenPart - voddPart;
		VSTU(&src[I(index, D3Q19_W)],VSUB(VSUB(vpdf_E,vevenPart),voddPart));
		//src[I(index, D3Q19_E)]  =[UA] vpdf_W - vevenPart + voddPart;
		VSTU(&src[I(index, D3Q19_E)],VADD(VSUB(vpdf_W,vevenPart),voddPart));

		//vui = vuz;
		vui = vuz;
		//vevenPart = vomegaEven * (VONE_HALF * (vpdf_T + vpdf_B) - vui * vui * vw_1_nine_half - vw_1_indep);
		vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_T,vpdf_B)),VMUL(vui,VMUL(vui,vw_1_nine_half))),vw_1_indep));
		//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_T - vpdf_B) - vui * vw_1_x3);
		voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_T,vpdf_B)),VMUL(vui,vw_1_x3)));
		//src[I(index, D3Q19_B)]  =[UA] vpdf_T - vevenPart - voddPart;
		VSTU(&src[I(index, D3Q19_B)],VSUB(VSUB(vpdf_T,vevenPart),voddPart));
		//src[I(index, D3Q19_T)]  =[UA] vpdf_B - vevenPart + voddPart;
		VSTU(&src[I(index, D3Q19_T)],VADD(VSUB(vpdf_B,vevenPart),voddPart));

		//vw_2_indep = vw_2 * vdir_indep_trm;
		vw_2_indep = VMUL(vw_2,vdir_indep_trm);

		//vui = vuy - vux;
		vui = VSUB(vuy,vux);
		//vevenPart = vomegaEven * (VONE_HALF * (vpdf_NW + vpdf_SE) - vui * vui * vw_2_nine_half - vw_2_indep);
		vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_NW,vpdf_SE)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
		//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_NW - vpdf_SE) - vui * vw_2_x3);
		voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_NW,vpdf_SE)),VMUL(vui,vw_2_x3)));
		//src[I(index, D3Q19_SE)] =[UA] vpdf_NW - vevenPart - voddPart;
		VSTU(&src[I(index, D3Q19_SE)],VSUB(VSUB(vpdf_NW,vevenPart),voddPart));
		//src[I(index, D3Q19_NW)] =[UA] vpdf_SE - vevenPart + voddPart;
		VSTU(&src[I(index, D3Q19_NW)],VADD(VSUB(vpdf_SE,vevenPart),voddPart));

		//vui = vux + vuy;
		vui = VADD(vux,vuy);
		//vevenPart = vomegaEven * (VONE_HALF * (vpdf_NE + vpdf_SW) - vui * vui * vw_2_nine_half - vw_2_indep);
		vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_NE,vpdf_SW)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
		//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_NE - vpdf_SW) - vui * vw_2_x3);
		voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_NE,vpdf_SW)),VMUL(vui,vw_2_x3)));
		//src[I(index, D3Q19_SW)] =[UA] vpdf_NE - vevenPart - voddPart;
		VSTU(&src[I(index, D3Q19_SW)],VSUB(VSUB(vpdf_NE,vevenPart),voddPart));
		//src[I(index, D3Q19_NE)] =[UA] vpdf_SW - vevenPart + voddPart;
		VSTU(&src[I(index, D3Q19_NE)],VADD(VSUB(vpdf_SW,vevenPart),voddPart));

		//vui = vuz - vux;
		vui = VSUB(vuz,vux);
		//vevenPart = vomegaEven * (VONE_HALF * (vpdf_TW + vpdf_BE) - vui * vui * vw_2_nine_half - vw_2_indep);
		vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TW,vpdf_BE)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
		//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_TW - vpdf_BE) - vui * vw_2_x3);
		voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TW,vpdf_BE)),VMUL(vui,vw_2_x3)));
		//src[I(index, D3Q19_BE)] =[UA] vpdf_TW - vevenPart - voddPart;
		VSTU(&src[I(index, D3Q19_BE)],VSUB(VSUB(vpdf_TW,vevenPart),voddPart));
		//src[I(index, D3Q19_TW)] =[UA] vpdf_BE - vevenPart + voddPart;
		VSTU(&src[I(index, D3Q19_TW)],VADD(VSUB(vpdf_BE,vevenPart),voddPart));

		//vui = vux + vuz;
		vui = VADD(vux,vuz);
		//vevenPart = vomegaEven * (VONE_HALF * (vpdf_TE + vpdf_BW) - vui * vui * vw_2_nine_half - vw_2_indep);
		vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TE,vpdf_BW)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
		//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_TE - vpdf_BW) - vui * vw_2_x3);
		voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TE,vpdf_BW)),VMUL(vui,vw_2_x3)));
		//src[I(index, D3Q19_BW)] =[UA] vpdf_TE - vevenPart - voddPart;
		VSTU(&src[I(index, D3Q19_BW)],VSUB(VSUB(vpdf_TE,vevenPart),voddPart));
		//src[I(index, D3Q19_TE)] =[UA] vpdf_BW - vevenPart + voddPart;
		VSTU(&src[I(index, D3Q19_TE)],VADD(VSUB(vpdf_BW,vevenPart),voddPart));

		//vui = vuz - vuy;
		vui = VSUB(vuz,vuy);
		//vevenPart = vomegaEven * (VONE_HALF * (vpdf_TS + vpdf_BN) - vui * vui * vw_2_nine_half - vw_2_indep);
		vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TS,vpdf_BN)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
		//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_TS - vpdf_BN) - vui * vw_2_x3);
		voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TS,vpdf_BN)),VMUL(vui,vw_2_x3)));
		//src[I(index, D3Q19_BN)] =[UA] vpdf_TS - vevenPart - voddPart;
		VSTU(&src[I(index, D3Q19_BN)],VSUB(VSUB(vpdf_TS,vevenPart),voddPart));
		//src[I(index, D3Q19_TS)] =[UA] vpdf_BN - vevenPart + voddPart;
		VSTU(&src[I(index, D3Q19_TS)],VADD(VSUB(vpdf_BN,vevenPart),voddPart));

		//vui = vuy + vuz;
		vui = VADD(vuy,vuz);
		//vevenPart = vomegaEven * (VONE_HALF * (vpdf_TN + vpdf_BS) - vui * vui * vw_2_nine_half - vw_2_indep);
		vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TN,vpdf_BS)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
		//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_TN - vpdf_BS) - vui * vw_2_x3);
		voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TN,vpdf_BS)),VMUL(vui,vw_2_x3)));
		//src[I(index, D3Q19_BS)] =[UA] vpdf_TN - vevenPart - voddPart;
		VSTU(&src[I(index, D3Q19_BS)],VSUB(VSUB(vpdf_TN,vevenPart),voddPart));
		//src[I(index, D3Q19_TN)] =[UA] vpdf_BS - vevenPart + voddPart;
		VSTU(&src[I(index, D3Q19_TN)],VADD(VSUB(vpdf_BS,vevenPart),voddPart));

	} // loop over fluid nodes

	for (int index = indexStopVec; index < indexStop; ++index) {

		#define X(name, idx, idxinv, _x, _y, _z)	JOIN(pdf_,name) = src[I(index, idx)];
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
		src[I(index, D3Q19_C)             ]  = pdf_C - omegaEven*(pdf_C - w_0*dir_indep_trm);

		// direction: w_1
		w_1_indep = w_1*dir_indep_trm;

		ui = uy;
		evenPart = omegaEven*( 0.5*(pdf_N + pdf_S) - ui*ui*w_1_nine_half - w_1_indep );
		oddPart = omegaOdd*(0.5*(pdf_N - pdf_S) - ui*w_1_x3 );
		src[I(index, D3Q19_S)]  = pdf_N - evenPart - oddPart;
		src[I(index, D3Q19_N)]  = pdf_S - evenPart + oddPart;

		ui = ux;
		evenPart = omegaEven*( 0.5*(pdf_E + pdf_W) - ui*ui*w_1_nine_half - w_1_indep );
		oddPart = omegaOdd*(0.5*(pdf_E - pdf_W) - ui*w_1_x3 );
		src[I(index, D3Q19_W)]  = pdf_E - evenPart - oddPart;
		src[I(index, D3Q19_E)]  = pdf_W - evenPart + oddPart;

		ui = uz;
		evenPart = omegaEven*( 0.5*(pdf_T + pdf_B) - ui*ui*w_1_nine_half - w_1_indep );
		oddPart = omegaOdd*(0.5*(pdf_T - pdf_B) - ui*w_1_x3 );
		src[I(index, D3Q19_B)]  = pdf_T - evenPart - oddPart;
		src[I(index, D3Q19_T)]  = pdf_B - evenPart + oddPart;

		// direction: w_2
		w_2_indep = w_2*dir_indep_trm;

		ui = -ux + uy;
		evenPart = omegaEven*( 0.5*(pdf_NW + pdf_SE) - ui*ui*w_2_nine_half - w_2_indep );
		oddPart = omegaOdd*(0.5*(pdf_NW - pdf_SE) - ui*w_2_x3 );
		src[I(index, D3Q19_SE)] = pdf_NW - evenPart - oddPart;
		src[I(index, D3Q19_NW)] = pdf_SE - evenPart + oddPart;

		ui = ux + uy;
		evenPart = omegaEven*( 0.5*(pdf_NE + pdf_SW) - ui*ui*w_2_nine_half - w_2_indep );
		oddPart = omegaOdd*(0.5*(pdf_NE - pdf_SW) - ui*w_2_x3 );
		src[I(index, D3Q19_SW)] = pdf_NE - evenPart - oddPart;
		src[I(index, D3Q19_NE)] = pdf_SW - evenPart + oddPart;

		ui = -ux + uz;
		evenPart = omegaEven*( 0.5*(pdf_TW + pdf_BE) - ui*ui*w_2_nine_half - w_2_indep );
		oddPart = omegaOdd*(0.5*(pdf_TW - pdf_BE) - ui*w_2_x3 );
		src[I(index, D3Q19_BE)] = pdf_TW - evenPart - oddPart;
		src[I(index, D3Q19_TW)] = pdf_BE - evenPart + oddPart;

		ui = ux + uz;
		evenPart = omegaEven*( 0.5*(pdf_TE + pdf_BW) - ui*ui*w_2_nine_half - w_2_indep );
		oddPart = omegaOdd*(0.5*(pdf_TE - pdf_BW) - ui*w_2_x3 );
		src[I(index, D3Q19_BW)] = pdf_TE - evenPart - oddPart;
		src[I(index, D3Q19_TE)] = pdf_BW - evenPart + oddPart;

		ui = -uy + uz;
		evenPart = omegaEven*( 0.5*(pdf_TS + pdf_BN) - ui*ui*w_2_nine_half - w_2_indep );
		oddPart = omegaOdd*(0.5*(pdf_TS - pdf_BN) - ui*w_2_x3 );
		src[I(index, D3Q19_BN)] = pdf_TS - evenPart - oddPart;
		src[I(index, D3Q19_TS)] = pdf_BN - evenPart + oddPart;

		ui = uy + uz;
		evenPart = omegaEven*( 0.5*(pdf_TN + pdf_BS) - ui*ui*w_2_nine_half - w_2_indep );
		oddPart = omegaOdd*(0.5*(pdf_TN - pdf_BS) - ui*w_2_x3 );
		src[I(index, D3Q19_BS)] = pdf_TN - evenPart - oddPart;
		src[I(index, D3Q19_TN)] = pdf_BS - evenPart + oddPart;

	} // loop over fluid nodes

	#undef I

	return;
}

static void KernelOdd(LatticeDesc * ld, KernelData * kernelData, CaseData * cd)
{

	Assert(ld != NULL);
	Assert(kernelData != NULL);
	Assert(cd != NULL);

	Assert(cd->Omega > 0.0);
	Assert(cd->Omega < 2.0);

	KernelData        * kd   = (KernelData *)kernelData;
	KernelDataList    * kdl  = KDL(kernelData);
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

	PdfT ui;

	PdfT ux, uy, uz;
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
	#define X(name, idx, idxinv, x, y, z)	PdfT JOIN(pdf_,name); VPDFT JOIN(vpdf_,name);
	D3Q19_LIST
	#undef X

	// Declare pointers to pdfs ppdf_N, ppdf_E, ppdf_S, ppdf_W, ...
	#define X(name, idx, idxinv, x, y, z)	PdfT * JOIN(ppdf_,name) = NULL;
	D3Q19_LIST
	#undef X

	uint32_t nConsecNodes = kdlr->nConsecNodes;
	uint32_t * consecNodes = kdlr->ConsecNodes;
	uint32_t consecIndex = 0;
	uint32_t consecValue = 0;

#ifndef DEBUG
	UNUSED(nConsecNodes);
#endif

	PdfT * src = kd->Pdfs[0];

	int nCells = kdl->nCells;

	uint32_t adjListIndex;
	uint32_t * adjList = kdl->AdjList;

	int threadId = 0;

	#ifdef _OPENMP
		threadId = omp_get_thread_num();
	#endif

	consecIndex = kdlr->ConsecThreadIndices[threadId];
	consecValue = 0;

	int * threadIndices = kdlr->FluidNodeThreadIndices;

	int nFluidThread = threadIndices[threadId + 1] - threadIndices[threadId];

	int indexStart = threadIndices[threadId];
	int indexStop  = threadIndices[threadId] + nFluidThread;

	#define I(index, dir)	P_INDEX_3((nCells), (index), (dir))

	#define ADJ_LIST(dir) adjList[adjListIndex + (dir)]

	int pointerOffset = 1;

	for (int index = indexStart; index < indexStop; index += 1) {

		if (consecValue > 0) {
			--consecValue;
			// Increment all pdf pointers by an offset. If the previous iteration was
			// scalar, increment only by one. If the previous iteration was vectorized,
			// increment by the vector width. These offsets are set in the corresponding
			// if branches.
			#define X(name, idx, idxinv, _x, _y, _z)	JOIN(ppdf_,name) += pointerOffset;
			D3Q19_LIST
			#undef X
		}
		else {
			Assert(consecIndex < nConsecNodes);

			consecValue = consecNodes[consecIndex] - 1;
			// Load new pointers to PDFs of local cell:

			adjListIndex = index * N_D3Q19_IDX;

			#define X(name, idx, idxinv, _x, _y, _z)	JOIN(ppdf_,name) = &(src[adjList[adjListIndex + idxinv]]);
			D3Q19_LIST_WO_C
			#undef X

			ppdf_C = &(src[P_INDEX_3(nCells, index, D3Q19_C)]);
			++consecIndex;
		}

		#define X(name, idx, idxinv, _x, _y, _z)	JOIN(pdf_,name) = *JOIN(ppdf_,name);
		D3Q19_LIST
		#undef X

		if (consecValue >= (VSIZE - 1)) {
			// Vectorized part.

			#define X(name, idx, idxinv, _x, _y, _z)	JOIN(vpdf_,name) = VLDU(JOIN(ppdf_,name));
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

			adjListIndex = index * N_D3Q19_IDX;

			//src[I(index, D3Q19_C)]  =[UA] vpdf_C - vomegaEven * (vpdf_C - vw_0 * vdir_indep_trm);
			VSTU(&src[I(index, D3Q19_C)],VSUB(vpdf_C,VMUL(vomegaEven,VSUB(vpdf_C,VMUL(vw_0,vdir_indep_trm)))));

			//vw_1_indep = vw_1 * vdir_indep_trm;
			vw_1_indep = VMUL(vw_1,vdir_indep_trm);

			//vui = vuy;
			vui = vuy;
			//vevenPart = vomegaEven * (VONE_HALF * (vpdf_N + vpdf_S) - vui * vui * vw_1_nine_half - vw_1_indep);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_N,vpdf_S)),VMUL(vui,VMUL(vui,vw_1_nine_half))),vw_1_indep));
			//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_N - vpdf_S) - vui * vw_1_x3);
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_N,vpdf_S)),VMUL(vui,vw_1_x3)));
			//src[ADJ_LIST(D3Q19_N)]  =[UA] vpdf_N - vevenPart - voddPart;
			VSTU(ppdf_S, VSUB(VSUB(vpdf_N,vevenPart),voddPart));
			//src[ADJ_LIST(D3Q19_S)]  =[UA] vpdf_S - vevenPart + voddPart;
			VSTU(ppdf_N, VADD(VSUB(vpdf_S,vevenPart),voddPart));

			//vui = vux;
			vui = vux;
			//vevenPart = vomegaEven * (VONE_HALF * (vpdf_E + vpdf_W) - vui * vui * vw_1_nine_half - vw_1_indep);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_E,vpdf_W)),VMUL(vui,VMUL(vui,vw_1_nine_half))),vw_1_indep));
			//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_E - vpdf_W) - vui * vw_1_x3);
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_E,vpdf_W)),VMUL(vui,vw_1_x3)));
			//src[ADJ_LIST(D3Q19_E)]  =[UA] vpdf_E - vevenPart - voddPart;
			VSTU(ppdf_W, VSUB(VSUB(vpdf_E,vevenPart),voddPart));
			//src[ADJ_LIST(D3Q19_W)]  =[UA] vpdf_W - vevenPart + voddPart;
			VSTU(ppdf_E, VADD(VSUB(vpdf_W,vevenPart),voddPart));

			//vui = vuz;
			vui = vuz;
			//vevenPart = vomegaEven * (VONE_HALF * (vpdf_T + vpdf_B) - vui * vui * vw_1_nine_half - vw_1_indep);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_T,vpdf_B)),VMUL(vui,VMUL(vui,vw_1_nine_half))),vw_1_indep));
			//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_T - vpdf_B) - vui * vw_1_x3);
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_T,vpdf_B)),VMUL(vui,vw_1_x3)));
			//src[ADJ_LIST(D3Q19_T)]  =[UA] vpdf_T - vevenPart - voddPart;
			VSTU(ppdf_B, VSUB(VSUB(vpdf_T,vevenPart),voddPart));
			//src[ADJ_LIST(D3Q19_B)]  =[UA] vpdf_B - vevenPart + voddPart;
			VSTU(ppdf_T, VADD(VSUB(vpdf_B,vevenPart),voddPart));

			//vw_2_indep = vw_2 * vdir_indep_trm;
			vw_2_indep = VMUL(vw_2,vdir_indep_trm);

			//vui = vuy - vux;
			vui = VSUB(vuy,vux);
			//vevenPart = vomegaEven * (VONE_HALF * (vpdf_NW + vpdf_SE) - vui * vui * vw_2_nine_half - vw_2_indep);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_NW,vpdf_SE)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_NW - vpdf_SE) - vui * vw_2_x3);
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_NW,vpdf_SE)),VMUL(vui,vw_2_x3)));
			//src[ADJ_LIST(D3Q19_NW)] =[UA] vpdf_NW - vevenPart - voddPart;
			VSTU(ppdf_SE, VSUB(VSUB(vpdf_NW,vevenPart),voddPart));
			//src[ADJ_LIST(D3Q19_SE)] =[UA] vpdf_SE - vevenPart + voddPart;
			VSTU(ppdf_NW, VADD(VSUB(vpdf_SE,vevenPart),voddPart));

			//vui = vux + vuy;
			vui = VADD(vux,vuy);
			//vevenPart = vomegaEven * (VONE_HALF * (vpdf_NE + vpdf_SW) - vui * vui * vw_2_nine_half - vw_2_indep);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_NE,vpdf_SW)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_NE - vpdf_SW) - vui * vw_2_x3);
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_NE,vpdf_SW)),VMUL(vui,vw_2_x3)));
			//src[ADJ_LIST(D3Q19_NE)] =[UA] vpdf_NE - vevenPart - voddPart;
			VSTU(ppdf_SW, VSUB(VSUB(vpdf_NE,vevenPart),voddPart));
			//src[ADJ_LIST(D3Q19_SW)] =[UA] vpdf_SW - vevenPart + voddPart;
			VSTU(ppdf_NE, VADD(VSUB(vpdf_SW,vevenPart),voddPart));

			//vui = vuz - vux;
			vui = VSUB(vuz,vux);
			//vevenPart = vomegaEven * (VONE_HALF * (vpdf_TW + vpdf_BE) - vui * vui * vw_2_nine_half - vw_2_indep);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TW,vpdf_BE)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_TW - vpdf_BE) - vui * vw_2_x3);
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TW,vpdf_BE)),VMUL(vui,vw_2_x3)));
			//src[ADJ_LIST(D3Q19_TW)] =[UA] vpdf_TW - vevenPart - voddPart;
			VSTU(ppdf_BE, VSUB(VSUB(vpdf_TW,vevenPart),voddPart));
			//src[ADJ_LIST(D3Q19_BE)] =[UA] vpdf_BE - vevenPart + voddPart;
			VSTU(ppdf_TW, VADD(VSUB(vpdf_BE,vevenPart),voddPart));

			//vui = vux + vuz;
			vui = VADD(vux,vuz);
			//vevenPart = vomegaEven * (VONE_HALF * (vpdf_TE + vpdf_BW) - vui * vui * vw_2_nine_half - vw_2_indep);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TE,vpdf_BW)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_TE - vpdf_BW) - vui * vw_2_x3);
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TE,vpdf_BW)),VMUL(vui,vw_2_x3)));
			//src[ADJ_LIST(D3Q19_TE)] =[UA] vpdf_TE - vevenPart - voddPart;
			VSTU(ppdf_BW, VSUB(VSUB(vpdf_TE,vevenPart),voddPart));
			//src[ADJ_LIST(D3Q19_BW)] =[UA] vpdf_BW - vevenPart + voddPart;
			VSTU(ppdf_TE, VADD(VSUB(vpdf_BW,vevenPart),voddPart));

			//vui = vuz - vuy;
			vui = VSUB(vuz,vuy);
			//vevenPart = vomegaEven * (VONE_HALF * (vpdf_TS + vpdf_BN) - vui * vui * vw_2_nine_half - vw_2_indep);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TS,vpdf_BN)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_TS - vpdf_BN) - vui * vw_2_x3);
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TS,vpdf_BN)),VMUL(vui,vw_2_x3)));
			//src[ADJ_LIST(D3Q19_TS)] =[UA] vpdf_TS - vevenPart - voddPart;
			VSTU(ppdf_TN, VSUB(VSUB(vpdf_TS,vevenPart),voddPart));
			//src[ADJ_LIST(D3Q19_BN)] =[UA] vpdf_BN - vevenPart + voddPart;
			VSTU(ppdf_BN, VADD(VSUB(vpdf_BN,vevenPart),voddPart));

			//vui = vuy + vuz;
			vui = VADD(vuy,vuz);
			//vevenPart = vomegaEven * (VONE_HALF * (vpdf_TN + vpdf_BS) - vui * vui * vw_2_nine_half - vw_2_indep);
			vevenPart = VMUL(vomegaEven,VSUB(VSUB(VMUL(VONE_HALF,VADD(vpdf_TN,vpdf_BS)),VMUL(vui,VMUL(vui,vw_2_nine_half))),vw_2_indep));
			//voddPart  = vomegaOdd  * (VONE_HALF * (vpdf_TN - vpdf_BS) - vui * vw_2_x3);
			voddPart = VMUL(vomegaOdd,VSUB(VMUL(VONE_HALF,VSUB(vpdf_TN,vpdf_BS)),VMUL(vui,vw_2_x3)));
			//src[ADJ_LIST(D3Q19_TN)] =[UA] vpdf_TN - vevenPart - voddPart;
			VSTU(ppdf_BS, VSUB(VSUB(vpdf_TN,vevenPart),voddPart));
			//src[ADJ_LIST(D3Q19_BS)] =[UA] vpdf_BS - vevenPart + voddPart;
			VSTU(ppdf_TN, VADD(VSUB(vpdf_BS,vevenPart),voddPart));

			consecValue   -= (VSIZE - 1);
			index         += (VSIZE - 1);
			pointerOffset  = VSIZE;

		}
		else {
			// Scalar part.

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

			adjListIndex = index * N_D3Q19_IDX;

			// direction: w_0
			src[I(index, D3Q19_C)             ]  = pdf_C - omegaEven * (pdf_C - w_0 * dir_indep_trm);

			// direction: w_1
			w_1_indep = w_1 * dir_indep_trm;

			ui = uy;
			evenPart = omegaEven * (0.5 * (pdf_N + pdf_S) - ui * ui * w_1_nine_half - w_1_indep);
			oddPart  = omegaOdd  * (0.5 * (pdf_N - pdf_S) - ui * w_1_x3);
			src[ADJ_LIST(D3Q19_N)]  = pdf_N - evenPart - oddPart;
			src[ADJ_LIST(D3Q19_S)]  = pdf_S - evenPart + oddPart;

			ui = ux;
			evenPart = omegaEven * (0.5 * (pdf_E + pdf_W) - ui * ui * w_1_nine_half - w_1_indep);
			oddPart  = omegaOdd  * (0.5 * (pdf_E - pdf_W) - ui * w_1_x3);
			src[ADJ_LIST(D3Q19_E)]  = pdf_E - evenPart - oddPart;
			src[ADJ_LIST(D3Q19_W)]  = pdf_W - evenPart + oddPart;

			ui = uz;
			evenPart = omegaEven * (0.5 * (pdf_T + pdf_B) - ui * ui * w_1_nine_half - w_1_indep);
			oddPart  = omegaOdd  * (0.5 * (pdf_T - pdf_B) - ui * w_1_x3);
			src[ADJ_LIST(D3Q19_T)]  = pdf_T - evenPart - oddPart;
			src[ADJ_LIST(D3Q19_B)]  = pdf_B - evenPart + oddPart;

			// direction: w_2
			w_2_indep = w_2 * dir_indep_trm;

			ui = -ux + uy;
			evenPart = omegaEven * (0.5 * (pdf_NW + pdf_SE) - ui * ui * w_2_nine_half - w_2_indep);
			oddPart  = omegaOdd  * (0.5 * (pdf_NW - pdf_SE) - ui * w_2_x3);
			src[ADJ_LIST(D3Q19_NW)] = pdf_NW - evenPart - oddPart;
			src[ADJ_LIST(D3Q19_SE)] = pdf_SE - evenPart + oddPart;

			ui = ux + uy;
			evenPart = omegaEven * (0.5 * (pdf_NE + pdf_SW) - ui * ui * w_2_nine_half - w_2_indep);
			oddPart  = omegaOdd  * (0.5 * (pdf_NE - pdf_SW) - ui * w_2_x3);
			src[ADJ_LIST(D3Q19_NE)] = pdf_NE - evenPart - oddPart;
			src[ADJ_LIST(D3Q19_SW)] = pdf_SW - evenPart + oddPart;

			ui = -ux + uz;
			evenPart = omegaEven * (0.5 * (pdf_TW + pdf_BE) - ui * ui * w_2_nine_half - w_2_indep);
			oddPart  = omegaOdd  * (0.5 * (pdf_TW - pdf_BE) - ui * w_2_x3);
			src[ADJ_LIST(D3Q19_TW)] = pdf_TW - evenPart - oddPart;
			src[ADJ_LIST(D3Q19_BE)] = pdf_BE - evenPart + oddPart;

			ui = ux + uz;
			evenPart = omegaEven * (0.5 * (pdf_TE + pdf_BW) - ui * ui * w_2_nine_half - w_2_indep);
			oddPart  = omegaOdd  * (0.5 * (pdf_TE - pdf_BW) - ui * w_2_x3);
			src[ADJ_LIST(D3Q19_TE)] = pdf_TE - evenPart - oddPart;
			src[ADJ_LIST(D3Q19_BW)] = pdf_BW - evenPart + oddPart;

			ui = -uy + uz;
			evenPart = omegaEven * (0.5 * (pdf_TS + pdf_BN) - ui * ui * w_2_nine_half - w_2_indep);
			oddPart  = omegaOdd  * (0.5 * (pdf_TS - pdf_BN) - ui * w_2_x3);
			src[ADJ_LIST(D3Q19_TS)] = pdf_TS - evenPart - oddPart;
			src[ADJ_LIST(D3Q19_BN)] = pdf_BN - evenPart + oddPart;

			ui = uy + uz;
			evenPart = omegaEven * (0.5 * (pdf_TN + pdf_BS) - ui * ui * w_2_nine_half - w_2_indep);
			oddPart  = omegaOdd  * (0.5 * (pdf_TN - pdf_BS) - ui * w_2_x3);
			src[ADJ_LIST(D3Q19_TN)] = pdf_TN - evenPart - oddPart;
			src[ADJ_LIST(D3Q19_BS)] = pdf_BS - evenPart + oddPart;

			pointerOffset = 1;
		}

	} // loop over fluid nodes

	#undef ADJ_LIST
	#undef I
}
