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
#include "BenchKernelD3Q19ListAaCommon.h"

#include "Memory.h"
#include "Vtk.h"
#include "LikwidIf.h"

#include <inttypes.h>
#include <math.h>

void FNAME(D3Q19ListAaKernel)(LatticeDesc * ld, KernelData * kernelData, CaseData * cd)
{
	Assert(ld != NULL);
	Assert(kernelData != NULL);
	Assert(cd != NULL);

	Assert(cd->Omega > F(0.0));
	Assert(cd->Omega < F(2.0));

	KernelData * kd = (KernelData *)kernelData;
	KernelDataList * kdl = (KernelDataList *)kernelData;

	PdfT omega = cd->Omega;
	PdfT omegaEven = omega;
	//  1/4: best stability;
	// 1/12: removes third-order advection error (best advection);
	//  1/6: removes fourth-order diffusion error (best diffusion);
	// 3/16: exact location of bounce back for poiseuille flow
	PdfT magicParam = F(1.0)/F(12.0);
	PdfT omegaOdd = F(1.0)/( F(0.5) + magicParam/(F(1.0)/omega - F(0.5)) );

	PdfT evenPart = F(0.0);
	PdfT oddPart = F(0.0);
	PdfT dir_indep_trm = F(0.0);

	const PdfT w_0 = F(1.0) /  F(3.0);
	const PdfT w_1 = F(1.0) / F(18.0);
	const PdfT w_2 = F(1.0) / F(36.0);

	const PdfT w_1_x3 = w_1 * F(3.0);	const PdfT w_1_nine_half = w_1 * F(9.0)/F(2.0);	PdfT w_1_indep = F(0.0);
	const PdfT w_2_x3 = w_2 * F(3.0);	const PdfT w_2_nine_half = w_2 * F(9.0)/F(2.0);	PdfT w_2_indep = F(0.0);

	PdfT ui;

	PdfT ux, uy, uz;
	PdfT dens;

	// Declare pdf_N, pdf_E, pdf_S, pdf_W, ...
	#define X(name, idx, idxinv, x, y, z)	PdfT JOIN(pdf_,name);
	D3Q19_LIST
	#undef X

	PdfT * src = kd->Pdfs[0];

	int maxIterations = cd->MaxIterations;
	int nFluid = kdl->nFluid;
	int nCells = kdl->nCells;

	uint32_t adjListIndex;
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

	// TODO: outer openmp parallel
	for(int iter = 0; iter < maxIterations; iter += 2) {


		// --------------------------------------------------------------------
		// even time step

		X_LIKWID_START("list-aa-even");

#ifdef _OPENMP
		#pragma omp parallel for default(none) \
				shared(nFluid, nCells, kd, kdl, adjList, omegaOdd, omegaEven, src) \
				private(ux, uy, uz, dens, adjListIndex, evenPart, oddPart, dir_indep_trm, w_1_indep, w_2_indep, ui,\
					pdf_C, \
				  	pdf_N, pdf_E, pdf_S, pdf_W, \
				  	pdf_NE, pdf_SE, pdf_SW, pdf_NW, \
				  	pdf_T, pdf_TN, pdf_TE, pdf_TS, pdf_TW, \
				  	pdf_B, pdf_BN, pdf_BE, pdf_BS, pdf_BW)
#endif
#ifdef INTEL_OPT_DIRECTIVES
		#pragma ivdep
		#pragma vector always
		#pragma simd
#endif
		for (int index = 0; index < nFluid; ++index) {


			#define I(index, dir)	P_INDEX_3((nCells), (index), (dir))

			// Load PDFs of local cell: pdf_N = src[I(x, y, z, D3Q19_N)]; ...
			#define X(name, idx, idxinv, _x, _y, _z)	JOIN(pdf_,name) = src[I(index, idx)];
			D3Q19_LIST
			#undef X

// #define LID_DRIVEN_CAVITY

#ifdef LID_DRIVEN_CAVITY
			int nX = kd->Dims[0];
			int nY = kd->Dims[1];
			int nZ = kd->Dims[2];

			int x = kdl->Coords[C_INDEX_X(index)];
			int y = kdl->Coords[C_INDEX_Y(index)];
			int z = kdl->Coords[C_INDEX_Z(index)];

			if (z == nZ - 4 && x > 3 && x < (nX - 4) && y > 3 && y < (nY - 4)) {
				ux = 0.1 * 0.577;
				uy = 0.0;
				uz = 0.0;
			} else {
#endif
				ux = pdf_E + pdf_NE + pdf_SE + pdf_TE + pdf_BE -
					 pdf_W - pdf_NW - pdf_SW - pdf_TW - pdf_BW;
				uy = pdf_N + pdf_NE + pdf_NW + pdf_TN + pdf_BN -
					 pdf_S - pdf_SE - pdf_SW - pdf_TS - pdf_BS;
				uz = pdf_T + pdf_TE + pdf_TW + pdf_TN + pdf_TS -
					 pdf_B - pdf_BE - pdf_BW - pdf_BN - pdf_BS;
#ifdef LID_DRIVEN_CAVITY
			}
#endif

			dens = pdf_C +
				   pdf_N  + pdf_E  + pdf_S  + pdf_W  +
				   pdf_NE + pdf_SE + pdf_SW + pdf_NW +
				   pdf_T  + pdf_TN + pdf_TE + pdf_TS + pdf_TW +
				   pdf_B  + pdf_BN + pdf_BE + pdf_BS + pdf_BW;

			dir_indep_trm = dens - (ux * ux + uy * uy + uz * uz)*F(3.0)/F(2.0);

			// direction: w_0
			src[I(index, D3Q19_C)             ]  = pdf_C - omegaEven*(pdf_C - w_0*dir_indep_trm);

			// direction: w_1
			w_1_indep = w_1*dir_indep_trm;

			ui = uy;
			evenPart = omegaEven*( F(0.5)*(pdf_N + pdf_S) - ui*ui*w_1_nine_half - w_1_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_N - pdf_S) - ui*w_1_x3 );
			src[I(index, D3Q19_S)]  = pdf_N - evenPart - oddPart;
			src[I(index, D3Q19_N)]  = pdf_S - evenPart + oddPart;

			ui = ux;
			evenPart = omegaEven*( F(0.5)*(pdf_E + pdf_W) - ui*ui*w_1_nine_half - w_1_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_E - pdf_W) - ui*w_1_x3 );
			src[I(index, D3Q19_W)]  = pdf_E - evenPart - oddPart;
			src[I(index, D3Q19_E)]  = pdf_W - evenPart + oddPart;

			ui = uz;
			evenPart = omegaEven*( F(0.5)*(pdf_T + pdf_B) - ui*ui*w_1_nine_half - w_1_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_T - pdf_B) - ui*w_1_x3 );
			src[I(index, D3Q19_B)]  = pdf_T - evenPart - oddPart;
			src[I(index, D3Q19_T)]  = pdf_B - evenPart + oddPart;

			// direction: w_2
			w_2_indep = w_2*dir_indep_trm;

			ui = -ux + uy;
			evenPart = omegaEven*( F(0.5)*(pdf_NW + pdf_SE) - ui*ui*w_2_nine_half - w_2_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_NW - pdf_SE) - ui*w_2_x3 );
			src[I(index, D3Q19_SE)] = pdf_NW - evenPart - oddPart;
			src[I(index, D3Q19_NW)] = pdf_SE - evenPart + oddPart;

			ui = ux + uy;
			evenPart = omegaEven*( F(0.5)*(pdf_NE + pdf_SW) - ui*ui*w_2_nine_half - w_2_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_NE - pdf_SW) - ui*w_2_x3 );
			src[I(index, D3Q19_SW)] = pdf_NE - evenPart - oddPart;
			src[I(index, D3Q19_NE)] = pdf_SW - evenPart + oddPart;

			ui = -ux + uz;
			evenPart = omegaEven*( F(0.5)*(pdf_TW + pdf_BE) - ui*ui*w_2_nine_half - w_2_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_TW - pdf_BE) - ui*w_2_x3 );
			src[I(index, D3Q19_BE)] = pdf_TW - evenPart - oddPart;
			src[I(index, D3Q19_TW)] = pdf_BE - evenPart + oddPart;

			ui = ux + uz;
			evenPart = omegaEven*( F(0.5)*(pdf_TE + pdf_BW) - ui*ui*w_2_nine_half - w_2_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_TE - pdf_BW) - ui*w_2_x3 );
			src[I(index, D3Q19_BW)] = pdf_TE - evenPart - oddPart;
			src[I(index, D3Q19_TE)] = pdf_BW - evenPart + oddPart;

			ui = -uy + uz;
			evenPart = omegaEven*( F(0.5)*(pdf_TS + pdf_BN) - ui*ui*w_2_nine_half - w_2_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_TS - pdf_BN) - ui*w_2_x3 );
			src[I(index, D3Q19_BN)] = pdf_TS - evenPart - oddPart;
			src[I(index, D3Q19_TS)] = pdf_BN - evenPart + oddPart;

			ui = uy + uz;
			evenPart = omegaEven*( F(0.5)*(pdf_TN + pdf_BS) - ui*ui*w_2_nine_half - w_2_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_TN - pdf_BS) - ui*w_2_x3 );
			src[I(index, D3Q19_BS)] = pdf_TN - evenPart - oddPart;
			src[I(index, D3Q19_TN)] = pdf_BS - evenPart + oddPart;

		} // loop over fluid nodes

		X_LIKWID_STOP("list-aa-even");

		// save current iteration
		kdl->Iteration = iter;
		#ifdef VERIFICATION
			kd->PdfsActive = src;
			KernelAddBodyForce(kd, ld, cd);
		#endif

		// --------------------------------------------------------------------
		// odd time step

		X_LIKWID_START("list-aa-odd");

#ifdef _OPENMP
		#pragma omp parallel for default(none) \
				shared(nFluid, nCells, kd, kdl, adjList, omegaOdd, omegaEven, src) \
				private(ux, uy, uz, dens, adjListIndex, evenPart, oddPart, dir_indep_trm, w_1_indep, w_2_indep, ui,\
					pdf_C, \
				  	pdf_N, pdf_E, pdf_S, pdf_W, \
				  	pdf_NE, pdf_SE, pdf_SW, pdf_NW, \
				  	pdf_T, pdf_TN, pdf_TE, pdf_TS, pdf_TW, \
				  	pdf_B, pdf_BN, pdf_BE, pdf_BS, pdf_BW)
#endif
#ifdef INTEL_OPT_DIRECTIVES
		#pragma ivdep
#endif
		for (int index = 0; index < nFluid; ++index) {


			adjListIndex = index * N_D3Q19_IDX;

			// Load PDFs of local cell: pdf_N = src[adjList[adjListIndex + D3Q19_S]]; ...
			pdf_C = src[P_INDEX_3(nCells, index, D3Q19_C)];

			#define X(name, idx, idxinv, _x, _y, _z)	JOIN(pdf_,name) = src[adjList[adjListIndex + idxinv]];
			D3Q19_LIST_WO_C
			#undef X

#ifdef LID_DRIVEN_CAVITY
			int nX = kd->Dims[0];
			int nY = kd->Dims[1];
			int nZ = kd->Dims[2];

			int x = kdl->Coords[C_INDEX_X(index)];
			int y = kdl->Coords[C_INDEX_Y(index)];
			int z = kdl->Coords[C_INDEX_Z(index)];

			if (z == nZ - 4 && x > 3 && x < (nX - 4) && y > 3 && y < (nY - 4)) {
				ux = F(0.1) * F(0.5)77;
				uy = F(0.0);
				uz = F(0.0);
			} else {
#endif
				ux = pdf_E + pdf_NE + pdf_SE + pdf_TE + pdf_BE -
					 pdf_W - pdf_NW - pdf_SW - pdf_TW - pdf_BW;
				uy = pdf_N + pdf_NE + pdf_NW + pdf_TN + pdf_BN -
					 pdf_S - pdf_SE - pdf_SW - pdf_TS - pdf_BS;
				uz = pdf_T + pdf_TE + pdf_TW + pdf_TN + pdf_TS -
					 pdf_B - pdf_BE - pdf_BW - pdf_BN - pdf_BS;
#ifdef LID_DRIVEN_CAVITY
			}
#endif

			dens = pdf_C +
				   pdf_N  + pdf_E  + pdf_S  + pdf_W  +
				   pdf_NE + pdf_SE + pdf_SW + pdf_NW +
				   pdf_T  + pdf_TN + pdf_TE + pdf_TS + pdf_TW +
				   pdf_B  + pdf_BN + pdf_BE + pdf_BS + pdf_BW;

			dir_indep_trm = dens - (ux * ux + uy * uy + uz * uz)*F(3.0)/F(2.0);

			adjListIndex = index * N_D3Q19_IDX;

			// direction: w_0
			src[I(index, D3Q19_C)             ]  = pdf_C - omegaEven*(pdf_C - w_0*dir_indep_trm);

			// direction: w_1
			w_1_indep = w_1*dir_indep_trm;

			ui = uy;
			evenPart = omegaEven*( F(0.5)*(pdf_N + pdf_S) - ui*ui*w_1_nine_half - w_1_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_N - pdf_S) - ui*w_1_x3 );
			src[adjList[adjListIndex + D3Q19_N]]  = pdf_N - evenPart - oddPart;
			src[adjList[adjListIndex + D3Q19_S]]  = pdf_S - evenPart + oddPart;

			ui = ux;
			evenPart = omegaEven*( F(0.5)*(pdf_E + pdf_W) - ui*ui*w_1_nine_half - w_1_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_E - pdf_W) - ui*w_1_x3 );
			src[adjList[adjListIndex + D3Q19_E]]  = pdf_E - evenPart - oddPart;
			src[adjList[adjListIndex + D3Q19_W]]  = pdf_W - evenPart + oddPart;

			ui = uz;
			evenPart = omegaEven*( F(0.5)*(pdf_T + pdf_B) - ui*ui*w_1_nine_half - w_1_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_T - pdf_B) - ui*w_1_x3 );
			src[adjList[adjListIndex + D3Q19_T]]  = pdf_T - evenPart - oddPart;
			src[adjList[adjListIndex + D3Q19_B]]  = pdf_B - evenPart + oddPart;

			// direction: w_2
			w_2_indep = w_2*dir_indep_trm;

			ui = -ux + uy;
			evenPart = omegaEven*( F(0.5)*(pdf_NW + pdf_SE) - ui*ui*w_2_nine_half - w_2_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_NW - pdf_SE) - ui*w_2_x3 );
			src[adjList[adjListIndex + D3Q19_NW]] = pdf_NW - evenPart - oddPart;
			src[adjList[adjListIndex + D3Q19_SE]] = pdf_SE - evenPart + oddPart;

			ui = ux + uy;
			evenPart = omegaEven*( F(0.5)*(pdf_NE + pdf_SW) - ui*ui*w_2_nine_half - w_2_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_NE - pdf_SW) - ui*w_2_x3 );
			src[adjList[adjListIndex + D3Q19_NE]] = pdf_NE - evenPart - oddPart;
			src[adjList[adjListIndex + D3Q19_SW]] = pdf_SW - evenPart + oddPart;

			ui = -ux + uz;
			evenPart = omegaEven*( F(0.5)*(pdf_TW + pdf_BE) - ui*ui*w_2_nine_half - w_2_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_TW - pdf_BE) - ui*w_2_x3 );
			src[adjList[adjListIndex + D3Q19_TW]] = pdf_TW - evenPart - oddPart;
			src[adjList[adjListIndex + D3Q19_BE]] = pdf_BE - evenPart + oddPart;

			ui = ux + uz;
			evenPart = omegaEven*( F(0.5)*(pdf_TE + pdf_BW) - ui*ui*w_2_nine_half - w_2_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_TE - pdf_BW) - ui*w_2_x3 );
			src[adjList[adjListIndex + D3Q19_TE]] = pdf_TE - evenPart - oddPart;
			src[adjList[adjListIndex + D3Q19_BW]] = pdf_BW - evenPart + oddPart;

			ui = -uy + uz;
			evenPart = omegaEven*( F(0.5)*(pdf_TS + pdf_BN) - ui*ui*w_2_nine_half - w_2_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_TS - pdf_BN) - ui*w_2_x3 );
			src[adjList[adjListIndex + D3Q19_TS]] = pdf_TS - evenPart - oddPart;
			src[adjList[adjListIndex + D3Q19_BN]] = pdf_BN - evenPart + oddPart;

			ui = uy + uz;
			evenPart = omegaEven*( F(0.5)*(pdf_TN + pdf_BS) - ui*ui*w_2_nine_half - w_2_indep );
			oddPart = omegaOdd*(F(0.5)*(pdf_TN - pdf_BS) - ui*w_2_x3 );
			src[adjList[adjListIndex + D3Q19_TN]] = pdf_TN - evenPart - oddPart;
			src[adjList[adjListIndex + D3Q19_BS]] = pdf_BS - evenPart + oddPart;

			#undef I
		} // loop over fluid nodes

		X_LIKWID_STOP("list-aa-odd");

		// save current iteration
		kdl->Iteration = iter + 1;

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


	} // for (int iter = 0; ...


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
