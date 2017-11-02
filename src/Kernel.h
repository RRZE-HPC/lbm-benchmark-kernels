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
#ifndef __KERNEL_H__
#define __KERNEL_H__

#include "Base.h"
#include "Lattice.h"

#ifdef DATA_LAYOUT_NAME
	#error DATA_LAYOUT_NAME must not be defined here.
#endif

#ifdef PROP_MODEL_NAME
	#error PROP_MODEL_NAME must not be defined here.
#endif


#ifdef DATA_LAYOUT_SOA
	#define DATA_LAYOUT_NAME SoA
#endif

#ifdef DATA_LAYOUT_AOS
	#define DATA_LAYOUT_NAME AoS
#endif

#ifdef PROP_MODEL_PUSH
	#define PROP_MODEL_NAME Push
#endif

#ifdef PROP_MODEL_PULL
	#define PROP_MODEL_NAME Pull
#endif

#ifdef PROP_MODEL_AA
	#define PROP_MODEL_NAME Aa
#endif



typedef double 	PdfT;



#define D3Q19


#define N_D3Q19		19

#define D3Q19_N		0
#define D3Q19_S		1
#define D3Q19_E		2
#define D3Q19_W		3

#define D3Q19_NE	4
#define D3Q19_SE	5
#define D3Q19_NW	6
#define D3Q19_SW	7

#define D3Q19_T		8
#define D3Q19_TN	9
#define D3Q19_TE	10
#define D3Q19_TW	11
#define D3Q19_TS	12

#define D3Q19_B		13
#define D3Q19_BS	14
#define D3Q19_BN	15
#define D3Q19_BW	16
#define D3Q19_BE	17

#define D3Q19_C		18		// IMPORTANT: Center particle must be the last one.

// ---

#ifdef X
	#error X is not allowed to be defined here!
#endif

// The following list must be sorted ascending according
// to the index of the direction, i.e. D3Q19_N, D3Q19_S, ...
#define D3Q19_LIST \
	X(N,  D3Q19_N, 	D3Q19_S, 	  0,  1,  0) \
	X(S,  D3Q19_S,	D3Q19_N, 	  0, -1,  0) \
	X(E,  D3Q19_E,	D3Q19_W, 	  1,  0,  0) \
	X(W,  D3Q19_W,	D3Q19_E, 	 -1,  0,  0) \
	X(NE, D3Q19_NE,	D3Q19_SW, 	  1,  1,  0) \
	X(SE, D3Q19_SE,	D3Q19_NW,	  1, -1,  0) \
	X(NW, D3Q19_NW,	D3Q19_SE, 	 -1,  1,  0) \
	X(SW, D3Q19_SW,	D3Q19_NE, 	 -1, -1,  0) \
	X(T,  D3Q19_T,	D3Q19_B, 	  0,  0,  1) \
	X(TN, D3Q19_TN,	D3Q19_BS,	  0,  1,  1) \
	X(TE, D3Q19_TE,	D3Q19_BW,     1,  0,  1) \
	X(TW, D3Q19_TW,	D3Q19_BE,    -1,  0,  1) \
	X(TS, D3Q19_TS,	D3Q19_BN,     0, -1,  1) \
	X(B,  D3Q19_B,	D3Q19_T,	  0,  0, -1) \
	X(BS, D3Q19_BS,	D3Q19_TN,     0, -1, -1) \
	X(BN, D3Q19_BN,	D3Q19_TS,     0,  1, -1) \
	X(BW, D3Q19_BW,	D3Q19_TE,    -1,  0, -1) \
	X(BE, D3Q19_BE,	D3Q19_TW,     1,  0, -1) \
	X(C,  D3Q19_C,	D3Q19_C,	  0,  0,  0)

#define D3Q19_LIST_WO_C \
	X(N,  D3Q19_N, 	D3Q19_S, 	  0,  1,  0) \
	X(S,  D3Q19_S,	D3Q19_N, 	  0, -1,  0) \
	X(E,  D3Q19_E,	D3Q19_W, 	  1,  0,  0) \
	X(W,  D3Q19_W,	D3Q19_E, 	 -1,  0,  0) \
	X(NE, D3Q19_NE,	D3Q19_SW, 	  1,  1,  0) \
	X(SE, D3Q19_SE,	D3Q19_NW,	  1, -1,  0) \
	X(NW, D3Q19_NW,	D3Q19_SE, 	 -1,  1,  0) \
	X(SW, D3Q19_SW,	D3Q19_NE, 	 -1, -1,  0) \
	X(T,  D3Q19_T,	D3Q19_B, 	  0,  0,  1) \
	X(TN, D3Q19_TN,	D3Q19_BS,	  0,  1,  1) \
	X(TE, D3Q19_TE,	D3Q19_BW,     1,  0,  1) \
	X(TW, D3Q19_TW,	D3Q19_BE,    -1,  0,  1) \
	X(TS, D3Q19_TS,	D3Q19_BN,     0, -1,  1) \
	X(B,  D3Q19_B,	D3Q19_T,	  0,  0, -1) \
	X(BS, D3Q19_BS,	D3Q19_TN,     0, -1, -1) \
	X(BN, D3Q19_BN,	D3Q19_TS,     0,  1, -1) \
	X(BW, D3Q19_BW,	D3Q19_TE,    -1,  0, -1) \
	X(BE, D3Q19_BE,	D3Q19_TW,     1,  0, -1)


extern int D3Q19_X[N_D3Q19];
extern int D3Q19_Y[N_D3Q19];
extern int D3Q19_Z[N_D3Q19];
extern int D3Q19_INV[N_D3Q19];

extern const char * D3Q19_NAMES[N_D3Q19];



typedef struct CaseData_ {
	PdfT Omega;
	PdfT RhoIn;
	PdfT RhoOut;
	PdfT Ux;
	PdfT Uy;
	PdfT Uz;
	PdfT XForce;
	int MaxIterations;
	int VtkOutput;
	int VtkModulus;
	int StatisticsModulus;
} CaseData;


typedef struct KernelData_ {
	PdfT * Pdfs[2];
	PdfT * SrcPdfs;
	PdfT * DstPdfs;
	PdfT * PdfsActive;
	int Dims[3];
	int GlobalDims[3];
	int Offsets[3];
	int * ObstIndices;
	int nObstIndices;
	int * BounceBackPdfsSrc;
	int * BounceBackPdfsDst;
	int nBounceBackPdfs;

	void (* BoundaryConditionsGetPdf)(struct KernelData_ * kd, int x, int y, int z, int dir, PdfT * pdf);
	void (* BoundaryConditionsSetPdf)(struct KernelData_ * kd, int x, int y, int z, int dir, PdfT pdf);

	void (* GetNode)(struct KernelData_ * kd, int x, int y, int z, PdfT * pdfs);
	void (* SetNode)(struct KernelData_ * kd, int x, int y, int z, PdfT * pdfs);

	void (* Kernel)(LatticeDesc * ld, struct KernelData_ * kd, CaseData * cd);

} KernelData;

typedef struct Parameters_ {
	int nArgs;
	char ** Args;
	int nKernelArgs;
	char ** KernelArgs;
} Parameters;

void KernelComputeBoundaryConditions(KernelData * kd, LatticeDesc * ld, CaseData * cd);

PdfT KernelDensity(KernelData * kd, LatticeDesc * ld);

void KernelStatistics(KernelData * kd, LatticeDesc * ld, CaseData * cd, int iteration);
void KernelStatisticsAdv(KernelData * kd, LatticeDesc * ld, CaseData * cd, int iteration, int forceOutput);


void KernelSetInitialDensity (LatticeDesc * ld, KernelData * kd, CaseData * cd);
void KernelSetInitialVelocity(LatticeDesc * ld, KernelData * kd, CaseData * cd);

void KernelVerifiy(LatticeDesc * ld, KernelData * kd, CaseData * cd, PdfT * errorNorm);

void KernelAddBodyForce(KernelData * kd, LatticeDesc * ld, CaseData * cd);

#endif // __KERNEL_H__
