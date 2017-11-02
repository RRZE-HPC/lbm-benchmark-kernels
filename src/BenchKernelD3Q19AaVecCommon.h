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
#ifndef __BENCH_KERNEL_D3Q19_AA_VEC_COMMON_H__
#define __BENCH_KERNEL_D3Q19_AA_VEC_COMMON_H__

// #include "BenchKernelD3Q19AaCommon.h"

#include "Kernel.h"

typedef struct KernelDataAa_
{
	KernelData kd;
	int Blk[3];			// Blocking in X, Y, and Z direction, value of 0 disables blocking.
	int Iteration;		// Current iteration number.
} KernelDataAa;

// Macro for casting KernelData * to KernelDataAa *.
#define KDA(_x_)	((KernelDataAa *)(_x_))

// Build a function name extended by the propagation model name and the data layout.
// FNANEM(test) will be expanded to test_PushSoA if DATA_LAYOUT_NAME is defined
// as SoA and PROP_MODEL is defined as Push.
#define FNAME(functionName)	JOIN(JOIN(functionName,_),JOIN(PROP_MODEL_NAME,DATA_LAYOUT_NAME))

#ifndef DATA_LAYOUT_SOA
	#error Only DATA_LAYOUT_SOA is supported.
#endif

#ifndef PROP_MODEL_AA
	#error Only PROP_MODEL_AA is supported.
#endif

// -----------------------------------------------------------------------
// Index function for accesssing PDF array for different data layouts.

#define P_INDEX_5		FNAME(PINDEX5)

static inline int FNAME(PINDEX5)(int dims[3], int x, int y, int z, int d)
{
	Assert(dims[0] > 0);
	Assert(dims[1] > 0);
	Assert(dims[2] > 0);

	Assert(x >= 0);
	Assert(x < dims[0]);
	Assert(y >= 0);
	Assert(y < dims[1]);
	Assert(z >= 0);
	Assert(z < dims[2]);
	Assert(d >= 0);
#ifdef D3Q19
	Assert(d < N_D3Q19);
#else
	#error Not implemented for this discretization.
#endif

#ifdef DATA_LAYOUT_SOA
	return d * dims[0] * dims[1] * dims[2] + x * dims[1] * dims[2] + y * dims[2] + z;
// #elif DATA_LAYOUT_AOS
// 	return x * dims[1] * dims[2] * N_D3Q19 + y * dims[2] * N_D3Q19 + z * N_D3Q19 + d;
#else
	#error P_INDEX_5 function no implemented for chosen data layout.
#endif
}



#endif // __BENCH_KERNEL_D3Q19_AA_VEC_COMMON_H__

