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
#ifndef __BENCH_KERNEL_D3Q19_COMMON_H__
#define __BENCH_KERNEL_D3Q19_COMMON_H__

#include "Kernel.h"

typedef struct KernelDataEx_
{
	KernelData kd;
	int Blk[3];			// Blocking in X, Y, and Z direction, value of 0 disables blocking.
} KernelDataEx;


// Build a function name extended by the propagation model name and the data layout.
// FNANEM(test) will be expanded to test_PushSoA if DATA_LAYOUT_NAME is defined
// as SoA and PROP_MODEL is defined as Push.
#define FNAME(functionName)	JOIN(JOIN(functionName,_),JOIN(PROP_MODEL_NAME,DATA_LAYOUT_NAME))

#ifndef DATA_LAYOUT_NAME
	#error DATA_LAYOUT_NAME must be defined
#endif

#ifndef PROP_MODEL_NAME
	#error PROP_MODEL_NAME must be defined
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
	return d * dims[0] * dims[1] * dims[2] + z * dims[0] * dims[1] + y * dims[0] + x;
#elif DATA_LAYOUT_AOS
	return z * dims[0] * dims[1] * N_D3Q19 + y * dims[0] * N_D3Q19 + x * N_D3Q19 + d;
#else
#error P_INDEX_5 function no implemented for chosen data layout.
#endif
}

#endif // __BENCH_KERNEL_D3Q19_COMMON_H__

