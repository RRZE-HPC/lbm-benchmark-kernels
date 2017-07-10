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
#ifndef __BENCH_KERNEL_D3Q19_LIST_AA_RIA_COMMON_H__
#define __BENCH_KERNEL_D3Q19_LIST_AA_RIA_COMMON_H__

#if !defined(DATA_LAYOUT_SOA)
	#error List AA Ria works only with DATA_LAYOUT_SOA
#endif

#include "BenchKernelD3Q19ListAaCommon.h"

typedef struct KernelDataListRia_ {
	KernelDataList kdl;

	// Array contains information of how many adjacent nodes share the same access pattern.
	uint32_t * ConsecNodes;
	uint32_t nConsecNodes;  // Number of entries in ConsecNodes array.

	// Array contains (for each thread) an index into ConsecNodes.
	uint32_t * ConsecThreadIndices;
	// Number of entries in ConsecThreadIndices.
	uint32_t nConsecThreadIndices;

	// Array contains fluid node indices for each thread where to start in the
	// vector of fluid nodes.
	int * FluidNodeThreadIndices;
	// Number of entries in FluidNodeThreadIndices.
	int nFluidNodeThreadIndices;

} KernelDataListRia;

// Macro for casting KernelData * to KernelDataList *.
#define KDLR(_x_)	((KernelDataListRia *)(_x_))


#endif // __BENCH_KERNEL_D3Q19_LIST_AA_RIA_COMMON_H__

