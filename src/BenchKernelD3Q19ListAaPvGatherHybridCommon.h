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
#ifndef __BENCH_KERNEL_D3Q19_LIST_AA_PV_GATHER_HYBRID_COMMON_H__
#define __BENCH_KERNEL_D3Q19_LIST_AA_PV_GATHER_HYBRID_COMMON_H__

#if !defined(DATA_LAYOUT_SOA) && !defined(DATA_LAYOUT_AOSOA)
	#error List Gather Hybrid works only with DATA_LAYOUT_SOA or DATA_LAYOUT_AOSOA
#endif

#include "BenchKernelD3Q19ListAaCommon.h"

typedef struct KernelDataListRia_ {
	KernelDataList kdl;

	// Array contains information for loop start indices with the following scheme for every thread:
	// scalar peel start | (vectorized load store | vectorized gather scatter) ... | scalar remainder.
	// Example for 3 threads with ! indicating thread boundaries 
	// [sp,vls,vgs,...,vls,vgs,sr ! sp, vls, vgs, ..., vls, vgs, sr ! sp, vls, vgs, ..., vls, vgs, sr]
	int * loopStartIndices;
	int nLoopStartIndices;  // Number of entries in loopStartIndices array.

	// Array contains (for each thread) an index into loopStartIndices.
	int * oddKernelThreadStartIndices;
	// Number of entries in threadStartIndices
	int nOddKernelThreadStartIndices;

} KernelDataListRia;

// Macro for casting KernelData * to KernelDataList *.
#define KDLR(_x_)	((KernelDataListRia *)(_x_))


#endif // __BENCH_KERNEL_D3Q19_LIST_AA_PV_GATHER_HYBRID_COMMON_H__

