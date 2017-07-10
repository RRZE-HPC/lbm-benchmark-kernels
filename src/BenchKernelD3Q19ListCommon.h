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
#ifndef __BENCH_KERNEL_D3Q19_LIST_COMMON_H__
#define __BENCH_KERNEL_D3Q19_LIST_COMMON_H__


#include "Kernel.h"

#include <inttypes.h>

#define N_D3Q19_IDX		18

typedef struct KernelDataList_
{
	KernelData kd;
	uint32_t * AdjList;	// Stores PDF indices, which are the destination for propagation.
						// Determine the destination for node index n and direction d via:
						// (n * N_D3Q19_IDX) + d
	uint32_t * Grid;	// Stores the node indices; use L_INDEX_4 macro for access.
	uint32_t * Coords;	// Map node indices to coordiantes; use C_INDEX_* macro for access.
	int nFluid;			// Number of fluid nodes allocated, i.e. length of adjList * N_D3Q19_IDX.
	int nCells;			// Total number of nodes allocated, including nodes for padding!
} KernelDataList;


// Macro for casting KernelData * to KernelDataList *.
#define KDL(_x_)	((KernelDataList *)(_x_))




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

#define P_INDEX_3		FNAME(PINDEX3)

static inline int FNAME(PINDEX3)(int nCells, int cellIndex, int d)
{
	Assert(nCells > 0);
	Assert(cellIndex >= 0);
	Assert(cellIndex < nCells);

	Assert(d >= 0);
	#ifdef D3Q19
		Assert(d < N_D3Q19);
	#else
		#error Not implemented for this discretization.
	#endif

#ifdef DATA_LAYOUT_SOA
	return d * nCells + cellIndex;
#elif  DATA_LAYOUT_AOS
	return cellIndex * N_D3Q19 + d;
#else
	#error P_INDEX_3 function not implemented for chosen data layout.
#endif
}

#define P_INDEX_5		FNAME(PINDEX5)

static inline int FNAME(PINDEX5)(KernelDataList * kdl, int x, int y, int z, int d)
{
	Assert(kdl != NULL);
#ifdef DEBUG
	uint32_t * grid = kdl->Grid;
	int * dims = kdl->kd.Dims;

	Assert(grid != NULL);
	Assert(dims != NULL);
	Assert(dims[0] > 0);
	Assert(dims[1] > 0);
	Assert(dims[2] > 0);
	Assert(x >= 0 && x < dims[0]);
	Assert(y >= 0 && y < dims[1]);
	Assert(z >= 0 && z < dims[2]);
	Assert(d >= 0 && d < N_D3Q19);
#endif

	return P_INDEX_3(kdl->nCells, kdl->Grid[L_INDEX_4(kdl->kd.Dims, x, y, z)], d);
}

// -----------------------------------------------------------------------
// Macros for accessing coord array

#define C_INDEX_X(cellIndex)	C_INDEX(cellIndex, 0)
#define C_INDEX_Y(cellIndex)	C_INDEX(cellIndex, 1)
#define C_INDEX_Z(cellIndex)	C_INDEX(cellIndex, 2)

static inline int C_INDEX(int cellIndex, int xyz)
{
	Assert(cellIndex >= 0);
	Assert(xyz >= 0);
	Assert(xyz < 3);

	return cellIndex * 3 + xyz;
}


#endif // __BENCH_KERNEL_D3Q19_LIST_COMMON_H__
