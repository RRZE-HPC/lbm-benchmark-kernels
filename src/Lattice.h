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
#ifndef __LATTICE_H__
#define __LATTICE_H__

#include "Base.h"


typedef int		LatticeT;

typedef enum LAT_CELL_TYPES_ {
	LAT_CELL_OBSTACLE 	= 0,
	LAT_CELL_FLUID 		= 1,
	LAT_CELL_INLET 		= 2,
	LAT_CELL_OUTLET 	= 4
} LAT_CELL_TYPES;


typedef struct LatticeDesc_ {
	int Dims[3];
	LatticeT * Lattice;
	int nCells;			// Total number of cells (Obstacles + Fluids).
	int nFluid;			// Total number of fluid cells. Fluid cells are fluid, inlet, and outlet.
	int nObst;			// Total number of obstacles in the lattice.
	int nInlet;			// Total number of inlets in the lattice.
	int nOutlet;		// Total number of outlets.
	int PeriodicX;		// Periodic in X direction.
	int PeriodicY;		// Periodic in Y direction.
	int PeriodicZ;		// Periodic in Z direction.
} LatticeDesc;


// #define L_INDEX_4(dims, x, y, z)	((z) * (dims[0]) * (dims[1]) + (y) * (dims[0]) + (x))

static inline int L_INDEX_4(int dims[3], int x, int y, int z)
{
	Assert(dims != NULL);

	Assert(dims[0] > 0);
	Assert(dims[1] > 0);
	Assert(dims[2] > 0);

	Assert(x >= 0);
	Assert(x < dims[0]);
	Assert(y >= 0);
	Assert(y < dims[1]);
	Assert(z >= 0);
	Assert(z < dims[2]);

	return z * dims[0] * dims[1] + y * dims[0] + x;
}


#endif // __LATTICE_H__
