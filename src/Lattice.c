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
#include "Lattice.h"

// Dumps the layers [zStart, zStop] of lattice as ASCII.
// Specify zStart = -1 and zStop = -1 as begin/end of lattice.

void LatDumpAscii(LatticeDesc * ld, int zStart, int zStop)
{
	Assert(ld != NULL);

	const char strLatCellType[] = "X.IxO"; // X = Obstacle, . = Fluid, I = inlet, O = outlet

	int localZStart = zStart;
	int localZStop  = zStop;

	int * dims = ld->Dims;
	LatticeT * lattice = ld->Lattice;

	if (localZStart == -1) localZStart = 0;
	if (localZStop  == -1) localZStop  = dims[2] - 1;

	for (int z = localZStop; z >= localZStart; --z) {
		printf("plane % 2d\n", z);

		for (int y = dims[1] - 1; y >= 0; --y) {
			printf(" %2d  ", y);
			for (int x = 0; x < dims[0]; ++x) {
				printf("%c", strLatCellType[lattice[L_INDEX_4(dims, x, y, z)]]);
			}
			printf("\n");
		}
	}
}

