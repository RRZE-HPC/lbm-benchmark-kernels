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
#include "Geometry.h"
#include "Memory.h"

#include <strings.h>
#include <math.h>

#include <errno.h>

void GeoCreateByStr(const char * geometryType, int dims[3], int periodic[3], LatticeDesc * ld)
{
	int type = -1;
	void * typeDetails = NULL;
	int tmp;

	if (strncasecmp("channel", geometryType, 7) == 0) {
		type = GEO_TYPE_CHANNEL;
	}
	else if (strncasecmp("box", geometryType, 3) == 0) {
		type = GEO_TYPE_BOX;
	}
	else if (strncasecmp("pipe", geometryType, 4) == 0) {
		type = GEO_TYPE_PIPE;
	}
	else if (strncasecmp("blocks", geometryType, 6) == 0) {
		type = GEO_TYPE_BLOCKS;

		// Default block size
		tmp = 8;

		if (strlen(geometryType) > 7) {
			int blockSize = atoi(&geometryType[7]);

			int dimMin = dims[0];

			if (dims[1] < dimMin) dimMin = dims[1];
			if (dims[2] < dimMin) dimMin = dims[2];

			if (blockSize < 0 || blockSize > dimMin / 2) {
				printf("ERROR: block size for geometry must be > 0 and smaller than half of the smalest dimension.\n");
				// TODO: find a better solution for handling errors in here.
				Verify(0);
			}

			tmp = blockSize;
		}

		typeDetails = &tmp;
	}
	else {
		printf("ERROR: unknown geometry specified.\n");
		Verify(0);
	}

	GeoCreateByType(type, typeDetails, dims, periodic, ld);

	return;
}

void GeoCreateByType(GEO_TYPES type, void * typeDetails, int dims[3], int periodic[3], LatticeDesc * ld)
{
	Assert(dims != NULL);
	Assert(dims[0] > 0);
	Assert(dims[1] > 0);
	Assert(dims[2] > 0);

	Assert(periodic != NULL);
	Assert(periodic[0] >= 0);
	Assert(periodic[1] >= 0);
	Assert(periodic[2] >= 0);

	Assert(ld != NULL);

	Assert(type >= GEO_TYPE_MIN);
	Assert(type <= GEO_TYPE_MAX);

	const char * geoTypeStr[] = { "box", "channel", "pipe", "blocks" };

	printf("# geometry: %d x %d x %d nodes, type %d %s\n", dims[0], dims[1], dims[2], type, geoTypeStr[type]);

	ld->Dims[0] = dims[0];
	ld->Dims[1] = dims[1];
	ld->Dims[2] = dims[2];
	ld->nCells = dims[0] * dims[1] * dims[2];
	ld->PeriodicX = periodic[0];
	ld->PeriodicY = periodic[1];
	ld->PeriodicZ = periodic[2];

	LatticeT * lattice;
	MemAlloc((void **)&lattice, sizeof(LatticeT) * dims[0] * dims[1] * dims[2]);

	ld->Lattice = lattice;

	for (int z = 0; z < dims[2]; ++z) {
		for (int y = 0; y < dims[1]; ++y) {
			for (int x = 0; x < dims[0]; ++x) {
				lattice[L_INDEX_4(dims, x, y, z)] = LAT_CELL_FLUID;
			}
		}
	}

	if (type == GEO_TYPE_CHANNEL || type == GEO_TYPE_BLOCKS || type == GEO_TYPE_PIPE) {
		periodic[0] = 1;
	}

	// Walls or periodic on first and last x plane.
	for (int z = 0; z < dims[2]; ++z) {
		for (int y = 0; y < dims[1]; ++y) {
			if(periodic[0]){
				lattice[L_INDEX_4(dims, 0, y, z)] 				= LAT_CELL_FLUID;
				lattice[L_INDEX_4(dims, dims[0] - 1, y, z)] 	= LAT_CELL_FLUID;
			} else {
				lattice[L_INDEX_4(dims, 0, y, z)] 				= LAT_CELL_OBSTACLE;
				lattice[L_INDEX_4(dims, dims[0] - 1, y, z)] 	= LAT_CELL_OBSTACLE;
			}
		}
	}

	// Walls or periodic on first and last y plane.
	for (int z = 0; z < dims[2]; ++z) {
		for (int x = 0; x < dims[0]; ++x) {
			if(periodic[1]){
				lattice[L_INDEX_4(dims, x, 0, z)] 				= LAT_CELL_FLUID;
				lattice[L_INDEX_4(dims, x, dims[1] - 1, z)] 	= LAT_CELL_FLUID;
			} else {
				lattice[L_INDEX_4(dims, x, 0, z)] 				= LAT_CELL_OBSTACLE;
				lattice[L_INDEX_4(dims, x, dims[1] - 1, z)] 	= LAT_CELL_OBSTACLE;
			}
		}
	}

	// Walls or periodic on first and last z plane.
	for (int y = 0; y < dims[1]; ++y) {
		for (int x = 0; x < dims[0]; ++x) {
			if(periodic[2]){
				lattice[L_INDEX_4(dims, x, y, 0)] 				= LAT_CELL_FLUID;
				lattice[L_INDEX_4(dims, x, y, dims[2] - 1)] 	= LAT_CELL_FLUID;
			} else {
				lattice[L_INDEX_4(dims, x, y, 0)] 				= LAT_CELL_OBSTACLE;
				lattice[L_INDEX_4(dims, x, y, dims[2] - 1)] 	= LAT_CELL_OBSTACLE;
			}
		}
	}

	if (type == GEO_TYPE_CHANNEL) {
		periodic[0] = 1;
	}
	else if (type == GEO_TYPE_PIPE) {
		#define SQR(a) ((a)*(a))
		double centerZ = dims[2] / 2.0 - 0.5;
		double centerY = dims[1] / 2.0 - 0.5;
		double minDiameter = MIN(dims[1], dims[2]);
		double minRadiusSquared = SQR(minDiameter / 2 - 1);

		for (int z = 0; z < dims[2]; ++z) {
			for (int y = 0; y < dims[1]; ++y) {
				if((SQR(z - centerZ) + SQR(y - centerY)) >= minRadiusSquared) {
					for (int x = 0; x < dims[0]; ++x) {
						lattice[L_INDEX_4(dims, x, y, z)] 	= LAT_CELL_OBSTACLE;
					}
				}
			}
		}
		#undef SQR
	}
	else if (type == GEO_TYPE_BLOCKS) {

		int blockSize = *((int *)typeDetails);

		if (blockSize == 0) {
			blockSize = 8;
		}

		int dimMin = dims[0];

		if (dims[1] < dimMin) dimMin = dims[1];
		if (dims[2] < dimMin) dimMin = dims[2];

		if (blockSize < 0 || blockSize > dimMin / 2) {
			printf("ERROR: block size for geometry must be > 0 and smaller than half of the smalest dimension.\n");
			// TODO: find a better solution for handling errors in here.
			Verify(0);
		}

		// Number of blocks in x, y, and z direction.
		int nbx = blockSize, nby = blockSize, nbz = blockSize;

		for (int z = 0; z < dims[2]; ++z) {
				if ((z % (2 * nbz)) < nbz) continue;

			for (int y = 0; y < dims[1]; ++y) {
				if ((y % (2 * nby)) < nby) continue;

				for (int x = 0; x < dims[0]; ++x) {

					if ((x % (2 * nbx)) >= nbx) {
						lattice[L_INDEX_4(dims, x, y, z)] 	= LAT_CELL_OBSTACLE;
					}
				}
			}
		}
	}

// 	if (latticeDumpAscii) {
// 		const char strLatCellType[] = "X.IxO"; // X = Obstacle, . = Fluid, I = inlet, O = outlet
// 		for (int z = dims[2] - 1; z >= 0; --z) {
// 			printf("plane % 2d\n", z);
//
// 			for (int y = dims[1] - 1; y >= 0; --y) {
// 				printf(" %2d  ", y);
// 				for (int x = 0; x < dims[0]; ++x) {
// 					printf("%c", strLatCellType[lattice[L_INDEX_4(dims, x, y, z)]]);
// 				}
// 				printf("\n");
// 			}
// 		}
//	}

// Lattice Helper Function

	ld->nObst = 0;
	ld->nFluid = 0;
	ld->nInlet = 0;
	ld->nOutlet = 0;

	for (int z = 0; z < dims[2]; ++z) {
		for (int y = 0; y < dims[1]; ++y) {
			for (int x = 0; x < dims[0]; ++x) {
				switch (lattice[L_INDEX_4(dims, x, y, z)]) {
					case LAT_CELL_OBSTACLE:					ld->nObst++; break;
					case LAT_CELL_FLUID:					ld->nFluid++; break;
					case LAT_CELL_INLET:					ld->nInlet++;  ld->nFluid++; break;
					case LAT_CELL_OUTLET:					ld->nOutlet++; ld->nFluid++; break;
					default:
						Verify(0);
				}
			}
		}
	}

	return;
}
