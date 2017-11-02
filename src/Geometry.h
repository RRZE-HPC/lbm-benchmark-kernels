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
#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include "Lattice.h"


typedef enum GEO_TYPES_ {
	GEO_TYPE_MIN = 0,
	GEO_TYPE_BOX = 0,
	GEO_TYPE_CHANNEL = 1,
	GEO_TYPE_PIPE = 2,
	GEO_TYPE_BLOCKS = 3,	// Expects a pointer to an integer, holding the
						    // value of the block size as type detail.
	GEO_TYPE_FLUID = 4,
	GEO_TYPE_MAX = 4
} GEO_TYPES;


void GeoCreateByType(GEO_TYPES type, void * typeDetails, int dims[3], int periodic[3], LatticeDesc * ld);

void GeoCreateByStr(const char * geometryType, int dims[3], int periodic[3], LatticeDesc * ld);

#endif // __GEOMETRY_H__
