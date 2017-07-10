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
#ifndef __LIKWID_IF_H__
#define __LIKWID_IF_H__

#ifdef HAVE_LIKWID
	#include <likwid.h>


	#define X_LIKWID_INIT()						likwid_markerInit()
	#define X_LIKWID_DEINIT()					likwid_markerClose()
    #define X_LIKWID_START(region_name)			likwid_markerStartRegion(region_name)
    #define X_LIKWID_STOP(region_name)			likwid_markerStopRegion(region_name)
#else

	#define X_LIKWID_INIT()
	#define X_LIKWID_DEINIT()
    #define X_LIKWID_START(region_name)
    #define X_LIKWID_STOP(region_name)
#endif

#endif // __LIKWID_IF_H__
