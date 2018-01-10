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
#include "BenchKernelD3Q19AaVecSlCommon.h"
#include "BenchKernelD3Q19AaVec.h"


#include "Memory.h"
#include "Vtk.h"
#include "Vector.h"

#include <inttypes.h>
#include <math.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

// Forward definition.
void FNAME(D3Q19AaVecSlKernel)(LatticeDesc * ld, struct KernelData_ * kd, CaseData * cd);

void FNAME(D3Q19AaVecSlInit)(LatticeDesc * ld, KernelData ** kd, Parameters * params)
{
	FNAME(D3Q19AaVecInit)(ld, kd, params);

	(*kd)->Kernel = FNAME(D3Q19AaVecSlKernel);

	return;
}

void FNAME(D3Q19AaVecSlDeinit)(LatticeDesc * ld, KernelData ** kd)
{
	FNAME(D3Q19AaVecDeinit)(ld, kd);

	return;
}

