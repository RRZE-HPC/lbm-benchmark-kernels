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
#ifndef __BENCH_KERNEL_D3Q19_LIST__
#define __BENCH_KERNEL_D3Q19_LIST__

#include "Kernel.h"

void D3Q19ListInit_PushSoA(LatticeDesc * ld, KernelData ** kernelData, Parameters * params);
void D3Q19ListInit_PushAoS(LatticeDesc * ld, KernelData ** kernelData, Parameters * params);

void D3Q19ListDeinit_PushSoA(LatticeDesc * ld, KernelData ** kernelData);
void D3Q19ListDeinit_PushAoS(LatticeDesc * ld, KernelData ** kernelData);

void D3Q19ListInit_PullSoA(LatticeDesc * ld, KernelData ** kernelData, Parameters * params);
void D3Q19ListInit_PullAoS(LatticeDesc * ld, KernelData ** kernelData, Parameters * params);

void D3Q19ListDeinit_PullSoA(LatticeDesc * ld, KernelData ** kernelData);
void D3Q19ListDeinit_PullAoS(LatticeDesc * ld, KernelData ** kernelData);

#endif // __BENCH_KERNEL_D3Q19_LIST__
