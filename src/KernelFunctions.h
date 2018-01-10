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
#ifndef __KERNEL_FUNCTIONS_H__
#define __KERNEL_FUNCTIONS_H__

#include "BenchKernelD3Q19.h"
#include "BenchKernelD3Q19Aa.h"
#include "BenchKernelD3Q19AaVec.h"
#include "BenchKernelD3Q19AaVecSl.h"
#include "BenchKernelD3Q19List.h"
#include "BenchKernelD3Q19ListAa.h"
#include "BenchKernelD3Q19ListAaRia.h"
#include "BenchKernelD3Q19ListAaPv.h"
#include "BenchKernelD3Q19ListPullSplitNt.h"

typedef struct KernelFunctions_
{
	char Name[128];
	void (* Init)(LatticeDesc * ld, KernelData ** kernelData, Parameters * params);
	void (* Deinit)(LatticeDesc * ld, KernelData ** kernelData);
} KernelFunctions;

KernelFunctions g_kernels[] =
{
	{
		.Name   = "list-aa-pv-soa",
		.Init   = D3Q19ListAaPvInit_PushSoA,
		.Deinit = D3Q19ListAaPvDeinit_PushSoA
	},
	{
		.Name   = "list-aa-ria-soa",
		.Init   = D3Q19ListAaRiaInit_PushSoA,
		.Deinit = D3Q19ListAaRiaDeinit_PushSoA
	},
	{
		.Name   = "list-aa-soa",
		.Init   = D3Q19ListAaInit_PushSoA,
		.Deinit = D3Q19ListAaDeinit_PushSoA
	},
	{
		.Name   = "list-aa-aos",
		.Init   = D3Q19ListAaInit_PushAoS,
		.Deinit = D3Q19ListAaDeinit_PushAoS
	},
	{
		.Name   = "list-pull-split-nt-1s-soa",
		.Init   = D3Q19ListPullSplitNt1SInit_PullSoA,
		.Deinit = D3Q19ListPullSplitNtDeinit_PullSoA
	},
	{
		.Name   = "list-pull-split-nt-2s-soa",
		.Init   = D3Q19ListPullSplitNt2SInit_PullSoA,
		.Deinit = D3Q19ListPullSplitNtDeinit_PullSoA
	},
	{
		.Name   = "list-push-soa",
		.Init   = D3Q19ListInit_PushSoA,
		.Deinit = D3Q19ListDeinit_PushSoA
	},
	{
		.Name   = "list-push-aos",
		.Init   = D3Q19ListInit_PushAoS,
		.Deinit = D3Q19ListDeinit_PushAoS
	},
	{
		.Name   = "list-pull-soa",
		.Init   = D3Q19ListInit_PullSoA,
		.Deinit = D3Q19ListDeinit_PullSoA
	},
	{
		.Name = "list-pull-aos",
		.Init = D3Q19ListInit_PullAoS,
		.Deinit = D3Q19ListDeinit_PullAoS
	},
	{
		.Name   = "push-soa",
		.Init   = D3Q19Init_PushSoA,
		.Deinit = D3Q19Deinit_PushSoA
	},
	{
		.Name   = "push-aos",
		.Init   = D3Q19Init_PushAoS,
		.Deinit = D3Q19Deinit_PushAoS
	},
	{
		.Name   = "pull-soa",
		.Init   = D3Q19Init_PullSoA,
		.Deinit = D3Q19Deinit_PullSoA
	},
	{
		.Name   = "pull-aos",
		.Init   = D3Q19Init_PullAoS,
		.Deinit = D3Q19Deinit_PullAoS
	},
	{
		.Name   = "blk-push-soa",
		.Init   = D3Q19BlkInit_PushSoA,
		.Deinit = D3Q19BlkDeinit_PushSoA
	},
	{
		.Name   = "blk-push-aos",
		.Init   = D3Q19BlkInit_PushAoS,
		.Deinit = D3Q19BlkDeinit_PushAoS
	},
	{
		.Name   = "blk-pull-soa",
		.Init   = D3Q19BlkInit_PullSoA,
		.Deinit = D3Q19BlkDeinit_PullSoA
	},
	{
		.Name   = "blk-pull-aos",
		.Init   = D3Q19BlkInit_PullAoS,
		.Deinit = D3Q19BlkDeinit_PullAoS
	},
	{
		.Name	= "aa-aos",
		.Init	= D3Q19AaInit_AaAoS,
		.Deinit	= D3Q19AaDeinit_AaAoS
	},
	{
		.Name	= "aa-soa",
		.Init	= D3Q19AaInit_AaSoA,
		.Deinit	= D3Q19AaDeinit_AaSoA
	},
	{
		.Name	= "aa-vec-soa",
		.Init	= D3Q19AaVecInit_AaSoA,
		.Deinit	= D3Q19AaVecDeinit_AaSoA
	},
	{
		.Name	= "aa-vec-sl-soa",
		.Init	= D3Q19AaVecSlInit_AaSoA,
		.Deinit	= D3Q19AaVecSlDeinit_AaSoA
	}


};

#endif // __KERNEL_FUNCTIONS_H__
