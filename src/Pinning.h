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
#ifndef __PINNING_H__
#define __PINNING_H__


int PinCurrentThreadToCore(int coreNumber);

int PinParseCpuList(const char * cpuList,
					int mpiRank, int nodeRank, int threadNumber);

int PinCurrentThreadByEnvVar(const char * envVarName,
	int mpiRank, int nodeRank, int threadNumber);

int PinCurrentThreadByCpuList(const char * cpuList,
	int mpiRank, int nodeRank, int threadNumber);

int PinCurrentCore();

char * PinCpuListAsString();


#endif // __PINNING_H__
