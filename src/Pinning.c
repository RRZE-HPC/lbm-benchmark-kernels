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
#ifndef _GNU_SOURCE
	#define _GNU_SOURCE
#endif
#include <sched.h>
#include <errno.h>


#include "Base.h"
#include "Pinning.h"

// -----------------------------------------------------------------------
//
// Binds the calling thread to specified core.
//
// Return value: 0 = success, else error.
//
// -----------------------------------------------------------------------

int PinCurrentThreadToCore(int coreNumber)
{
	int error = 0;

	cpu_set_t cpu_set;
	CPU_ZERO(&cpu_set);
	CPU_SET(coreNumber, &cpu_set);

	error = sched_setaffinity((pid_t)0, sizeof(cpu_set_t), &cpu_set);

	if (error != 0) {
		Error("pinning thread to core %d failed (%d): %s\n",
				coreNumber, error, strerror(error));
	}

	return error;
}


// -----------------------------------------------------------------------
//
// Binds the calling thread to specified core by a cpu list specified
// in the given environment variable.
//
// Return value: 0 = success, else error.
//
// -----------------------------------------------------------------------

int PinCurrentThreadByEnvVar(const char * envVarName, int threadNumber)
{
	const char * envVarValue;
	int core;

	envVarValue = getenv(envVarName);

	if (envVarValue == NULL) {
		printf("WARNING: skip pinning: env var %s not set\n", envVarName);

		return 0;
	}

	core = PinParseCpuList(envVarValue, threadNumber);

	if (core < 0) {
		return core;
	}

	return PinCurrentThreadToCore(core);
}


// -----------------------------------------------------------------------
//
// Binds the calling thread to a core specified in the CPU list.
//
// Return value: 0 = success, else error.
//
// -----------------------------------------------------------------------

int PinCurrentThreadByCpuList(const char * cpuList, int threadNumber)
{
	int core;

	if (cpuList == NULL) {
		printf("ERROR: cpu list is NULL.\n");

		exit(1);
	}

	core = PinParseCpuList(cpuList, threadNumber);

	if (core < 0) {
		return core;
	}

	return PinCurrentThreadToCore(core);
}


// -----------------------------------------------------------------------
//
// Parses the provided cpu list and returns the core number for the
// specified thread.
//
// The cpu list has for example a format of: 0,1,2
//
// -----------------------------------------------------------------------

int PinParseCpuList(const char * cpuList, int threadNumber)
{
	int cpu = -1;

	if (cpuList == NULL) {
		return -1;
	}

	const char * c = cpuList;

	// Ensure only valid characters are in the cpu list.
	// Cpu list is in the format of "0,1,2,3,4,5".
	while (((*c >= '0' && *c <= '9') || *c == ',')) {
		++c;
	}

	if (*c != 0x00) {
    	// Invalid character detected.
		return -2;
  	}

	c = cpuList;

	// Now find the core for the specified thread.

	int t = 0;

	while (t < threadNumber && *c != 0x00) {
		if (*c == ',') {
			++t;
		}
		++c;
	}

	if (t != threadNumber || *c < '0' || *c > '9') {
		// Cpu for this threadNumber not found.
		return -4;
	}

	cpu = atoi(c);

	return cpu;
}



// -----------------------------------------------------------------------
//
// Returns the first core from the calling thread's affinity set.
//
// On error a value < 0 is returned.
//
// -----------------------------------------------------------------------

int PinCurrentCore()
{
	int core = -1;
	int err;

	cpu_set_t cpu_set;
	CPU_ZERO(&cpu_set);

	err = sched_getaffinity((pid_t)0, sizeof(cpu_set_t), &cpu_set);

	// constant CPU_SETSIZE is one larger than the maximum CPU
	// number that can be stored in a CPU set
	for (int i = 0; i < CPU_SETSIZE; ++i) {
		if (CPU_ISSET(i, &cpu_set)) {
			core = i;
			break;
		}
	}

	if (err != 0) {
		Error("getting thread affinty failed (%d): %s\n", err, strerror(err));
		return -1;
	}

	return core;
}



// -----------------------------------------------------------------------
//
// Returns the all cores from the calling thread's affinity set.
//
// On error a value < 0 is returned.
//
// -----------------------------------------------------------------------

typedef cpu_set_t CpuSet;


static CpuSet PinCurrentCores()
{
	CpuSet cpuSet;
	int err;

	cpu_set_t cpu_set;
	CPU_ZERO(&cpu_set);

	err = sched_getaffinity((pid_t)0, sizeof(cpu_set_t), &cpu_set);

	cpuSet = cpu_set;

	if (err != 0) {
		Error("getting thread affinty failed (%d): %s\n", err, strerror(err));
		return cpuSet;
	}

	return cpuSet;
}

static char * CpuSetToString(cpu_set_t * cpu_set)
{
	int previousSetCore = -2;
	int rangeBeginCore = -2;

	char * buffer1 = (char *)malloc(1024);
	Assert(buffer1 != NULL);
	char * buffer2 = (char *)malloc(1024);
	Assert(buffer2 != NULL);

	buffer1[0] = 0x00;
	buffer2[0] = 0x00;

	char * buffer = buffer1;
	char * bufferOld = buffer2;

	const char * empty = "";
	const char * realComma = ",";
	const char * comma = empty;

	// TODO: use snprintf
	// TODO: increase allocated buffer if necessary

	for (int i = 0; i < CPU_SETSIZE; ++i) {
		if (!CPU_ISSET(i, cpu_set)) {
			continue;
		}

		if (i == previousSetCore + 1) {
			previousSetCore = i;
			continue;
		}

		// Now we reached the end of a range.
		// The range can also consist of only one core.
		// Be aware, that this core is not part of the range.

		// TODO: this code is repeated below -> use it only once
		if (rangeBeginCore >= 0 && previousSetCore >= 0) {
			char * tmp;

			tmp = buffer;
			buffer = bufferOld;
			bufferOld = tmp;

			if (rangeBeginCore < previousSetCore) {
				sprintf(buffer, "%s%s%d-%d", bufferOld, comma, rangeBeginCore, previousSetCore);
			}
			else {
				sprintf(buffer, "%s%s%d", bufferOld, comma, previousSetCore);
			}

			comma = realComma;
		}

		// With this core a new range begins.
		rangeBeginCore = i;
		previousSetCore = i;
	}

	if (rangeBeginCore >= 0 && previousSetCore >= 0) {
		char * tmp;

		tmp = buffer;
		buffer = bufferOld;
		bufferOld = tmp;

		if (rangeBeginCore < previousSetCore) {
			sprintf(buffer, "%s%s%d-%d", bufferOld, comma, rangeBeginCore, previousSetCore);
		}
		else {
			sprintf(buffer, "%s%s%d", bufferOld, comma, previousSetCore);
		}
	}

	free(bufferOld); bufferOld = NULL;

	return buffer;
}

char * PinCpuListAsString()
{
	CpuSet cpuSet = PinCurrentCores();

	return CpuSetToString(&cpuSet);
}

#ifdef TEST

int main(int argc, char * argv[])
{
	char * cpuList = PinCpuListAsString();

	printf("pinned to cores: %s\n", cpuList);

	free(cpuList); cpuList = NULL;

	return 0;
}

#endif // TEST

