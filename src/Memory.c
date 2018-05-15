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
//   Michael Hussnaetter, 2017-2018
//   University of Erlangen-Nuremberg, Germany
//   michael.hussnaetter -at- fau.de
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
// TODO: make configurable
#define HAVE_HUGE_PAGES


#ifdef HAVE_HUGE_PAGES
#define _BSD_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  // strerror
#include <errno.h>


#ifdef HAVE_MEMKIND
#include <hbwmalloc.h>
#endif

#ifdef HAVE_HUGE_PAGES
#include <sys/mman.h> // madvise
#endif

#include "Base.h"
#include "Memory.h"


int MemAlloc(void ** ptr, size_t bytesToAlloc)
{
	void * tmpPtr;

	tmpPtr = malloc(bytesToAlloc);

	if (tmpPtr == NULL) { //  && bytesToAlloc != 0) {
		Error("allocation of %lu bytes failed: %d - %s\n", bytesToAlloc, errno, strerror(errno));
		exit(1);
	}

	*ptr = tmpPtr;

	return 0;
}

int MemAllocAligned(void ** ptr, size_t bytesToAlloc, size_t alignmentBytes)
{
	int ret;

	ret = posix_memalign(ptr, alignmentBytes, bytesToAlloc);

	if (ret) {
		Error("allocation of %lu bytes aligned to %lu bytes failed: %d - %s\n", bytesToAlloc, alignmentBytes, errno, strerror(errno));
		exit(1);
	}

#ifdef HAVE_HUGE_PAGES

	if (alignmentBytes % 4096 == 0) {
		ret = madvise(*ptr, bytesToAlloc, MADV_HUGEPAGE);

		if (ret != 0) {
			DEBUG_BREAK_POINT();
			Error("madvise(%p, %lu, MADV_HUGEPAGE) failed: %d - %s.\n", *ptr, bytesToAlloc, errno, strerror(errno));
			exit(1);
		}
	}
#endif

	return 0;
}


int MemFree(void ** ptr)
{
	Assert(*ptr != NULL);

	free(*ptr);

	*ptr = NULL;

	return 0;
}

int MemZero(void * ptr, size_t bytesToZero)
{
	Assert(ptr != NULL);
	Assert(bytesToZero > 0);

	memset(ptr, 0, bytesToZero);

	return 0;
}

#ifdef HAVE_MEMKIND
int HbwAlloc(void ** ptr, size_t bytesToAlloc)
{
	void * tmpPtr;

	tmpPtr = hbw_malloc(bytesToAlloc);

	if (tmpPtr == NULL) { //  && bytesToAlloc != 0) {
		Error("allocation of %lu bytes in HBM failed: %d - %s\n", bytesToAlloc, errno, strerror(errno));
		exit(1);
	}

	*ptr = tmpPtr;

	return 0;
}

int HbwAllocAligned(void ** ptr, size_t bytesToAlloc, size_t alignmentBytes)
{
	int ret;

	ret = hbw_posix_memalign(ptr, alignmentBytes, bytesToAlloc);

	if (ret) {
		Error("allocation of %lu bytes in HBM aligned to %lu bytes failed: %d - %s\n", bytesToAlloc, alignmentBytes, errno, strerror(errno));
		exit(1);
	}

	return 0;
}


int HbwFree(void ** ptr)
{
	Assert(*ptr != NULL);

	hbw_free(*ptr);

	*ptr = NULL;

	return 0;
}

#endif	// HAVE_MEMKIND
