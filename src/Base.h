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
#ifndef __BASE_H__
#define __BASE_H__

#include "Config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>


static inline double Time()
{
	struct timeval time;

	gettimeofday(&time, NULL);

	return (time.tv_sec + 1e-6 * time.tv_usec);
}


#define STRINGIFYX(x)   #x
#define STRINGIFY(x)    STRINGIFYX(x)

// See top of BoostJoin.h for Boost Licence.
#include "BoostJoin.h"
#define JOIN(X, Y)      BOOST_DO_JOIN(X, Y)


// Some macro fu to remove the first comma.
// "x" is an empty macro agrument in EXPAND2
// before the first comma which is skipped
#ifndef EXPAND
	#define EXPAND2(x, ...)       __VA_ARGS__
	#define EXPAND(x, ...)        EXPAND2(x, ## __VA_ARGS__)
#endif

#ifdef DEBUG

    #define Assert(expression) \
        do { \
            if (!(expression)) { \
                Error("%s:%d assertion \"%s\" failed with code %d\n", \
                        __FILE__, __LINE__, \
                        #expression, expression); \
				 __asm__ ("int $3\n"); \
                exit(-1); \
            } \
        } while (0)

    #define AssertMsg(expression, formatString, ...) \
        do { \
            if (!(expression)) { \
                Error("%s:%d assertion \"%s\" failed with code %d\n", \
                        __FILE__, __LINE__, \
                        #expression, expression); \
                Error(formatString, ##__VA_ARGS__); \
				 __asm__ ("int $3\n"); \
                exit(-1); \
            } \
        } while (0)
#else

	#define Assert(expression)
	#define AssertMsg(expression, formatString, ...)

#endif

    #define Verify(expression) \
        do { \
            if (!(expression)) { \
                Error("%s:%d verify \"%s\" failed with code %d\n", \
                        __FILE__, __LINE__, \
                        #expression, expression); \
				 __asm__ ("int $3\n"); \
                exit(-1); \
            } \
        } while (0)

    #define VerifyMsg(expression, formatString, ...) \
        do { \
            if (!(expression)) { \
                Error("%s:%d verify \"%s\" failed with code %d\n", \
                        __FILE__, __LINE__, \
                        #expression, expression); \
                Error(formatString, ##__VA_ARGS__); \
				 __asm__ ("int $3\n"); \
                exit(-1); \
            } \
        } while (0)

    #define Print(formatString, ...) \
        fprintf(stdout, SHC_MAGENTA formatString SHC_NC, ##__VA_ARGS__)

    #define Warning(formatString, ...) \
        fprintf(stdout, SHC_BROWN "WARNING: " SHC_NC formatString, ##__VA_ARGS__)

    #define Error(formatString, ...) \
        fprintf(stderr, SHC_RED "ERROR: " formatString SHC_NC , ##__VA_ARGS__)
/*
    #define DebugPrint(formatString, ...) \
        fprintf(stderr, "DEBUG: " formatString, ##__VA_ARGS__)
*/
    #ifndef NO_SHELL_COLORS

        // or "\e"
        #define ESC             "\x1b"

        // No Color
        #define SHC_NC          ESC "[0m"

        #define SHC_BLACK       ESC "[0;30m"
        #define SHC_MAGENTA     ESC "[0;35m"
        #define SHC_RED         ESC "[0;31m"
        #define SHC_DARK_RED    ESC "[1;31m"
        #define SHC_CYAN        ESC "[0;36m"
        #define SHC_BROWN       ESC "[0;33m"
        #define SHC_DARK_GREEN  ESC "[1;32m"

    #else  // NO_SHELL_COLORS

        // No Color
        #define SHC_NC          ""

        #define SHC_BLACK       ""
        #define SHC_MAGENTA     ""
        #define SHC_RED         ""
        #define SHC_DARK_RED    ""
        #define SHC_CYAN        ""
        #define SHC_BROWN       ""
        #define SHC_DARK_GREEN  ""

    #endif  // NO_SHELL_COLORS


    #define N_ELEMS(x)  (sizeof(x) / sizeof((x)[0]))


	#define MIN(a, b)		((a) <= (b) ? (a) : (b))

static inline int MinI(int a, int b) { return a <= b ? a : b; }

// Raises a breakpoint if a debugger is attached, else SIG_TRAP is raised.
#define DEBUG_BREAK_POINT()		__asm__ ("int $3\n")

#define UNUSED(variable)                    (void)(variable)


static inline char * ByteToHuman(size_t bytes)
{
	static char buffer[256] = { 0 };

	if (bytes < 1024) {
		snprintf(buffer, sizeof(buffer), "%lu b", bytes);
		return buffer;
	}

	double KiB = bytes / 1024.0;

	if (KiB < 1024.0) {
		snprintf(buffer, sizeof(buffer), "%9.2e KiB", KiB);
		return buffer;
	}

	double MiB = KiB / 1024.0;
	if (MiB < 1024.0) {
		snprintf(buffer, sizeof(buffer), "%9.2e MiB", MiB);
		return buffer;
	}

	double GiB = MiB / 1024.0;
	snprintf(buffer, sizeof(buffer), "%9.2e GiB", GiB);
	return buffer;
}


#endif // __BASE_H__
