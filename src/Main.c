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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>  // strcasecmp

#include <inttypes.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Base.h"
#include "Kernel.h"
#include "Memory.h"

#include "Lattice.h"
#include "Geometry.h"
#include "Pinning.h"
#include "LikwidIf.h"

#include "KernelFunctions.h"

#ifdef __x86_64__
	#include <xmmintrin.h>


	#define MXCSR_DAZ					6
	#define MXCSR_FTZ 					15


	int FpIsMxCsrMaskSet(unsigned int mask)
	{
		unsigned int mxcsr;
		unsigned int mxcsrNew;

		mxcsr = _mm_getcsr();

		mxcsrNew = mxcsr & mask;

		return (mxcsrNew == mask);
	}

	int FpGetFtz()
	{
		return FpIsMxCsrMaskSet(1 << MXCSR_FTZ);
	}

	int FpGetDaz()
	{
		return FpIsMxCsrMaskSet(1 << MXCSR_DAZ);
	}
#endif


int ParseDimensions(const char * parameter, int * nX, int * nY, int * nZ)
{
	char * tmp;

	*nX = atoi(parameter);

	if (*nX <= 0) {
		printf("ERROR: parameter for X dimension must be > 0.\n");
		return 0;
	}

	tmp = strchr(parameter, 'x');

	if (tmp == NULL) {
		printf("ERROR: parameter for Y dimension is missing.\n");
		return 0;
	}

	*nY = atoi(tmp + 1);

	if (*nY <= 0) {
		printf("ERROR: parameter for Y dimension must be > 0.\n");
		return 0;
	}

	tmp = strchr(tmp + 1, 'x');

	if (tmp == NULL) {
		printf("ERROR: parameter for Z dimension is missing.\n");
		return 0;
	}

	*nZ = atoi(tmp + 1);

	if (*nZ <= 0) {
		printf("ERROR: parameter for Z dimension must be > 0.\n");
		return 0;
	}

	return 1;
}

int main(int argc, char * argv[])
{
	int dims[3] = { 20, 20, 20 };		// Dimensions in x, y, and z direction
	const char * geometryType = "channel";
	// int latticeDumpAscii = 0;
	int verify = 0; UNUSED(verify);
	char * kernelToUse = NULL;
	int nThreads = 1;
	const char * pinString = NULL;
	int periodic[3] = { 0 };

	CaseData cd;

	cd.MaxIterations 		= 10;
	cd.RhoIn		 		= F(1.0);
	cd.RhoOut 		 		= F(1.0);
	cd.Omega		 		= F(1.0);
	cd.VtkOutput	 		= 0;
	cd.VtkModulus    		= 100;
	cd.StatisticsModulus 	= 100;
	cd.XForce				= F(0.00001);
	kernelToUse				= "push-soa";

	Parameters p;
	p.nArgs        = argc;
	p.Args         = argv;
	p.nKernelArgs  = 0;
	p.KernelArgs   = NULL;

#define LBM_BENCH_KERNELS_VERSION_MAJOR 0
#define LBM_BENCH_KERNELS_VERSION_MINOR 1

    printf("Lattice Boltzmann Benchmark Kernels (LbmBenchKernels) Copyright (C) 2016, 2017 LSS, RRZE\n");
    printf("This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.\n");
	printf("This is free software, and you are welcome to redistribute it under certain conditions.\n");
	printf("\n");
	printf("# LBM Benchmark Kernels %d.%d, compiled %s %s, type: %s\n",
		LBM_BENCH_KERNELS_VERSION_MAJOR, LBM_BENCH_KERNELS_VERSION_MINOR, __DATE__, __TIME__,
#ifdef VERIFICATION
		"verification"
#else
		"benchmark"
#endif
	);

	// ----------------------------------------------------------------------
	// Parse command line arguments

	#define ARG_IS(param)                   (!strcmp(argv[i], param))
	#define NEXT_ARG_PRESENT() \
		do { \
			if (i + 1 >= argc) { \
				printf("ERROR: argument %s requires a parameter.\n", argv[i]); \
				return 1; \
			} \
		} while (0)

    for (int i = 1; i < argc; ++i) {

        if (ARG_IS("-dims") || ARG_IS("--dims")) {
            NEXT_ARG_PRESENT();


            if (!ParseDimensions(argv[++i], &dims[0], &dims[1], &dims[2])) {
                return 1;
            }
        }
		// else if (ARG_IS("-lattice-dump-ascii") || ARG_IS("--lattice-dump-ascii")) {
		// 	latticeDumpAscii = 1;
		// }
		else if (ARG_IS("-geometry") || ARG_IS("--geometry")) {
			NEXT_ARG_PRESENT();

			geometryType = argv[++i];
		}
		else if (ARG_IS("-iterations") ||ARG_IS("--iterations")) {
			NEXT_ARG_PRESENT();

			cd.MaxIterations = strtol(argv[++i], NULL, 0);

			if (cd.MaxIterations <= 0) {
				printf("ERROR: number of iterations must be > 0.\n");
				return 1;
			}
		}
		else if (ARG_IS("-rho-in") ||ARG_IS("--rho-in")) {
			NEXT_ARG_PRESENT();

			cd.RhoIn = F(strtod(argv[++i], NULL));
		}
		else if (ARG_IS("-rho-out") ||ARG_IS("--rho-out")) {
			NEXT_ARG_PRESENT();

			cd.RhoOut = F(strtod(argv[++i], NULL));
		}
		else if (ARG_IS("-omega") ||ARG_IS("--omega")) {
			NEXT_ARG_PRESENT();

			cd.Omega = F(strtod(argv[++i], NULL));
		}
		else if (ARG_IS("-x-force") ||ARG_IS("--x-force")) {
			NEXT_ARG_PRESENT();

			cd.XForce = F(strtod(argv[++i], NULL));
		}
		else if (ARG_IS("-verify") || ARG_IS("--verify")) {
#ifdef VERIFICATION

			// Choose this preset for verification. As geometry type "box" is
			// used but x and y direction are made periodic.
			// Everything else can be altered, but enough iterations should be
			// performed in order to receive a fully developed flow field.
			verify = 1;

			cd.Omega  = F(1.0);
			cd.RhoIn  = F(1.0);
			cd.RhoOut = F(1.0);
			geometryType = "box";
			dims[0] = 16;
			dims[1] = 16;
			dims[2] = 16;
			cd.XForce = F(0.00001);
			cd.MaxIterations = 1000;
			periodic[0] = 1;
			periodic[1] = 1;
			periodic[2] = 0;

			printf("#\n");
			printf("# VERIFICATION: verifying flow profile of channel flow.\n");
			printf("#\n");

			// TODO: this is not a good idea as we ignore all other options...

#else
			printf("ERROR: in order to use -verify VERIFICATION must be defined during compilation.\n");
			printf("       Recompile with VERIFICATION=on.\n");
			return 1;
#endif
		}
		else if (ARG_IS("-vtk") || ARG_IS("--vtk")) {
#ifdef VTK_OUTPUT

			cd.VtkOutput = 1;

			// If the next parameter is a number it is used as the itartion count,
			// if not it is probably another parameter.
			if (i + 1 < argc) {

				int vtkModulus = atoi(argv[i+1]);

				if (vtkModulus > 0) {
					cd.VtkModulus = vtkModulus;
					++i;
				}
			}
#else
			printf("ERROR: in order to use -vtk VTK_OUTPUT must be defined during compilation.\n");
			printf("       Recompile with VTK_OUTPUT=on.\n");
			return 1;
#endif
		}
		else if (ARG_IS("-statistics") || ARG_IS("--statistics")) {
#ifdef STATISTICS
			NEXT_ARG_PRESENT();

			cd.StatisticsModulus = atoi(argv[++i]);

			if (cd.StatisticsModulus <= 0) {
				printf("ERROR: the iteration count for -statistics must be > 0.\n");
				return 1;
			}
#else
			printf("ERROR: in order to use -statistics STATISTICS must be defined during compilation.\n");
			printf("       Recompile with STATISTICS=on.\n");
			return 1;
#endif
		}
		else if (ARG_IS("-kernel") || ARG_IS("--kernel")) {
			NEXT_ARG_PRESENT();

			kernelToUse = argv[++i];
		}
		else if (ARG_IS("-list") || ARG_IS("--list")) {
			printf("Available kernels to benchmark:\n");

			for (int j = 0; j < N_ELEMS(g_kernels); ++j) {
				printf("   %s\n", g_kernels[j].Name);
			}

			return 0;
		}
		else if (ARG_IS("-pin") || ARG_IS("--pin")) {
			NEXT_ARG_PRESENT();

			pinString = argv[++i];
		}
		else if (ARG_IS("-t") || ARG_IS("-threads") || ARG_IS("--threads")) {
#ifdef _OPENMP
			NEXT_ARG_PRESENT();

			nThreads = atoi(argv[++i]);

			if (nThreads <= 0) {
				printf("ERROR: number of threads must be > 0.\n");
				return 1;
			}
#else
			printf("ERROR: specifying number of threads is only available when compiled with OpenMP support.\n");
			return 1;
#endif
		}
		else if (ARG_IS("-periodic-x") || ARG_IS("--periodic-x")) {
			periodic[0] = 1;
		}
		else if (ARG_IS("-periodic-y") || ARG_IS("--periodic-y")) {
			periodic[1] = 1;
		}
		else if (ARG_IS("-periodic-z") || ARG_IS("--periodic-z")) {
			periodic[2] = 1;
		}
		else if (ARG_IS("-h") || ARG_IS("-help") || ARG_IS("--help")) {
			printf("ERROR: unknown argument: %s\n", argv[i]);
			printf("\n");
			printf("Usage:\n");
			printf("./lbmbenchk -list\n");
			printf("./lbmbenchk \n");
			printf("      [-dims XxYyZ] [-geometry box|channel|pipe|porosity[-value]] [-iterations <iterations>] [-lattice-dump-ascii]\n");
			printf("      [-rho-in <density>] [-rho-out <density] [-omega <omega>] [-kernel <kernel>]\n");
			printf("      [-periodic-x]\n");
#ifdef STATISTICS
			printf("      [-statistics <every-n-iteration>]\n");
#endif
#ifdef VTK_OUTPUT
			printf("      [-vtk [<every-n-iteration>]]\n");
#endif
#ifdef _OPENMP
			printf("      [-t <number of threads>]\n");
#endif
			printf("      [-pin core{,core}*]\n");
#ifdef VERIFICATION
			printf("      [-verify]\n");
#endif
			printf("      -- <kernel specific parameters>\n");
			printf("\n");
			printf("-list           List available kernels.\n");
			printf("\n");
			printf("-dims XxYxZ		Specify geometry dimensions.\n");
			printf("\n");
			printf("-geometry porosity-<value>\n");
			printf("                Geometetry with blocks of size <value> regularily layout out.\n");
			printf("\n");
			return 1;
		}
		else if (ARG_IS("--")) {
			// printf("# kernel args start with %s  these are %d args.\n", argv[i + 1], argc - i - 1);
			p.KernelArgs = &argv[++i];
			p.nKernelArgs = argc - i;
			break;
		}
		else {
			printf("ERROR: unknown parameter: %s.\n", argv[i]);
			exit(1);
		}
	}

	#undef ARG_IS
	#undef NEXT_ARG_PRESENT


	// ----------------------------------------------------------------------
	// Check if we exceed our index addressing PDFs.

	{
		uint64_t nPdfs = ((uint64_t)19) * dims[0] * dims[1] * dims[2];

		if (nPdfs > ((2LU << 31) - 1)) {
			printf("ERROR: number of PDFs exceed 2^31.\n");
			exit(1);
		}
	}

	// ----------------------------------------------------------------------

#ifdef _OPENMP
	omp_set_num_threads(nThreads);
#endif

	const char * defines[] = {
#ifdef DEBUG
	"DEBUG",
#endif
#ifdef VTK_OUTPUT
	"VTK_OUTPUT",
#endif
#ifdef STATISTICS
	"STATISTICS",
#endif
#ifdef VERIFICATION
	"VERIFICATION",
#endif
#ifdef _OPENMP
	"_OPENMP",
#endif
#ifdef HAVE_LIKWID
	"HAVE_LIKWID",
#endif
#ifdef INTEL_OPT_DIRECTIVES
	"INTEL_OPT_DIRECTIVES",
#endif
	};

	printf("#\n");

#ifdef PRECISION_DP
	printf("# - floating point:    double precision (%lu b, PRECISION_DP defined)\n", sizeof(PdfT));
#elif defined(PRECISION_SP)
	printf("# - floating point:    single precision (%lu b, PRECISION_SP defined)\n", sizeof(PdfT));
#else
	printf("# - floating point:    UNKNOWN (%lu b)\n", sizeof(PdfT));
#endif

#ifdef VECTOR_AVX
	printf("# - intrinsics:        AVX (VECTOR_AVX defined)\n");
#elif defined(VECTOR_SSE)
	printf("# - intrinsics:        SSE (VECTOR_SSE defined)\n");
#else
	printf("# - intrinsics:        UNKNOWN\n");
#endif

	printf("# - defines:           ");
	for (int j = 0; j < N_ELEMS(defines); ++j) {
		printf("%s ", defines[j]);
	}
	printf("\n");

#ifdef __x86_64__
	printf("# - fp status:         DAZ: %d  FTZ: %d\n", FpGetDaz(), FpGetFtz());
#endif

	printf("# - iterations:        %d\n", cd.MaxIterations);

 	LatticeDesc ld;

	GeoCreateByStr(geometryType, dims, periodic, &ld);

	printf("# - geometry:\n");
	printf("#   type:              %s\n", ld.Name);
	printf("#   dimensions:        %d x %d x %d (x, y, z)\n", ld.Dims[0], ld.Dims[1], ld.Dims[2]);

	printf("#   nodes total:       %d\n", ld.nObst + ld.nFluid);
	printf("#   nodes fluid:       %d (including inlet & outlet)\n", ld.nFluid);
	printf("#   nodes obstacles:   %d\n", ld.nObst);
	printf("#   nodes inlet:       %d\n", ld.nInlet);
	printf("#   nodes outlet:      %d\n", ld.nOutlet);
	printf("#   periodicity:       x: %d y: %d z: %d\n", ld.PeriodicX, ld.PeriodicY, ld.PeriodicZ);

#ifdef VTK_OUTPUT
	printf("# - VTK output:        %d (every %d iteration)\n", cd.VtkOutput, cd.VtkModulus);
#endif
#ifdef STATISTICS
	printf("# - statistics:        every %d iteration\n", cd.StatisticsModulus);
#endif

	printf("# - flow:\n");
	printf("#   omega:             %f\n", cd.Omega);
	printf("#   initial density at inlet/outlet:\n");
	printf("#     rho in:          %e\n", cd.RhoIn);
    printf("#     rho out:         %e\n", cd.RhoOut);

#ifdef _OPENMP
	printf("# - OpenMP threads:    %d\n", omp_get_max_threads());

	if (pinString != NULL) {
		#pragma omp parallel
		{
			int threadId = omp_get_thread_num();
			int err;

			err = PinCurrentThreadByCpuList(pinString, threadId);

			if (err) {
				printf("ERROR [thread %d]: pinning failed.\n", threadId);
				exit(1);
			}

			const char * cpuList = PinCpuListAsString();
			Assert(cpuList != NULL);

			// Not so nice hack to print the thread ids ordered.
			#pragma omp for ordered
			for (int i = 0; i < omp_get_num_threads(); ++i) {
				#pragma omp ordered
				printf("#   thread %2d  pinned to core(s):  %s\n", threadId, cpuList);
			}

			free((void *)cpuList);
		}
	}
#endif

	KernelData * kd;

	KernelFunctions * kf = NULL;

	if (kernelToUse == NULL) {
		kf = &g_kernels[0];
	}
	else {
		for (int j = 0; j < N_ELEMS(g_kernels); ++j) {

			if (!strcasecmp(kernelToUse, g_kernels[j].Name)) {
				kf = &g_kernels[j];
				break;
			}
		}
	}

	if (kf == NULL) {
		printf("ERROR: requested kernel \"%s\" not found.\n", kernelToUse);
		exit(1);
	}

	printf("#\n");
	printf("# - kernel:            %s\n", kf->Name);
	printf("#\n");

	// Initialize kernel by calling its own initialization function
	kf->Init(&ld, &kd, &p);

#ifdef VERIFICATION
	if (verify) {
		KernelSetInitialDensity( &ld, kd, &cd);
		KernelSetInitialVelocity(&ld, kd, &cd);
	}
#endif

	printf("# starting kernel...\n");

	X_LIKWID_INIT();

	double timeStart = Time();

	// Call the LBM kernel
	kd->Kernel(&ld, kd, &cd);

	double duration = Time() - timeStart;

	X_LIKWID_DEINIT();

	// Print some statistics...
	KernelStatisticsAdv(kd, &ld, &cd, cd.MaxIterations, 1 /* force output */);

#ifdef VERIFICATION
	PdfT errorNorm = -1.0;
	KernelVerifiy(&ld, kd, &cd, &errorNorm);
#endif

	// Deinitialize kernel by calling its own deinitialization function
	kf->Deinit(&ld, &kd);


	double perf = (double)ld.nFluid * (double)cd.MaxIterations / duration / 1.e6;

	printf("P:   %f MFLUP/s  t: %d  d: %f s  iter: %d  fnodes: %f x1e6  geo: %s  kernel: %s  %s  %s\n",
		perf, nThreads, duration, cd.MaxIterations, ld.nFluid / 1e6,
		geometryType, kernelToUse,
#ifdef VERIFICATION
		"VERIFICATION",
#else
		"B",
#endif
#ifdef PRECISION_DP
		"dp"
#elif defined(PRECISION_SP)
		"sp"
#else
		"unknown-precision"
#endif
	);

	int exitCode = 0;

#ifdef VERIFICATION

	if (verify) {
		printf("# VERIFICATION: deviation from analytical solution: %e\n", errorNorm);

		if (errorNorm > 0.1) {
			printf("# VERIFICATION FAILED.\n");
			exitCode = 1;
		}
		else {
			printf("# VERIFICATION SUCCEEDED.\n");
		}
	}
#else
//	printf("# VERIFICATION: deviation from analytical solution: %e\n", errorNorm);
//	printf("# VERIFICATION: this is only valid for pipe geometry with enough iterations performed.\n");
#endif

	MemFree((void **)&ld.Lattice);

	return exitCode;
}
