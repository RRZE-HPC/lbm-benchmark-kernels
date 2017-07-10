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
#include "Vtk.h"

#include <math.h>

// TODO: make this portable

// needed for stat & mkdir
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <string.h> // strerror

// TODO: make byteswap portable

#include <inttypes.h>
// glibc
#include <byteswap.h>

// macros for portability
// #define BS32(a)		bswap_32(*((uint32_t *)(&a)))
#define BS64(a)		bswap_64(*((uint64_t *)(&a)))


void VtkWrite(LatticeDesc * ld, KernelData * kd, CaseData * cd, int iteration)
{
	Assert(kd != NULL);
	Assert(ld != NULL);
	Assert(ld->Dims[0] > 0);
	Assert(ld->Dims[1] > 0);
	Assert(ld->Dims[2] > 0);

	// TODO: this should be made portable...
	// Check if subdirectory vtk exists, if not, create it.
	{
		int err;
		struct stat fileStatus;

		err = stat("vtk", &fileStatus);

		if (err) {
			// printf("ERROR: stat %d - %s\n", errno, strerror(errno));

			// Set default mask and hope mkdir applies umask...
			err = mkdir("vtk", 0700);

			if (err) {
				printf("ERROR: cannot create directory vtk - %d: %s\n", errno, strerror(errno));
				exit(1);
			}

			printf("# created directory vtk.\n");
		}
		else {

			if (!S_ISDIR(fileStatus.st_mode)) {
				printf("ERROR: cannot create subdirectory vtk as already a file with the same name exists.\n");
				exit(1);
			}

		}
	}


	char fileName[1024];

	snprintf(fileName, sizeof(fileName), "vtk/file-%04d.vtk", iteration);

	printf("# VTK: writing file %s\n", fileName);

	FILE * fh;

	fh = fopen(fileName, "w");

	if(fh == NULL) {
		printf("ERROR: opening file %s failed.\n", fileName);
		exit(1);
	}

	// http://www.vtk.org/pdf/file-formats.pdf
	int nX = ld->Dims[0];
	int nY = ld->Dims[1];
	int nZ = ld->Dims[2];
	int * lDims = ld->Dims;

	// Temporaries for endian conversion.
	uint64_t uDensity, uUx, uUy, uUz;

	PdfT pdfs[N_D3Q19];

	fprintf(fh, "# vtk DataFile Version 1.0\n");
	fprintf(fh, "Comment: lid driven cavity, iteration % 4d\n", iteration);
#ifdef VTK_OUTPUT_ASCII
	fprintf(fh, "ASCII\n");
#else
	fprintf(fh, "BINARY\n");
#endif
	fprintf(fh, "DATASET STRUCTURED_POINTS\n");
	fprintf(fh, "DIMENSIONS %d %d %d\n", nX, nY, nZ);
	fprintf(fh, "ORIGIN 0 0 0 \n");
	fprintf(fh, "SPACING 1 1 1\n");
	fprintf(fh, "POINT_DATA %d\n", nX * nY * nZ);

	// ----------------------------------------------------------------------
	// Flag field: obstacle = 0, fluid = 1, inlet = 2, outlet = 4

	fprintf(fh, "SCALARS NodesTypes unsigned_char 1\n");
	fprintf(fh, "LOOKUP_TABLE default\n");

	unsigned char c;

	for(int z = 0; z < nZ; ++z) {
		for(int y = 0; y < nY; ++y) {
			for(int x = 0; x < nX; ++x) {
#ifdef VTK_OUTPUT_ASCII
				fprintf(fh, "%d\n", ld->Lattice[L_INDEX_4(ld->Dims, x, y, z)]);
#else
				c = (unsigned char)ld->Lattice[L_INDEX_4(ld->Dims, x, y, z)];
				fwrite(&c, sizeof(unsigned char), 1, fh);
#endif
			}
		}
	}

	// ----------------------------------------------------------------------
	// Density field

	fprintf(fh, "SCALARS Density double\n");
	fprintf(fh, "LOOKUP_TABLE default\n");

	double density;

	for(int z = 0; z < nZ; ++z) {
		for(int y = 0; y < nY; ++y) {
			for(int x = 0; x < nX; ++x) {

				density = 0.0;
				if (ld->Lattice[L_INDEX_4(lDims, x, y, z)] != LAT_CELL_OBSTACLE) {
					kd->GetNode(kd, x, y, z, pdfs);

					for (int d = 0; d < N_D3Q19; ++d) {
						density += pdfs[d];
					}
				}

#ifdef VTK_OUTPUT_ASCII
				fprintf(fh, "%e\n", density);
#else
				uDensity = BS64(density);
				fwrite(&uDensity, sizeof(double), 1, fh);
#endif
			}
		}
	}

	// ----------------------------------------------------------------------
	// Velocity vectors: velocity in x, y, and z direction

	fprintf(fh, "VECTORS VelocityVectors double\n");

	// Declare pdf_N, pdf_E, pdf_S, pdf_W, ...
	#define X(name, idx, idxinv, x, y, z)	PdfT JOIN(pdf_,name);
	D3Q19_LIST
	#undef X

	double ux, uy, uz;

	for(int z = 0; z < nZ; ++z) {
		for(int y = 0; y < nY; ++y) {
			for(int x = 0; x < nX; ++x) {

				if (ld->Lattice[L_INDEX_4(lDims, x, y, z)] != LAT_CELL_OBSTACLE) {
					kd->GetNode(kd, x, y, z, pdfs);

// DETECT NANS
//					for (int d = 0; d < 19; ++d) {
//						if(isnan(pdfs[d])) {
//							printf("%d %d %d %d nan!\n", x, y, z, d);
//							for (int d2 = 0; d2 < 19; ++d2) {
//								printf("%d: %e\n", d2, pdfs[d2]);
//							}
//							exit(1);
//						}
//					}
					#define X(name, idx, idxinv, _x, _y, _z)	JOIN(pdf_,name) = pdfs[idx];
					D3Q19_LIST
					#undef X
					UNUSED(pdf_C);


					ux = pdf_E + pdf_NE + pdf_SE + pdf_TE + pdf_BE -
						 pdf_W - pdf_NW - pdf_SW - pdf_TW - pdf_BW;
					uy = pdf_N + pdf_NE + pdf_NW + pdf_TN + pdf_BN -
						 pdf_S - pdf_SE - pdf_SW - pdf_TS - pdf_BS;
					uz = pdf_T + pdf_TE + pdf_TW + pdf_TN + pdf_TS -
						 pdf_B - pdf_BE - pdf_BW - pdf_BN - pdf_BS;
					#ifdef VERIFICATION
					ux += 0.5 * cd->XForce;
					#endif
				}
				else {
					ux = 0.0; uy = 0.0; uz = 0.0;
				}

#ifdef VTK_OUTPUT_ASCII
				fprintf(fh, "%f %f %f\n", ux, uy, uz);
#else
				uUx = BS64(ux); uUy = BS64(uy); uUz = BS64(uz);
				fwrite(&uUx, sizeof(double), 1, fh);
				fwrite(&uUy, sizeof(double), 1, fh);
				fwrite(&uUz, sizeof(double), 1, fh);
#endif
			}
		}
	}

	fclose(fh);
}

