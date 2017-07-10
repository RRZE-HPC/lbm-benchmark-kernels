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
#include "Kernel.h"
#include "Lattice.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#define X(name, idx, idx_inv, x, y, z)	, x
int D3Q19_X[] = {
	EXPAND(D3Q19_LIST)
};
#undef X

#define X(name, idx, idx_inv, x, y, z)	, y
int D3Q19_Y[] = {
	EXPAND(D3Q19_LIST)
};
#undef X

#define X(name, idx, idx_inv, x, y, z)	, z
int D3Q19_Z[] = {
	EXPAND(D3Q19_LIST)
};
#undef X

#define X(name, idx, idxinv, x, y, z)	, idxinv
int D3Q19_INV[] = {
	EXPAND(D3Q19_LIST)
};
#undef X


#define X(name, idx, idxinv, x, y, z)	, STRINGIFY(name)
const char * D3Q19_NAMES[N_D3Q19] = {
	EXPAND(D3Q19_LIST)
};
#undef X

void KernelComputeBoundaryConditions(KernelData * kd, LatticeDesc * ld, CaseData * cd)
{
	Assert(kd != NULL);
	Assert(ld != NULL);
	Assert(cd != NULL);

	Assert(cd->RhoIn  > 0.0);
	Assert(cd->RhoOut > 0.0);

	PdfT rho_in         = cd->RhoIn;
	PdfT rho_out        = cd->RhoOut;
	PdfT rho_in_inv     = 1.0 / rho_in;
	PdfT rho_out_inv    = 1.0 / rho_out;
	PdfT indep_ux       = 0.0;

	PdfT dens;
	PdfT ux;

	const PdfT one_third  = 1.0 / 3.0;
	const PdfT one_fourth = 1.0 / 4.0;
	const PdfT one_sixth  = 1.0 / 6.0;

	PdfT pdfs[N_D3Q19];

	int nX = kd->Dims[0];
	int nY = kd->Dims[1];
	int nZ = kd->Dims[2];

	int x;
	int x_in  = 0;
	int x_out = nX - 1;

	double density_in = 0.0;
	double density_out = 0.0;

	// update inlet / outlet boundary conditions
	for (int z = 1; z < nZ - 1; ++z) {
		for (int y = 1; y < nY - 1; ++y) {


			// -----------------------------------------------------------------------------
			// update inlet conditions

			if (ld->Lattice[L_INDEX_4(ld->Dims, x_in, y, z)] == LAT_CELL_INLET) {

				x = x_in;

				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_C , pdfs + D3Q19_C);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_T , pdfs + D3Q19_T);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_B , pdfs + D3Q19_B);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_S , pdfs + D3Q19_S);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_N , pdfs + D3Q19_N);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_TS, pdfs + D3Q19_TS);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_BS, pdfs + D3Q19_BS);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_TN, pdfs + D3Q19_TN);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_BN, pdfs + D3Q19_BN);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_SW, pdfs + D3Q19_SW);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_TW, pdfs + D3Q19_TW);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_W , pdfs + D3Q19_W);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_BW, pdfs + D3Q19_BW);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_NW, pdfs + D3Q19_NW);

				dens = rho_in;

				ux = 1 - (pdfs[D3Q19_C]  +
						(pdfs[D3Q19_T]  + pdfs[D3Q19_B]  + pdfs[D3Q19_S]  + pdfs[D3Q19_N]) +
						(pdfs[D3Q19_TS] + pdfs[D3Q19_BS] + pdfs[D3Q19_TN] + pdfs[D3Q19_BN]) +
						2 * (pdfs[D3Q19_SW] + pdfs[D3Q19_TW] + pdfs[D3Q19_W] + pdfs[D3Q19_BW] + pdfs[D3Q19_NW])) * rho_in_inv;

				indep_ux = one_sixth * dens * ux;

				pdfs[D3Q19_E ] = pdfs[D3Q19_W]  + one_third  * dens * ux;
				pdfs[D3Q19_NE] = pdfs[D3Q19_SW] - one_fourth * (pdfs[D3Q19_N] - pdfs[D3Q19_S]) + indep_ux;
				pdfs[D3Q19_SE] = pdfs[D3Q19_NW] + one_fourth * (pdfs[D3Q19_N] - pdfs[D3Q19_S]) + indep_ux;
				pdfs[D3Q19_TE] = pdfs[D3Q19_BW] - one_fourth * (pdfs[D3Q19_T] - pdfs[D3Q19_B]) + indep_ux;
				pdfs[D3Q19_BE] = pdfs[D3Q19_TW] + one_fourth * (pdfs[D3Q19_T] - pdfs[D3Q19_B]) + indep_ux;


				kd->BoundaryConditionsSetPdf(kd, x, y, z, D3Q19_E , pdfs[D3Q19_E ]);
				kd->BoundaryConditionsSetPdf(kd, x, y, z, D3Q19_NE, pdfs[D3Q19_NE]);
				kd->BoundaryConditionsSetPdf(kd, x, y, z, D3Q19_SE, pdfs[D3Q19_SE]);
				kd->BoundaryConditionsSetPdf(kd, x, y, z, D3Q19_TE, pdfs[D3Q19_TE]);
				kd->BoundaryConditionsSetPdf(kd, x, y, z, D3Q19_BE, pdfs[D3Q19_BE]);

				for(int d = 0; d < N_D3Q19; ++d) {
					density_in += pdfs[d];
				}
			}

			// -----------------------------------------------------------------------------
			// update outlet conditions

			if (ld->Lattice[L_INDEX_4(ld->Dims, x_out, y, z)] == LAT_CELL_OUTLET) {
				// update outlet conditions

				x = x_out;

				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_C , pdfs + D3Q19_C );
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_T , pdfs + D3Q19_T );
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_B , pdfs + D3Q19_B );
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_S , pdfs + D3Q19_S );
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_N , pdfs + D3Q19_N );
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_TS, pdfs + D3Q19_TS);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_BS, pdfs + D3Q19_BS);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_TN, pdfs + D3Q19_TN);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_BN, pdfs + D3Q19_BN);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_NE, pdfs + D3Q19_NE);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_BE, pdfs + D3Q19_BE);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_E , pdfs + D3Q19_E );
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_TE, pdfs + D3Q19_TE);
				kd->BoundaryConditionsGetPdf(kd, x, y, z, D3Q19_SE, pdfs + D3Q19_SE);

				dens = rho_out;

				ux = -1 + (pdfs[D3Q19_C] +
						(pdfs[D3Q19_T]  + pdfs[D3Q19_B]  + pdfs[D3Q19_S]  + pdfs[D3Q19_N]) +
						(pdfs[D3Q19_TS] + pdfs[D3Q19_BS] + pdfs[D3Q19_TN] + pdfs[D3Q19_BN]) +
						2 * (pdfs[D3Q19_NE] + pdfs[D3Q19_BE] + pdfs[D3Q19_E] + pdfs[D3Q19_TE] + pdfs[D3Q19_SE])) * rho_out_inv;
				indep_ux = one_sixth * dens * ux;

				pdfs[D3Q19_W ] = pdfs[D3Q19_E] - one_third * dens * ux;
				pdfs[D3Q19_SW] = pdfs[D3Q19_NE] + one_fourth * (pdfs[D3Q19_N] - pdfs[D3Q19_S]) - indep_ux;
				pdfs[D3Q19_NW] = pdfs[D3Q19_SE] - one_fourth * (pdfs[D3Q19_N] - pdfs[D3Q19_S]) - indep_ux;
				pdfs[D3Q19_BW] = pdfs[D3Q19_TE] + one_fourth * (pdfs[D3Q19_T] - pdfs[D3Q19_B]) - indep_ux;
				pdfs[D3Q19_TW] = pdfs[D3Q19_BE] - one_fourth * (pdfs[D3Q19_T] - pdfs[D3Q19_B]) - indep_ux;

				kd->BoundaryConditionsSetPdf(kd, x, y, z, D3Q19_W , pdfs[D3Q19_W ]);
				kd->BoundaryConditionsSetPdf(kd, x, y, z, D3Q19_NW, pdfs[D3Q19_NW]);
				kd->BoundaryConditionsSetPdf(kd, x, y, z, D3Q19_SW, pdfs[D3Q19_SW]);
				kd->BoundaryConditionsSetPdf(kd, x, y, z, D3Q19_TW, pdfs[D3Q19_TW]);
				kd->BoundaryConditionsSetPdf(kd, x, y, z, D3Q19_BW, pdfs[D3Q19_BW]);

				for(int d = 0; d < N_D3Q19; ++d) {
					density_out += pdfs[d];
				}
			}
		}
	}

	// DEBUG: printf("# density inlet: %e  density outlet: %e\n", density_in, density_out);

}


PdfT KernelDensity(KernelData * kd, LatticeDesc * ld)
{
	Assert(kd != NULL);
	Assert(ld != NULL);

	Assert(ld->Lattice != NULL);
	Assert(ld->Dims    != NULL);

	Assert(ld->Dims[0] > 0);
	Assert(ld->Dims[1] > 0);
	Assert(ld->Dims[2] > 0);

	int * lDims = ld->Dims;
	int nX = lDims[0];
	int nY = lDims[1];
	int nZ = lDims[2];

	PdfT pdfs[N_D3Q19] = { -1.0 };
	PdfT density = 0.0;

	for(int z = 0; z < nZ; ++z) {
		for(int y = 0; y < nY; ++y) {
			for(int x = 0; x < nX; ++x) {

				if(ld->Lattice[L_INDEX_4(lDims, x, y, z)] != LAT_CELL_OBSTACLE) {

					kd->GetNode(kd, x, y, z, pdfs);

					for(int d = 0; d < N_D3Q19; ++d) {
// 						if (pdfs[d] < 0.0) {
// 							printf("# %d %d %d %d < 0 %e %s\n", x, y, z, d, pdfs[d], D3Q19_NAMES[d]);
// 	exit(1);
// 						}
						density += pdfs[d];
					}
				}

			}
		}
	}

	return density / ld->nFluid;
}


// prescribes a given density
void KernelSetInitialDensity(LatticeDesc * ld, KernelData * kd, CaseData * cd)
{
	int * lDims = ld->Dims;

	PdfT rho_in = cd->RhoIn;
	PdfT rho_out = cd->RhoOut;

	PdfT ux = 0.0;
	PdfT uy = 0.0;
	PdfT uz = 0.0;
	PdfT dens = 1.0;

	PdfT omega = cd->Omega;

	PdfT w_0 = 1.0 /  3.0;
	PdfT w_1 = 1.0 / 18.0;
	PdfT w_2 = 1.0 / 36.0;

	PdfT dir_indep_trm;
	PdfT omega_w0  = 3.0 * w_0 * omega;
	PdfT omega_w1  = 3.0 * w_1 * omega;
	PdfT omega_w2  = 3.0 * w_2 * omega;
	PdfT one_third = 1.0 / 3.0;

	int nX = lDims[0];
	int nY = lDims[1];
	int nZ = lDims[2];

	PdfT pdfs[N_D3Q19];

	#ifdef _OPENMP
		#pragma omp parallel for collapse(3)
	#endif
	for(int z = 0; z < nZ; ++z) { for(int y = 0; y < nY; ++y) { for(int x = 0; x < nX; ++x) {

		if (ld->Lattice[L_INDEX_4(ld->Dims, x, y, z)] != LAT_CELL_OBSTACLE) {
	// TODO: fix later.
	//		if((caseData->geoType == GEO_TYPE_CHANNEL) || (caseData->geoType == GEO_TYPE_RCHANNEL))
				dens = rho_in + (rho_out - rho_in)*(x)/(nX-1.0);

			#define SQR(a) ((a)*(a))

			dir_indep_trm = one_third * dens - 0.5 * (ux * ux + uy * uy + uz * uz);

			pdfs[D3Q19_C]  = omega_w0 * (dir_indep_trm);

			pdfs[D3Q19_NW] = omega_w2 * (dir_indep_trm - (ux - uy) + 1.5 * SQR(ux - uy));
			pdfs[D3Q19_SE] = omega_w2 * (dir_indep_trm + (ux - uy) + 1.5 * SQR(ux - uy));

			pdfs[D3Q19_NE] = omega_w2 * (dir_indep_trm + (ux + uy) + 1.5 * SQR(ux + uy));
			pdfs[D3Q19_SW] = omega_w2 * (dir_indep_trm - (ux + uy) + 1.5 * SQR(ux + uy));


			pdfs[D3Q19_TW] = omega_w2 * (dir_indep_trm - (ux - uz) + 1.5 * SQR(ux - uz));
			pdfs[D3Q19_BE] = omega_w2 * (dir_indep_trm + (ux - uz) + 1.5 * SQR(ux - uz));

			pdfs[D3Q19_TE] = omega_w2 * (dir_indep_trm + (ux + uz) + 1.5 * SQR(ux + uz));
			pdfs[D3Q19_BW] = omega_w2 * (dir_indep_trm - (ux + uz) + 1.5 * SQR(ux + uz));


			pdfs[D3Q19_TS] = omega_w2 * (dir_indep_trm - (uy - uz) + 1.5 * SQR(uy - uz));
			pdfs[D3Q19_BN] = omega_w2 * (dir_indep_trm + (uy - uz) + 1.5 * SQR(uy - uz));

			pdfs[D3Q19_TN] = omega_w2 * (dir_indep_trm + (uy + uz) + 1.5 * SQR(uy + uz));
			pdfs[D3Q19_BS] = omega_w2 * (dir_indep_trm - (uy + uz) + 1.5 * SQR(uy + uz));


			pdfs[D3Q19_N]  = omega_w1 * (dir_indep_trm + uy + 1.5 * SQR(uy));
			pdfs[D3Q19_S]  = omega_w1 * (dir_indep_trm - uy + 1.5 * SQR(uy));

			pdfs[D3Q19_E]  = omega_w1 * (dir_indep_trm + ux + 1.5 * SQR(ux));
			pdfs[D3Q19_W]  = omega_w1 * (dir_indep_trm - ux + 1.5 * SQR(ux));

			pdfs[D3Q19_T]  = omega_w1 * (dir_indep_trm + uz + 1.5 * SQR(uz));
			pdfs[D3Q19_B]  = omega_w1 * (dir_indep_trm - uz + 1.5 * SQR(uz));


			kd->SetNode(kd, x, y, z, pdfs);

			#undef SQR
		}
	} } }
}


// prescribes a given velocity
void KernelSetInitialVelocity(LatticeDesc * ld, KernelData * kd, CaseData * cd)
{

	int * lDims = ld->Dims;

	// TODO: ux is overriden below...
	PdfT ux = 0.09; // caseData->initUx;
	PdfT uy = 0.0; // caseData->initUy;
	PdfT uz = 0.0; // caseData->initUz;
	PdfT dens = 1.0;

	PdfT omega = cd->Omega;

	PdfT w_0 = 1.0 /  3.0;
	PdfT w_1 = 1.0 / 18.0;
	PdfT w_2 = 1.0 / 36.0;

	PdfT dir_indep_trm;
	PdfT omega_w0  = 3.0 * w_0 * omega;
	PdfT omega_w1  = 3.0 * w_1 * omega;
	PdfT omega_w2  = 3.0 * w_2 * omega;
	PdfT one_third = 1.0 / 3.0;

	int nX = lDims[0];
	int nY = lDims[1];
	int nZ = lDims[2];

	PdfT pdfs[N_D3Q19];

	PdfT density;

	#ifdef _OPENMP
		#pragma omp parallel for collapse(3)
	#endif
	for(int z = 0; z < nZ; ++z) { for(int y = 0; y < nY; ++y) { for(int x = 0; x < nX; ++x) {

		if (ld->Lattice[L_INDEX_4(ld->Dims, x, y, z)] == LAT_CELL_FLUID) {

			ux = 0.0;
			uy = 0.0;
			uz = 0.0;

			kd->GetNode(kd, x, y, z, pdfs);


			density = 0.0;

			#define X(name, idx, idxinv, _x, _y, _z)	density += pdfs[idx];
				D3Q19_LIST
			#undef X


			#define SQR(a) ((a)*(a))
			dir_indep_trm = one_third * dens - 0.5 * (ux * ux + uy * uy + uz * uz);

			pdfs[D3Q19_C]  = omega_w0 * (dir_indep_trm);

			pdfs[D3Q19_NW] = omega_w2 * (dir_indep_trm - (ux - uy) + 1.5 * SQR(ux - uy));
			pdfs[D3Q19_SE] = omega_w2 * (dir_indep_trm + (ux - uy) + 1.5 * SQR(ux - uy));

			pdfs[D3Q19_NE] = omega_w2 * (dir_indep_trm + (ux + uy) + 1.5 * SQR(ux + uy));
			pdfs[D3Q19_SW] = omega_w2 * (dir_indep_trm - (ux + uy) + 1.5 * SQR(ux + uy));


			pdfs[D3Q19_TW] = omega_w2 * (dir_indep_trm - (ux - uz) + 1.5 * SQR(ux - uz));
			pdfs[D3Q19_BE] = omega_w2 * (dir_indep_trm + (ux - uz) + 1.5 * SQR(ux - uz));

			pdfs[D3Q19_TE] = omega_w2 * (dir_indep_trm + (ux + uz) + 1.5 * SQR(ux + uz));
			pdfs[D3Q19_BW] = omega_w2 * (dir_indep_trm - (ux + uz) + 1.5 * SQR(ux + uz));


			pdfs[D3Q19_TS] = omega_w2 * (dir_indep_trm - (uy - uz) + 1.5 * SQR(uy - uz));
			pdfs[D3Q19_BN] = omega_w2 * (dir_indep_trm + (uy - uz) + 1.5 * SQR(uy - uz));

			pdfs[D3Q19_TN] = omega_w2 * (dir_indep_trm + (uy + uz) + 1.5 * SQR(uy + uz));
			pdfs[D3Q19_BS] = omega_w2 * (dir_indep_trm - (uy + uz) + 1.5 * SQR(uy + uz));


			pdfs[D3Q19_N]  = omega_w1 * (dir_indep_trm + uy + 1.5 * SQR(uy));
			pdfs[D3Q19_S]  = omega_w1 * (dir_indep_trm - uy + 1.5 * SQR(uy));

			pdfs[D3Q19_E]  = omega_w1 * (dir_indep_trm + ux + 1.5 * SQR(ux));
			pdfs[D3Q19_W]  = omega_w1 * (dir_indep_trm - ux + 1.5 * SQR(ux));

			pdfs[D3Q19_T]  = omega_w1 * (dir_indep_trm + uz + 1.5 * SQR(uz));
			pdfs[D3Q19_B]  = omega_w1 * (dir_indep_trm - uz + 1.5 * SQR(uz));

			#undef SQR


			kd->SetNode(kd, x, y, z, pdfs);
		}
	} } }

}

// Compute analytical x velocity for channel flow.
//
// Formula 7 from Kutay et al. "Laboratory validation of lattice Boltzmann method for modeling
// pore-scale flow in granular materials", doi:10.1016/j.compgeo.2006.08.002.
//
// also formula 10 from
// Pan et al. "An evaluation of lattice Boltzmann equation methods for simulating flow
// through porous media", doi:10.1016/S0167-5648(04)80040-6.
//
// calculate velocity in a pipe for a given radius
//
static PdfT CalcXVelForPipeProfile(PdfT maxRadiusSquared, PdfT curRadiusSquared, PdfT xForce, PdfT viscosity)
{
	return xForce*(maxRadiusSquared - curRadiusSquared) / (2.0*viscosity);
}

static void KernelGetXSlice(LatticeDesc * ld, KernelData * kd, CaseData * cd, PdfT * outputArray, int xPos)
{
	Assert(ld != NULL);
	Assert(kd != NULL);

	int nY = ld->Dims[1];
	int nZ = ld->Dims[2];

	Assert(xPos >= 0);
	Assert(xPos < ld->Dims[0]);


	PdfT ux = 0.0;

	// Declare pdf_N, pdf_E, pdf_S, pdf_W, ...
	#define X(name, idx, idxinv, x, y, z)	PdfT JOIN(pdf_,name);
	D3Q19_LIST
	#undef X
	PdfT pdfs[N_D3Q19];

	for(int z = 0; z < nZ; ++z) {
		for(int y = 0; y < nY; ++y) {

			if (ld->Lattice[L_INDEX_4(ld->Dims, xPos, y, z)] != LAT_CELL_OBSTACLE) {
				kd->GetNode(kd, xPos, y, z, pdfs);

				#define X(name, idx, idxinv, _x, _y, _z)	JOIN(pdf_,name) = pdfs[idx];
				D3Q19_LIST
				#undef X
				UNUSED(pdf_C); UNUSED(pdf_S); UNUSED(pdf_N); UNUSED(pdf_T); UNUSED(pdf_B);
				UNUSED(pdf_TN); UNUSED(pdf_BN); UNUSED(pdf_TS); UNUSED(pdf_BS);

				ux = pdf_E + pdf_NE + pdf_SE + pdf_TE + pdf_BE -
					 pdf_W - pdf_NW - pdf_SW - pdf_TW - pdf_BW;

				#ifdef VERIFICATION
				ux += 0.5 * cd->XForce;
				#endif

				outputArray[y * nZ + z] = ux;
			}
			else {
				outputArray[y * nZ + z] = 0.0;
			}
		}
	}

}

// Verification of channel profile with analytical solution.
// Taken from Kutay et al. "Laboratory validation of lattice Boltzmann method for modeling
// pore-scale flow in granular materials", doi:10.1016/j.compgeo.2006.08.002. and
// Pan et al. "An evaluation of lattice Boltzmann equation methods for simulating flow
// through porous media", doi:10.1016/S0167-5648(04)80040-6
//
void KernelVerifiy(LatticeDesc * ld, KernelData * kd, CaseData * cd, PdfT * errorNorm)
{
	Assert(ld != NULL);
	Assert(kd != NULL);
	Assert(cd != NULL);
	Assert(errorNorm != NULL);

	int nX = ld->Dims[0];
	int nY = ld->Dims[1];
	int nZ = ld->Dims[2];

	PdfT omega   = cd->Omega;
	PdfT viscosity        = (1.0 / omega - 0.5) / 3.0;

	// ux averaged across cross sections in x direction
	PdfT * outputArray = (PdfT *)malloc(nZ * nY * sizeof(PdfT));
	Verify(outputArray != NULL);

	memset(outputArray, -10, nZ*nY*sizeof(PdfT));

	// uncomment this to get values averaged along the x-axis
	//AveragePipeCrossSections(ld, kd, outputArray);
	KernelGetXSlice(ld, kd, cd, outputArray, (int)(nX/2));


	FILE * fh;
	char fileName[1024];
	PdfT tmpAvgUx = 0.0;
	PdfT tmpAnalyUx = 0.0;
	int flagEvenNy = 0;
	int y = 0;

	if (nY % 2 == 0)
		flagEvenNy = 1;

	y = (nY-flagEvenNy-1)/2;

	snprintf(fileName, sizeof(fileName), "flow-profile.dat");

	printf("# Kernel validation: writing profile to %s\n", fileName);

	fh = fopen(fileName, "w");

	if(fh == NULL) {
		printf("ERROR: opening file %s failed.\n", fileName);
		exit(1);
	}

	fprintf(fh, "# Flow profile in Z direction. Taken at the middle of the X length (= %d) of total length %d.\n", nZ / 2, nZ);
	// fprintf(fh, "# Snapshot taken at iteration %d.\n", iteration);
	fprintf(fh, "# Plot on terminal: gnuplot -e \"set terminal dumb; plot \\\"%s\\\" u 1:3 t \\\"analytical\\\", \\\"\\\" u 1:4 t \\\"simulation\\\";\"\n", fileName);
	fprintf(fh, "# Plot graphically: gnuplot -e \"plot \\\"%s\\\" u 1:3 w linesp t \\\"analytical\\\", \\\"\\\" u 1:4 w linesp t \\\"simulation\\\"; pause -1;\"\n", fileName);
	fprintf(fh, "# z coord., radius, analytic, simulation, diff abs, diff rel, undim_analytic, undim_sim\n");

	double deviation = 0.0;
	double curRadiusSquared;
	double center = nY / 2.0;
	double minDiameter = nY;
	#define SQR(a) ((a)*(a))
	double minRadiusSquared = SQR(minDiameter / 2.0 - 1.0);
	#undef SQR
	double u_max = cd->XForce*minRadiusSquared/(2.0*viscosity);

	for(int z = 0; z < nZ; ++z) {

		fprintf(fh, "%d\t", z);

		#define SQR(a) ((a)*(a))
		curRadiusSquared = SQR(z-center+0.5);


		// dimensionless radius
		fprintf(fh, "%e\t", (z-center+0.5)/center);

		// analytic profile
		if(curRadiusSquared >= minRadiusSquared)
			tmpAnalyUx = 0.0;
		else
			tmpAnalyUx = CalcXVelForPipeProfile(minRadiusSquared, curRadiusSquared, cd->XForce, viscosity);

		//averaged profile
		if(flagEvenNy == 1)
			tmpAvgUx = (outputArray[y*nZ + z] + outputArray[(y+1)*nZ + z])/2.0;
		else
			tmpAvgUx = outputArray[y*nZ + z];

		fprintf(fh, "%e\t", tmpAnalyUx);
		fprintf(fh, "%e\t", tmpAvgUx);

		fprintf(fh, "%e\t", fabs(tmpAnalyUx-tmpAvgUx));
		if (tmpAnalyUx != 0.0) {
			fprintf(fh, "%e\t", fabs(tmpAnalyUx - tmpAvgUx) / tmpAnalyUx);
			deviation += SQR(fabs(tmpAnalyUx - tmpAvgUx) / tmpAnalyUx);
		}
		else {
			fprintf(fh, "0.0\t");
		}

		fprintf(fh, "%e\t", tmpAnalyUx / u_max);
		fprintf(fh, "%e\t", tmpAvgUx / u_max);
		fprintf(fh, "\n");

		#undef SQR
	}

	*errorNorm = sqrt(deviation);

	printf("# Kernel validation: L2 error norm of relative error: %e\n", *errorNorm);


	fclose(fh);
	free(outputArray);


}


void KernelStatistics(KernelData * kd, LatticeDesc * ld, CaseData * cd, int iteration)
{
	KernelStatisticsAdv(kd, ld, cd, iteration, 0);
}

void KernelStatisticsAdv(KernelData * kd, LatticeDesc * ld, CaseData * cd, int iteration, int forceOutput)
{
	if (iteration % cd->StatisticsModulus == 0 || forceOutput) {
		printf("# iter: %4d   avg density: %e\n", iteration, KernelDensity(kd, ld));
	}

	if (iteration % 10 != 0 && !forceOutput) {
		return;
	}

	int nX = ld->Dims[0];
	int nY = ld->Dims[1];
	int nZ = ld->Dims[2];

	int x = nX / 2;

	PdfT pdfs[N_D3Q19];

	// ----------------------------------------------------------------------
	// velocity in x-direction in cross section appended for each iteration

	double density;
	double densitySum;
	double ux;
	double uxSum = 0.0;
	int nFluidNodes = 0;

	for (int y = 0; y < nY; ++y) {
		for (int z = 0; z < nZ; ++z) {

			if (ld->Lattice[L_INDEX_4(ld->Dims, x, y, z)] != LAT_CELL_OBSTACLE) {
				kd->GetNode(kd, x, y, z, pdfs);

				ux = pdfs[D3Q19_E] + pdfs[D3Q19_NE] + pdfs[D3Q19_SE] + pdfs[D3Q19_TE] + pdfs[D3Q19_BE] -
					 pdfs[D3Q19_W] - pdfs[D3Q19_NW] - pdfs[D3Q19_SW] - pdfs[D3Q19_TW] - pdfs[D3Q19_BW];

				uxSum += ux;
				++nFluidNodes;
			}
		}
	}

	const char * mode = "w";

	if (iteration > 0) {
		mode = "a";
	}

	const char * fileName = "ux-progress.dat";
	FILE * fh;

	fh = fopen(fileName, mode);

	if(fh == NULL) {
		printf("ERROR: opening file %s failed.\n", fileName);
		exit(1);
	}

	if (iteration == 0) {
		fprintf(fh, "# Average velocity in x direction of cross section in the middle (x = %d) of the geometry (NX = %d).\n", x, nX);
		fprintf(fh, "# Plot on terminal: gnuplot -e \"set terminal dumb; plot \\\"%s\\\";\"\n", fileName);
		fprintf(fh, "# iteration, avg ux\n");
	}

	fprintf(fh, "%d %e\n", iteration, uxSum / nFluidNodes);

	fclose(fh);

	// ----------------------------------------------------------------------
	// average velocity/density for each in cross section in x direction

	fileName = "density-ux.dat";

	fh = fopen(fileName, "w");

	if(fh == NULL) {
		printf("ERROR: opening file %s failed.\n", fileName);
		exit(1);
	}

	fprintf(fh, "# Average density and average x velocity over each cross section in x direction. Snapshot taken at iteration %d.\n", iteration);
	fprintf(fh, "# Plot on terminal: gnuplot -e \"set terminal dumb; plot \\\"%s\\\" u 1:2; plot \\\"%s\\\" u 1:3;\"\n", fileName, fileName);
//	fprintf(fh, "# Plot graphically: gnuplot -e \"plot \\\"%s\\\" u 1:3 w linesp t \\\"l\\\", \\\"\\\" u 1:4 w linesp t \\\"simulation\\\"; pause -1;"
	fprintf(fh, "# x, avg density, avg ux\n");

	for (x = 0; x < nX; ++x) {

		uxSum = 0.0;
		densitySum = 0.0;
		nFluidNodes = 0;

		for (int y = 0; y < nY; ++y) {
			for (int z = 0; z < nZ; ++z) {

				if (ld->Lattice[L_INDEX_4(ld->Dims, x, y, z)] == LAT_CELL_OBSTACLE) {
					continue;
				}

				kd->GetNode(kd, x, y, z, pdfs);

				density =
					pdfs[D3Q19_C] +
					pdfs[D3Q19_N]  + pdfs[D3Q19_E]  + pdfs[D3Q19_S]  + pdfs[D3Q19_W]  +
					pdfs[D3Q19_NE] + pdfs[D3Q19_SE] + pdfs[D3Q19_SW] + pdfs[D3Q19_NW] +
					pdfs[D3Q19_T]  + pdfs[D3Q19_TN] + pdfs[D3Q19_TE] + pdfs[D3Q19_TS] + pdfs[D3Q19_TW] +
					pdfs[D3Q19_B]  + pdfs[D3Q19_BN] + pdfs[D3Q19_BE] + pdfs[D3Q19_BS] + pdfs[D3Q19_BW];

				densitySum += density;

				ux =
					pdfs[D3Q19_E] + pdfs[D3Q19_NE] + pdfs[D3Q19_SE] + pdfs[D3Q19_TE] + pdfs[D3Q19_BE] -
					pdfs[D3Q19_W] - pdfs[D3Q19_NW] - pdfs[D3Q19_SW] - pdfs[D3Q19_TW] - pdfs[D3Q19_BW];

				uxSum += ux;

				++nFluidNodes;
			}
		}

		fprintf(fh, "%d  %e  %e\n", x, densitySum / nFluidNodes, uxSum / nFluidNodes);
	}

	fclose(fh);
}



void KernelAddBodyForce(KernelData * kd, LatticeDesc * ld, CaseData * cd)
{
	Assert(kd != NULL);
	Assert(ld != NULL);
	Assert(cd != NULL);

	int nX = kd->Dims[0];
	int nY = kd->Dims[1];
	int nZ = kd->Dims[2];

	PdfT w_0 = 1.0 /  3.0; // C
	PdfT w_1 = 1.0 / 18.0; // N,S,E,W,T,B
	PdfT w_2 = 1.0 / 36.0; // NE,NW,SE,SW,TE,TW,BE,BW,TN,TS,BN,BS
	PdfT w[] = {w_1,w_1,w_1,w_1,w_2,w_2,w_2,w_2,w_1,w_2,w_2,w_2,w_2,w_1,w_2,w_2,w_2,w_2,w_0};

	PdfT xForce = cd->XForce;

	PdfT pdfs[N_D3Q19];


	#ifdef _OPENMP
	#pragma omp parallel for collapse(3) default(none) \
			shared(nX,nY,nZ,ld,kd,w,xForce,D3Q19_X,cd) \
			private(pdfs)
	#endif
	for(int z = 0; z < nZ; ++z) {
		for(int y = 0; y < nY; ++y) {
			for(int x = 0; x < nX; ++x) {
 				if(ld->Lattice[L_INDEX_4(ld->Dims, x, y, z)] == LAT_CELL_OBSTACLE)
 					continue;

				// load pdfs into temp array
				kd->GetNode(kd, x, y, z, pdfs);

				// add body force in x direction ( method by Luo)
				for (int d = 0; d < N_D3Q19; ++d) {
					pdfs[d] = pdfs[d] + 3.0*w[d]*D3Q19_X[d]*xForce;
				}

				kd->SetNode(kd, x, y, z, pdfs);

			}
		}
	}
}
