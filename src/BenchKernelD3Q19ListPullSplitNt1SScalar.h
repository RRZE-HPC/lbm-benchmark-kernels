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

#ifndef INDEX_START
	#error INDEX_START must be defined
#endif

#ifndef INDEX_STOP
	#error INDEX_STOP must be defined
#endif

	#define I(index, dir)	P_INDEX_3((nCells), (index), (dir))

	for (int blockedIndex = (INDEX_START); blockedIndex < (INDEX_STOP); blockedIndex += nTmpArray) {

		indexMax = MinI(nTmpArray, (INDEX_STOP) - blockedIndex);
#ifdef DEBUG
		memset(tmpArray, -1, sizeof(PdfT) * nTmpArray * N_TMP);
#endif
		for (int index = 0; index < indexMax; ++index) {


			adjListIndex = (index + blockedIndex) * N_D3Q19_IDX;

			pdf_C = src[I(index + blockedIndex, D3Q19_C)];

			#define X(name, idx, idxinv, _x, _y, _z)	JOIN(pdf_,name) = src[adjList[adjListIndex + idx]]; tmpArray[TMP_INDEX(index, idx)] = JOIN(pdf_,name);
			D3Q19_LIST_WO_C
			#undef X

			ux = pdf_E + pdf_NE + pdf_SE + pdf_TE + pdf_BE -
				 pdf_W - pdf_NW - pdf_SW - pdf_TW - pdf_BW;
			uy = pdf_N + pdf_NE + pdf_NW + pdf_TN + pdf_BN -
				 pdf_S - pdf_SE - pdf_SW - pdf_TS - pdf_BS;
			uz = pdf_T + pdf_TE + pdf_TW + pdf_TN + pdf_TS -
				 pdf_B - pdf_BE - pdf_BW - pdf_BN - pdf_BS;

			tmpArray[TMP_INDEX(index, TMP_UX)] = ux;
			tmpArray[TMP_INDEX(index, TMP_UY)] = uy;
			tmpArray[TMP_INDEX(index, TMP_UZ)] = uz;

			dens = pdf_C +
				   pdf_N  + pdf_E  + pdf_S  + pdf_W  +
				   pdf_NE + pdf_SE + pdf_SW + pdf_NW +
				   pdf_T  + pdf_TN + pdf_TE + pdf_TS + pdf_TW +
				   pdf_B  + pdf_BN + pdf_BE + pdf_BS + pdf_BW;

			dir_indep_trm = dens - (ux * ux + uy * uy + uz * uz) * F(3.0) / F(2.0);

			w_1_indep = w_1 * dir_indep_trm;
			w_2_indep = w_2 * dir_indep_trm;

			tmpArray[TMP_INDEX(index, TMP_W1)] = w_1_indep;
			tmpArray[TMP_INDEX(index, TMP_W2)] = w_2_indep;

			dst[I(index + blockedIndex, D3Q19_C )]  = pdf_C - omegaEven * (pdf_C - w_0 * dir_indep_trm);
		}


		#define LOOP_1(_dir1, _dir2, _vel, _vel_tmp) \
			for (int index = 0; index < indexMax; ++index) { \
				_vel               = tmpArray[TMP_INDEX(index, JOIN(TMP_,_vel_tmp))]; \
				JOIN(pdf_,_dir1)   = tmpArray[TMP_INDEX(index, JOIN(D3Q19_,_dir1))]; \
				JOIN(pdf_,_dir2)   = tmpArray[TMP_INDEX(index, JOIN(D3Q19_,_dir2))]; \
				w_1_indep          = tmpArray[TMP_INDEX(index, TMP_W1)]; \
				\
				ui = _vel; \
				evenPart = omegaEven * (F(0.5) * (JOIN(pdf_,_dir1) + JOIN(pdf_,_dir2)) - ui * ui * w_1_nine_half - w_1_indep); \
				oddPart  = omegaOdd  * (F(0.5) * (JOIN(pdf_,_dir1) - JOIN(pdf_,_dir2)) - ui * w_1_x3); \
				dst[I(index + blockedIndex, JOIN(D3Q19_,_dir1) )]  = JOIN(pdf_,_dir1) - evenPart - oddPart; \
				tmpArray[TMP_INDEX(index, JOIN(D3Q19_,_dir2))]     = JOIN(pdf_,_dir2) - evenPart + oddPart; \
			} \
			for (int index = 0; index < indexMax; ++index) { \
				dst[I(index + blockedIndex, JOIN(D3Q19_,_dir2) )] = tmpArray[TMP_INDEX(index, JOIN(D3Q19_,_dir2))]; \
			}

		#define LOOP_2(_dir1, _dir2, _v1, _v2, _v1_tmp, _v2_tmp, _expr) \
			for (int index = 0; index < indexMax; ++index) { \
				_v1 = tmpArray[TMP_INDEX(index, JOIN(TMP_,_v1_tmp))]; \
				_v2 = tmpArray[TMP_INDEX(index, JOIN(TMP_,_v2_tmp))]; \
				JOIN(pdf_,_dir1) = tmpArray[TMP_INDEX(index, JOIN(D3Q19_,_dir1))]; \
				JOIN(pdf_,_dir2) = tmpArray[TMP_INDEX(index, JOIN(D3Q19_,_dir2))]; \
				w_2_indep = tmpArray[TMP_INDEX(index, TMP_W2)]; \
				\
				ui = _expr; \
				evenPart = omegaEven * (F(0.5) * (JOIN(pdf_,_dir1) + JOIN(pdf_,_dir2)) - ui * ui * w_2_nine_half - w_2_indep); \
				oddPart  = omegaOdd  * (F(0.5) * (JOIN(pdf_,_dir1) - JOIN(pdf_,_dir2)) - ui * w_2_x3); \
				dst[I(index + blockedIndex, JOIN(D3Q19_,_dir1))] = JOIN(pdf_,_dir1) - evenPart - oddPart; \
				tmpArray[TMP_INDEX(index, JOIN(D3Q19_,_dir2))]     = JOIN(pdf_,_dir2) - evenPart + oddPart; \
			} \
			for (int index = 0; index < indexMax; ++index) { \
				dst[I(index + blockedIndex, JOIN(D3Q19_,_dir2) )] = tmpArray[TMP_INDEX(index, JOIN(D3Q19_,_dir2))]; \
			}

		LOOP_1(N, S, uy, UY);
		LOOP_1(E, W, ux, UX);
		LOOP_1(T, B, uz, UZ);

		LOOP_2(NW, SE, uy, ux, UY, UX, uy - ux);
		LOOP_2(NE, SW, uy, ux, UY, UX, uy + ux);
		LOOP_2(TW, BE, ux, uz, UX, UZ, uz - ux);
		LOOP_2(TE, BW, ux, uz, UX, UZ, uz + ux);
		LOOP_2(TS, BN, uy, uz, UY, UZ, uz - uy);
		LOOP_2(TN, BS, uy, uz, UY, UZ, uz + uy);

		#undef LOOP_1
		#undef LOOP_2


	}

	#undef I

	#undef INDEX_START
	#undef INDEX_STOP
