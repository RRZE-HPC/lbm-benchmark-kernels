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
#ifndef __VECTOR_H__
#define __VECTOR_H__

#if !defined(VECTOR_AVX) && !defined(VECTOR_SSE)
	#warning Defining VECTOR_AVX as no ISA extension was selected.
	#define VECTOR_AVX
#endif

#if defined(VECTOR_AVX) && defined(VECTOR_SSE)
	#error Only VECTOR_AVX or VECTOR_SSE can be defined at the same time.
#endif

#ifdef VECTOR_AVX

	#include <immintrin.h>
	// Vector size in double-precision floatin-point numbers.
	#define VSIZE	4

	#define VPDFT				__m256d

	#define VSET(scalar)		_mm256_set1_pd(scalar)

	#define VLD(expr)			_mm256_load_pd(expr)
	#define VLDU(expr)			_mm256_loadu_pd(expr)

	#define VST(dst, src)		_mm256_store_pd(dst, src)
	#define VSTU(dst, src)		_mm256_storeu_pd(dst, src)
	#define VSTNT(dst, src)		_mm256_stream_pd(dst, src)

	#define VMUL(a, b)			_mm256_mul_pd(a, b)
	#define VADD(a, b)			_mm256_add_pd(a, b)
	#define VSUB(a, b)			_mm256_sub_pd(a, b)
#endif

#ifdef VECTOR_SSE
	#include <emmintrin.h>
	// Vector size in double-precision floatin-point numbers.
	#define VSIZE 2

	#define VPDFT				__m128d

	#define VSET(scalar)		_mm_set1_pd(scalar)

	#define VLD(expr)			_mm_load_pd(expr)
	#define VLDU(expr)			_mm_loadu_pd(expr)

	#define VST(dst, src)		_mm_store_pd(dst, src)
	#define VSTU(dst, src)		_mm_storeu_pd(dst, src)
	#define VSTNT(dst, src)		_mm_stream_pd(dst, src)

	#define VMUL(a, b)			_mm_mul_pd(a, b)
	#define VADD(a, b)			_mm_add_pd(a, b)
	#define VSUB(a, b)			_mm_sub_pd(a, b)
#endif


#endif // __VECTOR_H__
