// --------------------------------------------------------------------------
//
// Copyright
//   Markus Wittmann, 2016-2018
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
#ifndef __VECTOR_H__
#define __VECTOR_H__

#if !defined(VECTOR_AVX512) && !defined(VECTOR_AVX) && !defined(VECTOR_SSE)
	#warning Defining VECTOR_AVX as no ISA extension was selected.
	#define VECTOR_AVX
#endif

#if (defined(VECTOR_AVX512) && defined(VECTOR_AVX)) || \
	(defined(VECTOR_AVX512) && defined(VECTOR_SEE)) || \
	(defined(VECTOR_AVX) && defined(VECTOR_SSE))
	#error Only VECTOR_AVX512 or VECTOR_AVX or VECTOR_SSE can be defined at the same time.
#endif

#if !defined(PRECISION_DP) && !defined(PRECISION_SP)
	#error PRECISION_DP or PRECISION_SP must be defined.
#endif

#if defined(PRECISION_DP) && defined(PRECISION_SP)
	#error Only PRECISION_DP or PRECISION_SP can be defined at the same time.
#endif

#ifdef PRECISION_DP

	#ifdef VECTOR_AVX512

		#include <immintrin.h>
		// Vector size in double-precision floatin-point numbers.
		#define VSIZE	8

		#define VPDFT				__m512d

		#define VSET(scalar)		_mm512_set1_pd(scalar)
		#define VSETI32(scalar)     _mm256_set1_epi32(scalar)

		#define VLD(expr)			_mm512_load_pd(expr)
		#define VLDU(expr)			_mm512_loadu_pd(expr)
		#define VLIU(expr)			_mm256_loadu_si256((__m256i const *)expr)
		#define VLI64(expr)			_mm512_load_epi64(expr)

		#define VST(dst, src)		_mm512_store_pd(dst, src)
		#define VSTU(dst, src)		_mm512_storeu_pd(dst, src)
		#define VSTNT(dst, src)		_mm512_stream_pd(dst, src)

		#define VG32(offsets, base, scale)	_mm512_i32gather_pd(offsets, base, scale)
		#define VG64(offsets, base, scale)	_mm512_i64gather_pd(offsets, base, scale)

		#define VPG32(offsets, base, scale, hint) _mm512_prefetch_i32gather_pd(offsets, base, scale, hint)

		#define VS32(dst_base, dst_offsets, src, scale)	_mm512_i32scatter_pd(dst_base, dst_offsets, src, scale)
		#define VS64(dst_base, dst_offsets, src, scale)	_mm512_i64scatter_pd(dst_base, dst_offsets, src, scale)

		#define VPS32(dst_base, dst_offsets, scale, hint) _mm512_prefetch_i32scatter_pd(dst_base, dst_offsets, scale, hint)

		#define VMUL(a, b)			_mm512_mul_pd(a, b)
		#define VADD(a, b)			_mm512_add_pd(a, b)
		#define VADDI32(a,b)		_mm256_add_epi32(a,b)
		#define VMULI32(a,b)		_mm256_mul_epi32(a,b)
		#define VSUB(a, b)			_mm512_sub_pd(a, b)
	#endif

	#ifdef VECTOR_AVX

		#include <immintrin.h>
		// Vector size in double-precision floating-point numbers.
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
		// Vector size in double-precision floating-point numbers.
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

#elif defined(PRECISION_SP)

	#ifdef VECTOR_AVX512
		#error Single precision intrinsic kernels for AVX512 are currently not implemented.
	#endif

	#ifdef VECTOR_AVX

		#include <immintrin.h>
		// Vector size in double-precision floating-point numbers.
		#define VSIZE	8

		#define VPDFT				__m256

		#define VSET(scalar)		_mm256_set1_ps(scalar)

		#define VLD(expr)			_mm256_load_ps(expr)
		#define VLDU(expr)			_mm256_loadu_ps(expr)

		#define VST(dst, src)		_mm256_store_ps(dst, src)
		#define VSTU(dst, src)		_mm256_storeu_ps(dst, src)
		#define VSTNT(dst, src)		_mm256_stream_ps(dst, src)

		#define VMUL(a, b)			_mm256_mul_ps(a, b)
		#define VADD(a, b)			_mm256_add_ps(a, b)
		#define VSUB(a, b)			_mm256_sub_ps(a, b)
	#endif

	#ifdef VECTOR_SSE
		#include <emmintrin.h>
		// Vector size in double-precision floating-point numbers.
		#define VSIZE	4

		#define VPDFT				__m128

		#define VSET(scalar)		_mm_set1_ps(scalar)

		#define VLD(expr)			_mm_load_ps(expr)
		#define VLDU(expr)			_mm_loadu_ps(expr)

		#define VST(dst, src)		_mm_store_ps(dst, src)
		#define VSTU(dst, src)		_mm_storeu_ps(dst, src)
		#define VSTNT(dst, src)		_mm_stream_ps(dst, src)

		#define VMUL(a, b)			_mm_mul_ps(a, b)
		#define VADD(a, b)			_mm_add_ps(a, b)
		#define VSUB(a, b)			_mm_sub_ps(a, b)
	#endif

#endif  // PRECISION

#endif // __VECTOR_H__
