/* ---------------- */
/* --- vnrdef.h --- */
/* ---------------- */

/*
 * Copyright (c) 2004-2017, Lionel Lacassagne
*/

#ifndef __V_NR_DEF_H__
#define __V_NR_DEF_H__

#pragma message("  include  vnrdef.h")

#ifdef __cplusplus
#pragma message ("C++")
extern "C" {
#endif

#include <xmmintrin.h>
#include <emmintrin.h>
#include <mmintrin.h>

typedef __m128i vuint8;
typedef __m128i vsint8;

typedef __m128i vuint16;
typedef __m128i vsint16;

typedef __m128i vuint32;
typedef __m128i vsint32;

typedef __m128 vfloat32;

typedef __m128i vbool;    
    
typedef __m128i  vchar;
typedef __m128i vschar;
typedef __m128i vuchar;

typedef __m128i  vshort;
typedef __m128i vsshort;
typedef __m128i vushort;

typedef __m128i  vint;
typedef __m128i vsint;
typedef __m128i vuint;

typedef __m128 vfloat;

#ifdef __cplusplus
}
#endif

#endif // __V_NR_DEF_H__
