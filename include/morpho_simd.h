#pragma once
#include "nrdef.h"
#include "vnrdef.h"

void erosion3x3_SIMD(vuint8 **It0,vuint8 **It1,long vi0,long vi1,long vj0,long vj1);
void dilatation3x3_SIMD(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1);
void fermeture3x3_SIMD(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1);
void ouverture3x3_SIMD(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1);
void erosion5x5_SIMD(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1);
void dilatation5x5_SIMD(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1);
void fermeture5x5_SIMD(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1);
void ouverture5x5_SIMD(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1);
void erosion3x3_SIMD_B(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1);
void dilatation3x3_SIMD_B(vuint8 **It,vuint8 **It1,long vi0,long vi1,long vj0,long vj1);