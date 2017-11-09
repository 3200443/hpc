#pragma once
#include "nrdef.h"
#include "vnrdef.h"
void test_unitaire_SD_SSE2();
void MatScal2MatSIMD(vuint8 **vX1, uint8 **Itm1, int nrl, int nrh, int ncl, int nch);
void MatSIMD2MatScal(vuint8 **vX1, uint8 **Itm1, int nrl, int nrh, int ncl, int nch);
void test_routine_sigmaDelta_SSE2();
void test_routine_FrameDifference_SSE2(int seuil);
void test_routine_FrameDifference_SSE2M(int seuil);