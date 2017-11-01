#pragma once
#include "nrdef.h"
#include "vnrdef.h"
void init_tab(vuint8 **vX1, uint8 **Itm1, int nrl, int nrh, int ncl, int nch);
void init_tab_inv(vuint8 **vX1, uint8 **Itm1, int nrl, int nrh, int ncl, int nch);
void test_routine_FrameDifference_SSE2(int seuil);
