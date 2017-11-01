#pragma once
#include "nrdef.h"

void routine_FrameDifference_SSE2(vuint8 **It, vuint8 **Itm1, vuint8 **Et, long nrl,long nrh,long ncl,long nch, vuint8 seuil);
