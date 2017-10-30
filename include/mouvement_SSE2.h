#pragma once
#include "nrdef.h"

uint8 ** routine_FrameDifference_SSE2(uint8 **I0, uint8 **I1, uint8 **E0, long nrl,long nrh,long ncl,long nch, int seuil);
