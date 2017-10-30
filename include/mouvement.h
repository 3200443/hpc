#pragma once
#include "nrdef.h"

uint8 ** routine_FrameDifference(uint8 **I0, uint8 **I1, uint8 **E0, long nrl,long nrh,long ncl,long nch, int seuil);

void routine_SigmaDelta_step0(uint8** I, uint8 **M, uint8 **V, long nrl, long nrh, long ncl, long nch);

uint8** routine_SigmaDelta_1step(uint8 **I0, uint8 **I1, uint8**V0, uint8 **V1, uint8**M0, uint8 **M1, uint8 **E0,  long nrl, long nrh, long ncl, long nch );
