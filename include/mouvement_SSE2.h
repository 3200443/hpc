#pragma once
#include "nrdef.h"

void routine_FrameDifference_SSE2(vuint8 **It, vuint8 **Itm1, vuint8 **Et, long nrl,long nrh,long ncl,long nch, vuint8 seuil);
void routine_SigmaDelta_step0SSE2(vuint8** I, vuint8 **M, vuint8 **V, long vi0, long vi1, long vj0, long vj1);
void routine_SigmaDelta_1stepSSE2(vuint8 **It, vuint8 **Itm1, vuint8**Vt, vuint8 **Vtm1, vuint8**Mt, vuint8 **Mtm1, vuint8 **Et,  long vi0, long vi1, long vj0, long vj1 );
