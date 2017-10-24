#pragma once
#include "nrdef.h"

void dilatation3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
void erosion3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);