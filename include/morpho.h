#pragma once
#include "nrdef.h"

void dilatation3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
void erosion3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
void fermeture3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
void ouverture3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);

void dilatation5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
void erosion5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
void fermeture5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
void ouverture5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
