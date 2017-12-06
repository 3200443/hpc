#pragma once
#include "nrdef.h"

void dilatation3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
void erosion3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
void fermeture3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
void ouverture3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);

void erosion3x3_line(uint8** X,uint8** Y, int i,long ncl,long nch);
void dilatation3x3_line(uint8** X,uint8** Y, int i,long ncl,long nch);
void fermeture3x3_pipe(uint8** X,uint8** Y,uint8** O0, long nrl,long nrh,long ncl,long nch);
void ouverture3x3_pipe(uint8** X,uint8** Y,uint8** O0, long nrl,long nrh,long ncl,long nch);
void ouverture3x3_bin(ulong32** X,ulong32** Y, ulong32** O0, long bi0,long bi1,long bj0,long bj1);
void fermeture3x3_bin(ulong32** X,ulong32** Y, ulong32** O0, long bi0,long bi1,long bj0,long bj1);


void dilatation5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
void erosion5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
void fermeture5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
void ouverture5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch);
