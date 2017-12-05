#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "nrdef.h"
#include "nrutil.h"

#define N 2
#define VMIN 20
#define VMAX 255 //V est entre 2 et 2^m-1 avec m le nombre de bits des donnees ici 8 => https://hal.inria.fr/hal-01130889/document

__attribute__ ((always_inline))int inline minO(int a, int b)
{
	return a<b ? a:b;
}
int min(int a, int b)
{
	return a<b ? a:b;
}


__attribute__ ((always_inline))int inline maxO(int a, int b)
{
	return a>b ? a:b;
}

int max(int a, int b)
{
	return a>b ? a:b;
}

void routine_FrameDifference(uint8 **It, uint8 **Itm1, uint8 **Et, long nrl,long nrh,long ncl,long nch, int seuil)
{
    //m[nrl..nrh][ncl..nch]

	uint8 Ot;
	for(int i = nrl; i <= nrh; i++ )
	{
		for(int j = ncl; j <= nch; j++)
		{
			Ot = abs(It[i][j] - Itm1[i][j]);
			if(Ot < seuil)
				Et[i][j] = 0;
			else
				Et[i][j] = 255;
		}
	}
}

void routine_SigmaDelta_step0(uint8** I, uint8 **M, uint8 **V, long nrl, long nrh, long ncl, long nch)
{
	for(int i = nrl; i <= nrh; i++ )
	{
		for(int j = ncl; j <= nch; j++)
		{
			M[i][j] = I[i][j];
            V[i][j] = VMIN; 
        }
    }
}

void routine_SigmaDelta_1step(uint8 **It, uint8**Vt, uint8 **Vtm1, uint8**Mt, uint8 **Mtm1, uint8 **Et,  long nrl, long nrh, long ncl, long nch )
{
	uint8 Ot;
	uint8 tmpMtm1, tmpVtm1, tmpIt, tmpMt, tmpVt;
    for(int i = nrl; i <= nrh; i++ ) //Step1 Mt Estimation
    {
    	for(int j = ncl; j <= nch; j++)
    	{
    		tmpMtm1 = Mtm1[i][j];
    		tmpVtm1 = Vtm1[i][j];
    		tmpIt = It[i][j];

    		if(tmpMtm1 < tmpIt)
    			tmpMt  = tmpMtm1 + 1;

    		else if(tmpMtm1 > tmpIt)
    			tmpMt = tmpMtm1 - 1;

    		else
    			tmpMt = tmpMtm1;


    		//Step 2 difference Computation
    		Ot = abs(tmpMt - tmpIt);


    		//Step 3 Update and clamping
    		if(tmpVtm1 < N * Ot)
    			tmpVt = tmpVtm1 + 1;

    		else if(tmpVtm1 > N * Ot)
    			tmpVt = tmpVtm1 - 1;

    		else
    			tmpVt = tmpVtm1;

    		tmpVt = max( min(tmpVt, VMAX), VMIN);


    		//Step 4 Et estimation
    		if(Ot < tmpVt)
    			Et[i][j] = 0;
    		else
    			Et[i][j] = 255;

    		Mt[i][j] = tmpMt;
    		Vt[i][j] = tmpVt;
    	}
    }

}

void routine_SigmaDelta_1stepO(uint8 ** restrict It, uint8 ** restrict Itm1,uint8** restrict Vt, uint8 ** restrict Vtm1, 
	uint8** restrict Mt, uint8 ** restrict Mtm1,uint8 ** restrict Et,  long nrl, long nrh, long ncl, long nch )
{
	uint8 **Ot = ui8matrix(nrl, nrh, ncl, nch);
    for(int i = nrl; i <= nrh; i++ ) //Step1 Mt Estimation
    {
    	for(int j = ncl; j <= nch; j++)
    	{
    		if(Mtm1[i][j] < It[i][j])
    			Mt[i][j]  = Mtm1[i][j] + 1;

    		else if(Mtm1[i][j] > It[i][j])
    			Mt[i][j] = Mtm1[i][j] - 1;

    		else
    			Mt[i][j] = Mtm1[i][j];


    	}
    }

    for(int i = nrl; i <= nrh; i++)//Step 2 difference Computation
    {
    	for(int j = ncl; j <= nch; j++)
    	{
    		Ot[i][j] = abs(Mt[i][j] - It[i][j]);
    	}
    }


    for(int i = nrl; i <= nrh; i++)//Step 3 Update and clamping
    {
    	for(int j = ncl; j <= nch; j++)
    	{

    		if(Vtm1[i][j] < N * Ot[i][j])
    			Vt[i][j] = Vtm1[i][j] + 1;

    		else if(Vtm1[i][j] > N * Ot[i][j])
    			Vt[i][j] = Vtm1[i][j] - 1;

    		else
    			Vt[i][j] = Vtm1[i][j];


    		Vt[i][j] = maxO( minO(Vt[i][j], VMAX), VMIN);


    	}
    }

    for(int i = nrl; i <= nrh; i++)//Step 4 Et estimation
    {
    	for(int j = ncl; j <= nch; j++)
    	{
    		if(Ot[i][j] < Vt[i][j])
    			Et[i][j] = 0;
    		else
    			Et[i][j] = 255;
    	}
    }
    free_ui8matrix(Ot, nrl, nrh, ncl, nch);

}