#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "nrdef.h"
#include "nrutil.h"
#define BORD 2

void erosion3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
	int i,j,k,l;
	uint8 result;
	for(i=nrl;i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
			result = 255;
			for(k = i-1 ; k < i+2 ; k++)
			{
				for(l = j-1;l<j+2 ; l++)
				{
					result &= X[k][l];
				}
			}
			Y[i][j] = result;
		}
	}
}

void dilatation3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
	int i,j,k,l;
	uint8 result;
	for(i=nrl;i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
			result = 0;
			for(k = i-1 ; k < i+2 ; k++)
			{
				for(l = j-1;l<j+2 ; l++)
				{
					result |= X[k][l];
				}
			}
			Y[i][j] = result;
		}
	}
}

void fermeture3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
	uint8 **O0 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
	dilatation3x3(X, O0, nrl, nrh, ncl, nch);
	erosion3x3(O0, Y,nrl,nrh,ncl,nch);
	free_ui8matrix(O0, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}

void ouverture3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
	uint8 **O0 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
	erosion3x3(X, O0,nrl,nrh,ncl,nch);
	dilatation3x3(O0, Y, nrl, nrh, ncl, nch);
	free_ui8matrix(O0, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}

void erosion5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
	int i,j,k,l;
	uint8 result;
	for(i=nrl;i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
			result = 255;
			for(k = i-2 ; k < i+3 ; k++)
			{
				for(l = j-2;l<j+3 ; l++)
				{
					result &= X[k][l];
				}
			}
			Y[i][j] = result;
		}
	}
}

void dilatation5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
	int i,j,k,l;
	uint8 result;
	for(i=nrl;i<=nrh;i++)
	{
		for(j=ncl;j<=nch;j++)
		{
			result = 0;
			for(k = i-2 ; k < i+3 ; k++)
			{
				for(l = j-2;l<j+3 ; l++)
				{
					result |= X[k][l];
				}
			}
			Y[i][j] = result;
		}
	}
}

void fermeture5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
	uint8 **O0 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
	dilatation5x5(X, O0, nrl, nrh, ncl, nch);
	erosion5x5(O0, Y,nrl,nrh,ncl,nch);
	free_ui8matrix(O0, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}

void ouverture5x5(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
	uint8 **O0 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
	erosion5x5(X, O0,nrl,nrh,ncl,nch);
	dilatation5x5(O0, Y, nrl, nrh, ncl, nch);
	free_ui8matrix(O0, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
}
