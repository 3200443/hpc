#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "nrdef.h"
#include "nrutil.h"

void erosion3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{

	int i,j;
	uint8 result,x1,x2,x3,x4,x5,x6,x7,x8,x9;


	for(i=nrl;i<nrh;i++)
	{
		//Prologue:
		x1 = X[i-1][ncl-1];
		x2 = X[i][ncl-1];
		x3 = X[i+1][ncl-1];

		x4 = X[i-1][ncl];
		x5 = X[i][ncl];
		x6 = X[i+1][ncl];
		for(j=ncl;j<nch-2;j+=3)
		{
			x7 = X[i-1][j+1];
			x8 = X[i][j+1];
			x9 = X[i+1][j+1];
			result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
			Y[i][j] = result;

			x1 = X[i-1][j+2];
			x2 = X[i][j+2];
			x3 = X[i+1][j+2];
			result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
			Y[i][j+1] = result;

			x4 = X[i-1][j+3];
			x5 = X[i][j+3];
			x6 = X[i+1][j+3];
			result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
			Y[i][j+2] = result;
		}
		//Epilogue
		switch((nch-2)%3)
		{
			case 0 :
			{
				x7 = X[i-1][j+1];
				x8 = X[i][j+1];
				x9 = X[i+1][j+1];
				result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
				Y[i][j] = result;

				x1 = X[i-1][j+2];
				x2 = X[i][j+2];
				x3 = X[i+1][j+2];
				result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
				Y[i][j+1] = result;
				break;
			}
			case 1:
			{
				if(j>(nch-2))
				{
					x7 = X[i-1][j+1];
					x8 = X[i][j+1];
					x9 = X[i+1][j+1];
					result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
					Y[i][j] = result;
					break;
				}else
				{
					x7 = X[i-1][j+1];
					x8 = X[i][j+1];
					x9 = X[i+1][j+1];
					result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
					Y[i][j] = result;

					x1 = X[i-1][j+2];
					x2 = X[i][j+2];
					x3 = X[i+1][j+2];
					result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
					Y[i][j+1] = result;

					x4 = X[i-1][j+3];
					x5 = X[i][j+3];
					x6 = X[i+1][j+3];
					result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
					Y[i][j+2] = result;

					x4 = X[i-1][j+4];
					x5 = X[i][j+4];
					x6 = X[i+1][j+4];
					result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
					Y[i][j+3] = result;

					break;
				}
			}
			case 2:
			{
				if(j>(nch-2))
				{
					x7 = X[i-1][j+1];
					x8 = X[i][j+1];
					x9 = X[i+1][j+1];
					result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
					Y[i][j] = result;

					x1 = X[i-1][j+2];
					x2 = X[i][j+2];
					x3 = X[i+1][j+2];
					result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
					Y[i][j+1] = result;
					break;
				}else
				{
					x7 = X[i-1][j+1];
					x8 = X[i][j+1];
					x9 = X[i+1][j+1];
					result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
					Y[i][j] = result;

					x1 = X[i-1][j+2];
					x2 = X[i][j+2];
					x3 = X[i+1][j+2];
					result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
					Y[i][j+1] = result;

					x4 = X[i-1][j+3];
					x5 = X[i][j+3];
					x6 = X[i+1][j+3];
					result = x1&x2&x3&x4&x5&x6&x7&x8&x9;
					Y[i][j+2] = result;
				}
			}
		}
	}
}

void dilatation3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{

	int i,j;
	uint8 result,x1,x2,x3,x4,x5,x6,x7,x8,x9;


	for(i=nrl;i<nrh;i++)
	{
		//Prologue:
		x1 = X[i-1][ncl-1];
		x2 = X[i][ncl-1];
		x3 = X[i+1][ncl-1];

		x4 = X[i-1][ncl];
		x5 = X[i][ncl];
		x6 = X[i+1][ncl];
		for(j=ncl;j<nch-2;j+=3)
		{
			x7 = X[i-1][j+1];
			x8 = X[i][j+1];
			x9 = X[i+1][j+1];
			result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
			Y[i][j] = result;

			x1 = X[i-1][j+2];
			x2 = X[i][j+2];
			x3 = X[i+1][j+2];
			result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
			Y[i][j+1] = result;

			x4 = X[i-1][j+3];
			x5 = X[i][j+3];
			x6 = X[i+1][j+3];
			result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
			Y[i][j+2] = result;
		}
		//Epilogue
		switch((nch-2)%3)
		{
			case 0 :
			{
				x7 = X[i-1][j+1];
				x8 = X[i][j+1];
				x9 = X[i+1][j+1];
				result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
				Y[i][j] = result;

				x1 = X[i-1][j+2];
				x2 = X[i][j+2];
				x3 = X[i+1][j+2];
				result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
				Y[i][j+1] = result;
				break;
			}
			case 1:
			{
				if(j>(nch-2))
				{
					x7 = X[i-1][j+1];
					x8 = X[i][j+1];
					x9 = X[i+1][j+1];
					result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
					Y[i][j] = result;
					break;
				}else
				{
					x7 = X[i-1][j+1];
					x8 = X[i][j+1];
					x9 = X[i+1][j+1];
					result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
					Y[i][j] = result;

					x1 = X[i-1][j+2];
					x2 = X[i][j+2];
					x3 = X[i+1][j+2];
					result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
					Y[i][j+1] = result;

					x4 = X[i-1][j+3];
					x5 = X[i][j+3];
					x6 = X[i+1][j+3];
					result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
					Y[i][j+2] = result;

					x4 = X[i-1][j+4];
					x5 = X[i][j+4];
					x6 = X[i+1][j+4];
					result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
					Y[i][j+3] = result;

					break;
				}
			}
			case 2:
			{
				if(j>(nch-2))
				{
					x7 = X[i-1][j+1];
					x8 = X[i][j+1];
					x9 = X[i+1][j+1];
					result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
					Y[i][j] = result;

					x1 = X[i-1][j+2];
					x2 = X[i][j+2];
					x3 = X[i+1][j+2];
					result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
					Y[i][j+1] = result;
					break;
				}else
				{
					x7 = X[i-1][j+1];
					x8 = X[i][j+1];
					x9 = X[i+1][j+1];
					result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
					Y[i][j] = result;

					x1 = X[i-1][j+2];
					x2 = X[i][j+2];
					x3 = X[i+1][j+2];
					result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
					Y[i][j+1] = result;

					x4 = X[i-1][j+3];
					x5 = X[i][j+3];
					x6 = X[i+1][j+3];
					result = x1|x2|x3|x4|x5|x6|x7|x8|x9;
					Y[i][j+2] = result;
				}
			}
		}
	}
}

void ouverture3x3(uint8** X,uint8** Y, long nrl,long nrh,long ncl,long nch)
{
	
}