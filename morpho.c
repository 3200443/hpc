#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "nrdef.h"
#include "nrutil.h"

void erosion3x3(uint8** X,uint8** Y,uint32 taille)
{

	int i,j;
	uint8 result,x1,x2,x3,x4,x5,x6,x7,x8,x9;


	for(i=0;i<taille;i++)
	{
		//Prologue:
		x1 = X[i-1][-1];
		x2 = X[i][-1];
		x3 = X[i+1][-1];

		x4 = X[i-1][0];
		x5 = X[i][0];
		x6 = X[i+1][0];
		for(j=0;j<taille-2;j+=3)
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
		switch((taille-2)%3)
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
				if(j>(taille-2))
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
				if(j>(taille-2))
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

void dilatation3x3(uint8** X,uint8** Y,uint32 taille)
{

	int i,j;
	uint8 result,x1,x2,x3,x4,x5,x6,x7,x8,x9;


	for(i=0;i<taille;i++)
	{
		//Prologue:
		x1 = X[i-1][-1];
		x2 = X[i][-1];
		x3 = X[i+1][-1];

		x4 = X[i-1][0];
		x5 = X[i][0];
		x6 = X[i+1][0];
		for(j=0;j<taille-2;j+=3)
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
		switch((taille-2)%3)
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
				if(j>(taille-2))
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
				if(j>(taille-2))
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