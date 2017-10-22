#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "nrdef.h"
#include "nrutil.h"

void erosion3x3(uint8** entree,uint8 sortie,uint32 taille)
{
	int i = 0;
	int j = 0;

	for(i; i < taille; i++)
	{
		uint8 x1 = entree[i-1][j-1];
		uint8 x2 = entree[i-1][j];
		
		for(int j = 0; j<taille; j++)
		{

		}
	}
}