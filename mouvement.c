
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h> 
#include <string.h>

#include "nrdef.h"
#include "nrutil.h"



uint8 ** routine_FrameDifference(uint8 **m1, uint8 **m2,long nrl,long nrh,long ncl,long nch, int seuil){
	//m[nrl..nrh][ncl..nch]


	uint8 **EtiqBin =  ui8matrix(nrl, nrh, ncl, nch);
	uint8 **resm = ui8matrix(nrl, nrh, ncl, nch);
	for(int i = nrl; i < nrh; i++ ){
		for(int j = ncl; j < nch; j++){
			m2[i][j] = abs(m2[i][j] - m1[i][j]);

		}
	}
	for(int i = nrl; i < nrh; i++ ){
		for(int j = ncl; j < nch; j++){
			if(m2[i][j] < seuil)
				EtiqBin[i][j] = 0;
			else
				EtiqBin[i][j] = 255;

		}
	}
	return EtiqBin;
}

