
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h> 
#include <string.h>

#include "nrdef.h"
#include "nrutil.h"



void routine_FrameDifference(char* nomFichier1, char* nomFichier2, int seuil){
	//m[nrl..nrh][ncl..nch]
	long nrl, nrh, ncl, nch;
	uint8 **m1 =  LoadPGM_ui8matrix(nomFichier1, &nrl, &nrh, &ncl, &nch);
	uint8 **m2 =  LoadPGM_ui8matrix(nomFichier2, &nrl, &nrh, &ncl, &nch);


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
	SavePGM_ui8matrix(EtiqBin, nrl, nrh, ncl, nch, "test.pgm");
}

