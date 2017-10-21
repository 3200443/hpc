#include <stdio.h>
#include <stdlib.h>
#include "mouvement.h"





int main(int argc, char* argv[]){
	long nrl, nrh, ncl, nch;
	uint8 **m1 =  LoadPGM_ui8matrix("/home/souley/HPC/Projet/car3/car_3063.pgm", &nrl, &nrh, &ncl, &nch);
	uint8 **m2 =  LoadPGM_ui8matrix("/home/souley/HPC/Projet/car3/car_3064.pgm", &nrl, &nrh, &ncl, &nch);


	routine_FrameDifference(m1,m2, atoi(argv[1]));




	return 0;
}