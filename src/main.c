#include <stdio.h>
#include <stdlib.h>
#include "test_mouvement.h"
#include "test_mouvement_SSE2.h"
#define NBIMAGES 199
#include "nrdef.h"
#include "nrutil.h"



void creation_ppm(){
	char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    rgb8 **It = rgb8matrix(nrl, nrh, ncl, nch);
	for(int i = 0; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, Itm1);

        for(int k = nrl; k <= nrh; k++){
            for(int l = ncl; l <= nch; l++){
                It[k][l].r = Itm1[k][l];
                It[k][l].g = Itm1[k][l];
                It[k][l].b = Itm1[k][l];
            }
        }

        sprintf(nomImageSave, "/home/souley/Bureau/ppm/car_3%03d.ppm", i);
        SavePPM_rgb8matrix(It, nrl, nrh, ncl, nch, nomImageSave);
    }

}

#define OPTI 3 //1 pour optimisation 2 sans optimisation 3 pour tout 0 pour rien


int main(int argc, char* argv[])
{
#if OPTI & 0x1
    test_routine_FrameDifference_SSE2(10);
    test_routine_sigmaDelta_SSE2();
#endif
#if OPTI & 0x2
    test_routine_FrameDifference(10);
    test_routine_FrameDifferenceMorpho3x3ouverture(10);
    test_routine_FrameDifferenceMorpho3x3fermeture(10);
    test_routine_FrameDifferenceMorpho3x3ouvertureFermeture(10);
    test_routine_FrameDifferenceMorpho3x3fermetureOuverture(10);
    test_routine_sigmaDelta();
#endif

    //creation_ppm();
    return 0;
}
