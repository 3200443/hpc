#include <stdio.h>
#include <stdlib.h>
#include "test_mouvement.h"
#include "test_mouvement_SSE2.h"
#define NBIMAGES 199
#include "nrdef.h"
#include "nrutil.h"
#include "matric_roc.h"



void creation_ppm()
{

    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3Frame/car_3158.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = LoadPGM_ui8matrix("car3FrameSIMD/car_3158.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **Image = LoadPGM_ui8matrix("car3/car_3158.pgm",&nrl, &nrh, &ncl, &nch);
    uint8 **Image2 = LoadPGM_ui8matrix("car3/car_3104.pgm",&nrl, &nrh, &ncl, &nch);

    for(int k = nrl; k <= nrh; k++)
    {
        for(int l = ncl; l <= nch; l++)
        {
            if(Itm1[k][l] != It[k][l]){
                printf("Probleme ! i = %d, j = %d, ItNormal = %d, Itavant = %d\n",k , l, Image[k][l], Image2[k][l]);
            }
        }
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
    // test_routine_FrameDifferenceMorpho3x3ouverture(10);
    // test_routine_FrameDifferenceMorpho3x3fermeture(10);
    // test_routine_FrameDifferenceMorpho3x3ouvertureFermeture(10);
    test_routine_FrameDifferenceMorpho3x3fermetureOuverture(10);
    //test_routine_sigmaDelta();
#endif
    //creation_matrices_ROC("verite/car_3165.pgm", "car3Frame3x3FO/car_3165.pgm");
    creation_ppm();

    return 0;
}
