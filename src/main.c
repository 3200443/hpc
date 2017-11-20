#include <stdio.h>
#include <stdlib.h>
#include "test_mouvement.h"
#include "test_mouvement_SSE2.h"
#define NBIMAGES 199
#include "nrdef.h"
#include "nrutil.h"
#include "matric_roc.h"



void differenceImageScal_SIMD()
{

    long nrl, nrh, ncl, nch;
    uint8 **ItNonSIMD = LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **ItSIMD =  ui8matrix(nrl, nrh, ncl, nch);;
    uint8 **ItNormal = ui8matrix(nrl, nrh, ncl, nch);
    char nomImageLoadNonSIMD[100];//Pour initialiser nrl, etc
    char nomImageLoadSIMD[100];
    char nomImageLoadNormal[100];


    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoadNonSIMD,"car3Sigma/car_3%03d.pgm",i);//Image a t
        sprintf(nomImageLoadSIMD,"car3SigmaSIMD/car_3%03d.pgm",i);//Image a t
        sprintf(nomImageLoadNormal,"car3/car_3%03d.pgm",i);//Image a t

        MLoadPGM_ui8matrix(nomImageLoadNonSIMD, nrl, nrh, ncl, nch, ItNonSIMD);
        MLoadPGM_ui8matrix(nomImageLoadSIMD, nrl, nrh, ncl, nch, ItSIMD);
        MLoadPGM_ui8matrix(nomImageLoadNormal, nrl, nrh, ncl, nch, ItNormal);

        for(int k = nrl; k <= nrh; k++)
        {
            for(int l = ncl; l <= nch; l++)
            {
                if(ItNonSIMD[k][l] != ItSIMD[k][l]){
                    printf("Probleme ! i = %d, j = %d, ItNormal = %d\n",k , l, ItNormal[k][l]);
                }
            }
        }

    }
    
}

#define OPTI 2 //1 pour optimisation 2 sans optimisation 3 pour tout 0 pour rien


int main(int argc, char* argv[])
{
#if OPTI & 0x1
    //test_routine_FrameDifference_SSE2(20);
    //test_routine_sigmaDelta_SSE2();
    test_routine_FrameDifference_SSE2M(20);
#endif
#if OPTI & 0x2
    //test_routine_FrameDifference(20);
    // test_routine_FrameDifferenceMorpho3x3ouverture(20);
    test_routine_FrameDifferenceMorpho3x3fermeture(20);
    //test_routine_FrameDifferenceMorpho3x3ouvertureFermeture(20);
    //test_routine_FrameDifferenceMorpho3x3fermetureOuverture(20);
    //test_routine_FrameDifferenceMorpho3x3fermeturefermeture(20);
    //test_routine_sigmaDelta();
#endif
    //creation_matrices_ROC("verite/car_3165.pgm", "car3Frame3x3FO/car_3165.pgm");
    //differenceImageScal_SIMD();
    //test_unitaire_SD_SSE2();

    return 0;
}
