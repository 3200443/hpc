#include <stdio.h>
#include <stdlib.h>
#include "test_mouvement.h"
#include "test_mouvement_SSE2.h"
#define NBIMAGES 199
#include "nrdef.h"
#include "nrutil.h"
#include "matric_roc.h"
#include "test_morpho.h"



void difference2Images()
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
//        sprintf(nomImageLoadNonSIMD,"car3Sigma/car_3%03d.pgm",i);//Image a t
//        sprintf(nomImageLoadSIMD,"car3SigmaSIMD/car_3%03d.pgm",i);//Image a t
//        sprintf(nomImageLoadNormal,"car3/car_3%03d.pgm",i);//Image a t
        sprintf(nomImageLoadNonSIMD,"car3Frame3x3O_pipe/car_3%03d.pgm",i);//Image a t
        sprintf(nomImageLoadSIMD,"car3Frame3x3O_bin/car_3%03d.pgm",i);//Image a t
        sprintf(nomImageLoadNormal,"car3/car_3%03d.pgm",i);//Image a t

        MLoadPGM_ui8matrix(nomImageLoadNonSIMD, nrl, nrh, ncl, nch, ItNonSIMD);
        MLoadPGM_ui8matrix(nomImageLoadSIMD, nrl, nrh, ncl, nch, ItSIMD);
        MLoadPGM_ui8matrix(nomImageLoadNormal, nrl, nrh, ncl, nch, ItNormal);

        for(int k = nrl; k <= nrh; k++)
        {
            for(int l = ncl; l <= nch; l++)
            {
                if(ItNonSIMD[k][l] != ItSIMD[k][l])
                {
                    printf("Probleme ! i = %d, j = %d, ItSIMD = %d, image = %d\n",k, l, ItSIMD[k][l], i);
                }
            }
        }

    }

}

#define OPTI 3 //1 pour optimisation 2 sans optimisation 3 pour tout 0 pour rien


int main(int argc, char* argv[])
{
#if OPTI & 0x1
    test_routine_FrameDifference_SSE2(16);
    test_routine_sigmaDelta_SSE2();
    //test_routine_FrameDifference_SSE2M(16);

#endif
#if OPTI & 0x2
    test_routine_FrameDifference(16);
    test_routine_sigmaDelta();
//
    test_routine_FrameDifferenceMorpho3x3ouverture(16);
    test_routine_FrameDifferenceMorpho3x3fermeture(16);
    test_routine_FrameDifferenceMorpho3x3ouvertureFermeture(16);
    test_routine_FrameDifferenceMorpho3x3fermetureOuverture(16);
    test_routine_FrameDifferenceMorpho3x3fermeturefermeture(16);
    test_routine_FrameDifferenceMorpho3x3ouverture_pipe(16);
    test_routine_FrameDifferenceMorpho3x3fermeture_pipe(16);
    test_routine_FrameDifferenceMorpho3x3ouverture_bin(16);
    test_routine_FrameDifferenceMorpho3x3fermeture_bin(16);


#endif
   // creation_matrices_ROC(argv[1]);
   // difference2Images();
//    test_unitaire_SD_SSE2();

    /*ulong32 test = 1 << 31;
    printf("%lu et %lu", test, test>>1);*/
    return 0;
}
