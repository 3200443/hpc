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
    uint8 **ItImage1 = LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **ItImage2 =  ui8matrix(nrl, nrh, ncl, nch);;
    uint8 **ItNormal = ui8matrix(nrl, nrh, ncl, nch);
    char nomImage1[100];//Pour initialiser nrl, etc
    char nonImage2[100];
    char nomImageLoadNormal[100];


    for(int i = 1; i <= NBIMAGES; i++)
    {
//        sprintf(nomImage1,"car3Sigma/car_3%03d.pgm",i);//Image a t
//        sprintf(nonImage2,"car3SigmaSIMD/car_3%03d.pgm",i);//Image a t
//        sprintf(nomImageLoadNormal,"car3/car_3%03d.pgm",i);//Image a t
        sprintf(nomImage1,"car3Frame3x3O/car_3%03d.pgm",i);//Image a t
        sprintf(nonImage2,"car3Frame3x3O_pipe/car_3%03d.pgm",i);//Image a t
        sprintf(nomImageLoadNormal,"car3/car_3%03d.pgm",i);//Image a t

        MLoadPGM_ui8matrix(nomImage1, nrl, nrh, ncl, nch, ItImage1);
        MLoadPGM_ui8matrix(nonImage2, nrl, nrh, ncl, nch, ItImage2);
        MLoadPGM_ui8matrix(nomImageLoadNormal, nrl, nrh, ncl, nch, ItNormal);

        for(int k = nrl; k <= nrh; k++)
        {
            for(int l = ncl; l <= nch; l++)
            {
                if(ItImage1[k][l] != ItImage2[k][l])
                {
                    printf("Probleme ! i = %d, j = %d, ItImage2 = %d, image = %d\n",k, l, ItImage2[k][l], i);
                }
            }
        }

    }

}

#define OPTI 3 //1 pour optimisation 2 sans optimisation 3 pour tout 0 pour rien


int main(int argc, char* argv[])
{
#if OPTI & 0x1
    test_routine_FrameDifference_SSE2(20);
    test_routine_sigmaDelta_SSE2();
    test_routine_FrameDifference_SSE2M(20);

#endif
#if OPTI & 0x2
    test_routine_FrameDifference(20);
    test_routine_sigmaDelta();

    test_routine_FrameDifferenceMorpho3x3ouverture(20);
    test_routine_FrameDifferenceMorpho3x3fermeture(20);
    test_routine_FrameDifferenceMorpho3x3ouvertureFermeture(20);
    test_routine_FrameDifferenceMorpho3x3fermetureOuverture(20);
    test_routine_FrameDifferenceMorpho3x3fermeturefermeture(20);
    test_routine_FrameDifferenceMorpho3x3ouverture_pipe(20);
    test_routine_FrameDifferenceMorpho3x3fermeture_pipe(20);
    test_routine_FrameDifferenceMorpho3x3ouverture_bin(20);
    test_routine_FrameDifferenceMorpho3x3fermeture_bin(20);


#endif
    // creation_matrices_ROC(argv[1]);
    //difference2Images();

    return 0;
}
