#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrdef.h"
#include "nrutil.h"
#include "mouvement.h"
#include "morpho.h"
#define NBIMAGES 199
#define BORD 2
//I0 = It et I1 = It-1 : pareil pour tout




void test_routine_sigmaDelta()
{
    //m[nrl..nrh][ncl..nch]
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;


    sprintf(nomImageLoad,"car3/car_3000.pgm");//Image a t-1
    uint8 **Itm1 =  LoadPGM_ui8matrix(nomImageLoad, &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);


    uint8 **Vt = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Vtm1 =  ui8matrix(nrl, nrh, ncl, nch);

    uint8 **Mt = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Mtm1 =  ui8matrix(nrl, nrh, ncl, nch);

    routine_SigmaDelta_step0(Itm1, Mtm1, Vtm1, nrl, nrh, ncl, nch);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad,"car3/car_3%03d.pgm",i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);
        
        routine_SigmaDelta_1step(It, Itm1, Vt, Vtm1, Mt, Mtm1, Et, nrl, nrh, ncl, nch);
        //routine_SigmaDelta_1stepO(It, Itm1, Vt, Vtm1, Mt, Mtm1, Et, nrl, nrh, ncl, nch); // Pas si vraiment optimisé que ca quand on compile avec -O3...

        sprintf(nomImageSave,"car3Sigma/car_3%03d.pgm",i);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);//Copie de t a t-1
        copy_ui8matrix_ui8matrix(Mt, nrl, nrh, ncl, nch, Mtm1);
        copy_ui8matrix_ui8matrix(Vt, nrl, nrh, ncl, nch, Vtm1);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);

        /*memcpy(Mtm1[nrl], Mt[nrl], sizeof(uint8)*(nrow*ncol));
        memcpy(Vtm1[nrl], Vt[nrl], sizeof(uint8)*(nrow*ncol));
        memcpy(Itm1[nrl], It[nrl], sizeof(uint8)*(nrow*ncol));*/

    }
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(It, nrl, nrh, ncl, nch);
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch);
    free_ui8matrix(Vt, nrl, nrh, ncl, nch);
    free_ui8matrix(Vtm1, nrl, nrh, ncl, nch);
    free_ui8matrix(Mt, nrl, nrh, ncl, nch);
    free_ui8matrix(Mtm1, nrl, nrh, ncl, nch);





}

void test_routine_FrameDifference(int seuil)
{
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame/car_3%03d.pgm", i);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
}


void test_routine_FrameDifferenceMorpho3x3ouverture(int seuil)
{
    /*pour commit */
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **Et1 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3O/car_3%03d.pgm", i);
        ouverture3x3(Et,Et1, nrl,nrh,ncl,nch);
        SavePGM_ui8matrix(Et1, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(Et1, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );

}

void test_routine_FrameDifferenceMorpho3x3fermeture(int seuil)
{
    /*pour commit */
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **Et1 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3F/car_3%03d.pgm", i);
        fermeture3x3(Et,Et1, nrl,nrh,ncl,nch);
        SavePGM_ui8matrix(Et1, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(Et1, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );

}
void test_routine_FrameDifferenceMorpho3x3ouvertureFermeture(int seuil)
{
    /*pour commit */
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **Et1 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3OF/car_3%03d.pgm", i);
        ouverture3x3(Et,Et1, nrl,nrh,ncl,nch);
        fermeture3x3(Et1,Et, nrl,nrh,ncl,nch);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(Et1, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );

}

void test_routine_FrameDifferenceMorpho3x3fermetureOuverture(int seuil)
{
    /*pour commit */
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **Et1 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3FO/car_3%03d.pgm", i);
        fermeture3x3(Et,Et1, nrl,nrh,ncl,nch);
        ouverture3x3(Et1,Et, nrl,nrh,ncl,nch);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(Et1, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );

}

void test_routine_FrameDifferenceMorpho3x3fermeturefermeture(int seuil)
{
    /*pour commit */
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);
    uint8 **Et1 = ui8matrix(nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD);

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3FF/car_3%03d.pgm", i);
        ouverture3x3(Et,Et1, nrl,nrh,ncl,nch);
        fermeture3x3(Et1,Et, nrl,nrh,ncl,nch);
        fermeture3x3(Et,Et1, nrl,nrh,ncl,nch);
        fermeture3x3(Et1,Et, nrl,nrh,ncl,nch);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        copy_ui8matrix_ui8matrix(It, nrl, nrh, ncl, nch, Itm1);
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );
    free_ui8matrix(Et1, nrl-BORD, nrh+BORD, ncl-BORD, nch+BORD );

}
