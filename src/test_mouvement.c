#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nrdef.h"
#include "nrutil.h"
#include "mouvement.h"
#include "morpho.h"
#define NBIMAGES 199
//I0 = It et I1 = It-1 : pareil pour tout



void test_unitaire_FD()
{
    long nrl=0, nrh = 1, ncl = 0, nch = 7;
    uint8 It[2][8]={{245, 216, 145, 196, 255, 0, 147, 196},{3, 220, 206, 42, 237, 135, 24, 87}};
    uint8 Itm1[2][8]= {{235, 226, 134, 207, 255, 255, 138, 205},{84, 194, 135, 235, 177, 66, 212, 29}};
    uint8 Et[2][8];

    int nrow=nrh-nrl+1,ncol=nch-ncl+1;


    uint8 Rep[2][8] = {{255, 255, 255, 255, 0, 255, 0, 0 },{255, 255, 255,255, 255,255, 255,255}};
    for(int i = nrl; i <= nrh; i++ )
    {
      for(int j = ncl; j <= nch; j++)
      {
      }
  }

}



void test_routine_sigmaDelta()
{
    //m[nrl..nrh][ncl..nch]
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;


    sprintf(nomImageLoad,"car3/car_3000.pgm");//Image a t-1
    uint8 **Itm1 =  LoadPGM_ui8matrix(nomImageLoad, &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);

    int nrow=nrh-nrl+1,ncol=nch-ncl+1;

    uint8 **Et = ui8matrix(nrl, nrh, ncl, nch);


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
        memcpy(Mtm1[nrl], Mt[nrl], sizeof(uint8)*(nrow*ncol));
        memcpy(Vtm1[nrl], Vt[nrl], sizeof(uint8)*(nrow*ncol));
        memcpy(Itm1[nrl], It[nrl], sizeof(uint8)*(nrow*ncol));

    }
    free_ui8matrix(Et, nrl, nrh, ncl, nch);
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
    uint8 **Et = ui8matrix(nrl, nrh, ncl, nch);

    int nrow=nrh-nrl+1,ncol=nch-ncl+1;

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame/car_3%03d.pgm", i);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        memcpy(Itm1[nrl], It[nrl], sizeof(uint8)*(nrow*ncol));
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl, nrh, ncl, nch );
}


void test_routine_FrameDifferenceMorpho3x3ouverture(int seuil)
{
    /*pour commit */
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et1 = ui8matrix(nrl, nrh, ncl, nch);
    int nrow=nrh-nrl+1,ncol=nch-ncl+1;

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3O/car_3%03d.pgm", i);
        ouverture3x3(Et,Et1, nrl,nrh,ncl,nch);
        SavePGM_ui8matrix(Et1, nrl, nrh, ncl, nch, nomImageSave);
        memcpy(Itm1[nrl], It[nrl], sizeof(uint8)*(nrow*ncol));
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl, nrh, ncl, nch );
    free_ui8matrix(Et1, nrl, nrh, ncl, nch );

}

void test_routine_FrameDifferenceMorpho3x3fermeture(int seuil)
{
    /*pour commit */
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et1 = ui8matrix(nrl, nrh, ncl, nch);
    int nrow=nrh-nrl+1,ncol=nch-ncl+1;

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3F/car_3%03d.pgm", i);
        fermeture3x3(Et,Et1, nrl,nrh,ncl,nch);
        SavePGM_ui8matrix(Et1, nrl, nrh, ncl, nch, nomImageSave);
        memcpy(Itm1[nrl], It[nrl], sizeof(uint8)*(nrow*ncol));
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl, nrh, ncl, nch );
    free_ui8matrix(Et1, nrl, nrh, ncl, nch );

}
void test_routine_FrameDifferenceMorpho3x3ouvertureFermeture(int seuil)
{
    /*pour commit */
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et1 = ui8matrix(nrl, nrh, ncl, nch);
    int nrow=nrh-nrl+1,ncol=nch-ncl+1;

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3OF/car_3%03d.pgm", i);
        ouverture3x3(Et,Et1, nrl,nrh,ncl,nch);
        fermeture3x3(Et1,Et, nrl,nrh-1,ncl,nch);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        memcpy(Itm1[nrl], It[nrl], sizeof(uint8)*(nrow*ncol));
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl, nrh, ncl, nch );
    free_ui8matrix(Et1, nrl, nrh, ncl, nch );

}

void test_routine_FrameDifferenceMorpho3x3fermetureOuverture(int seuil)
{
    /*pour commit */
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et1 = ui8matrix(nrl, nrh, ncl, nch);
    int nrow=nrh-nrl+1,ncol=nch-ncl+1;

    for(int i = 1; i <= NBIMAGES; i++)
    {
        sprintf(nomImageLoad, "car3/car_3%03d.pgm", i);//Image a t
        MLoadPGM_ui8matrix(nomImageLoad, nrl, nrh, ncl, nch, It);

        routine_FrameDifference(It, Itm1, Et, nrl,nrh,ncl,nch, seuil);
        sprintf(nomImageSave, "car3Frame3x3FO/car_3%03d.pgm", i);
        fermeture3x3(Et,Et1, nrl,nrh,ncl,nch);
        ouverture3x3(Et1,Et, nrl,nrh,ncl,nch);
        SavePGM_ui8matrix(Et, nrl, nrh, ncl, nch, nomImageSave);
        memcpy(Itm1[nrl], It[nrl], sizeof(uint8)*(nrow*ncol));
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl, nrh, ncl, nch );
    free_ui8matrix(Et1, nrl, nrh, ncl, nch );

}

void test_routine_FrameDifferenceMorpho3x3fermeturefermeture(int seuil)
{
    /*pour commit */
    char nomImageLoad[50];// = "car3/car_3";
    char nomImageSave[50];// = "car3Sigma/car_3"
    long nrl, nrh, ncl, nch;
    uint8 **Itm1 =  LoadPGM_ui8matrix("car3/car_3000.pgm", &nrl, &nrh, &ncl, &nch);
    uint8 **It = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et = ui8matrix(nrl, nrh, ncl, nch);
    uint8 **Et1 = ui8matrix(nrl, nrh, ncl, nch);
    int nrow=nrh-nrl+1,ncol=nch-ncl+1;

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
        memcpy(Itm1[nrl], It[nrl], sizeof(uint8)*(nrow*ncol));
    }

    free_ui8matrix(It, nrl, nrh, ncl, nch );
    free_ui8matrix(Itm1, nrl, nrh, ncl, nch );
    free_ui8matrix(Et, nrl, nrh, ncl, nch );
    free_ui8matrix(Et1, nrl, nrh, ncl, nch );

}
